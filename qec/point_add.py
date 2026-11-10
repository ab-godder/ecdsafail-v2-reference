"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import (Mod, _Gates, _add_selected_tree, _doubling, _emit, _maj, add_const, adjacent_xor, cswap,
                    fan_into, load_const, lookahead_sweep, tree_temps)
from .builder import Circuit
from .curve import Curve
from .square import SquareFold
from .two_bit import TwoBitMod
from .model import tick_bound, tick_widths
from . import block0_code, tail_code, tape_block
from .tape_block import FREE, TICKS


def _fanout(C, src, copies):
    have, todo = [src], list(copies)
    while todo:
        nxt = []
        for h in list(have):
            if not todo:
                break
            t = todo.pop(0)
            C.cx(h, t)
            nxt.append(t)
        have += nxt

def _unfanout(C, src, copies):
    # The doubling tree's CX in reverse order.
    rec = []
    have, todo = [src], list(copies)
    while todo:
        nxt = []
        for h in list(have):
            if not todo:
                break
            t = todo.pop(0)
            rec.append((h, t))
            nxt.append(t)
        have += nxt
    for h, t in reversed(rec):
        C.cx(h, t)

def _par_rotate(C, two, V, cp):
    """Rotate V left by one under `two`: reverse V[0..L-1], then reverse V[0..L-2]."""
    L = len(V)
    # With fewer clean wires than L // 2 - 1, the copies are shared round-robin; the rotation is
    # the same permutation, only less parallel.
    ctl = [two] + cp
    _fanout(C, two, cp)
    for p in range(L // 2):
        cswap(C, ctl[p % len(ctl)], [V[p]], [V[L - 1 - p]])
    for p in range((L - 1) // 2):
        cswap(C, ctl[p % len(ctl)], [V[p]], [V[L - 2 - p]])
    _unfanout(C, two, cp)


RCAP = 6


def ring_layout(n):
    """Split n ring wires: h, r copies each of g0 and g1 (every bit's select reads both, and one
    wire for all bits serializes them), and m >= 2 carry-copy slots."""
    L = n - 1
    r = min(RCAP, max(0, (L - 2) // 3))
    return r, L - 2 * r

def shared_difference_u1(C, U, V, g0, g1, ring):
    w = len(U)
    assert w >= 2 and len(ring) >= 3
    r, m = ring_layout(len(ring))
    h = ring[0]
    G0, G1, slots = list(ring[1:1 + r]), list(ring[1 + r:1 + 2 * r]), list(ring[1 + 2 * r:])
    assert len(slots) == m >= 2
    (cin,) = C.alloc(1)
    carries = [cin] + U[:-1]
    for wire in V:
        C.x(wire)
    for c in G0:
        C.cx(g0, c)
    for i in range(w):
        if i == 0:
            C.cx(U[i], V[i])
            C.cx(U[i], carries[i])
            C.cx(V[i], U[i])
        else:
            _maj(C, carries[i], V[i], U[i])
    C.ccx(g0, U[-1], g1)
    for c in G1:
        C.cx(g1, c)
    C.cx(g0, h)
    C.cx(g1, h)
    C.x(h)
    K = {}
    a0 = lambda i: G0[i % r] if r else g0
    a1 = lambda i: G1[i % r] if r else g1

    def select(i):
        u, v, k = U[i], V[i], K[i]
        C.cx(u, v)
        C.ccx(a0(i), k, v)
        C.cx(u, k)
        C.cx(k, v)
        C.ccx(a1(i), v, u)
        C.cx(k, v)

    def finish(i):
        # unc(i) clears K[i + 1]; then the bit's final values
        u, v = U[i], V[i]
        if i + 1 in K:
            C.cx(a0(i), v)
            C.cx(u, v)
            C.ccx(K[i], v, K[i + 1])
            C.cx(u, v)
            C.ccx(u, v, K[i + 1])
            C.cx(h, v)
            del K[i + 1]
        else:
            C.x(v)
            C.cx(a1(i), v)
        C.cx(a1(i), u)

    for i in reversed(range(1, w)):
        K[i] = slots[i % m]
        C.cx(U[i - 1], K[i])
        C.ccx(U[i - 1], V[i], U[i])
        if i != w - 1:
            C.cx(U[i + 1], U[i])
            select(i + 1)
            finish(i + 1)
    C.cx(U[1], U[0])
    select(1)
    finish(1)
    C.cx(U[0], K[1])
    del K[1]
    assert not K
    # bit 0 as the tree with R
    C.cx(V[0], U[0])
    C.cx(g0, V[0])
    C.cx(carries[0], V[0])
    C.cx(U[0], carries[0])
    C.cx(carries[0], V[0])
    C.free([cin])
    C.x(V[0])
    C.x(h)
    C.cx(g1, h)
    C.cx(g0, h)
    for c in G1:
        C.cx(g1, c)
    for c in G0:
        C.cx(g0, c)


U6_RING = 193
# Give a step the bit-level tree's temps when that is more than U6_RING and the pass has them.
WIDE_RING = True
# Run the first Euclid step's control fans, constant add and rotation in log depth, on TICK0_SCRATCH
# wires of unwritten tape.
TICK0_TREE = True
TICK0_SCRATCH = 450
TICK0_COPIES = 31
# Clean wires a step may take beyond the pass's tape for the bit-level tree's ring. The steps short
# of a full ring take 42 at most, and the pass runs further than that below the circuit's peak, so
# they cost no wire.
RING_BORROW = 300
# The u6 step's copies of g0 and of g1 on wires of their own, this many each, instead of on the
# ring: the copies then wait for no tree temp and no tree temp waits for them. The Euclid passes run
# about 274 wires below the circuit's peak; these and the TWO_EARLY wires take 274.
U6_OWN_COPIES = 72
# Compute a step's second-strip flag `two` at the step's start instead of after its shared
# difference, and fan it onto pass-level copies while the shared difference runs. two = 1 implies
# g0 = 0, and the shared difference is the identity when g0 = 0, so `two` read from V[1] before it
# equals `two` read from V[0] after the strip. The pass alternates two wires for `two`, so a step's
# flag does not wait for the previous step's uncompute; the copies are shared.
TWO_EARLY = True
TWO_COPIES = 200


def _ring_size(w, avail):
    """The clean tape wires a step at width w takes as its ring, from avail."""
    k = min(max(U6_RING, (w - 1) // 2 - 1), avail)
    full = tree_temps(w + 1, 1)
    return full if WIDE_RING and k < full <= avail else k


def ring_layout_u6(n):
    """Copies of g0 and of g1, n // 2 each, and lookahead temps, all n, on the same ring wires at
    different times. None below 7."""
    if n < 7:
        return None
    return n // 2, n


def _fan(C, src, copies):
    """CX src into each clean copy by doubling: depth about log2(len(copies) + 1)."""
    have = [src]
    todo = list(copies)
    while todo:
        new = []
        for s in have:
            if not todo:
                break
            d = todo.pop(0)
            C.cx(s, d)
            new.append(d)
        have += new


def _unfan(C, src, copies):
    g = _Gates()
    _fan(g, src, copies)
    _emit(C, reversed(g.ops))


def shared_difference_u6(C, U, V, g0, g1, ring, peek=None, off=None):
    """The step's shared difference with its down-sweep as the up-sweep run backwards.
    On the step's domain (U odd, so a_0 = 1; g0 = V[0]; g1 = 0), with b' = ~b and s = a + b'
    mod 2^w, the tree's unwind leaves (a, b) if g0 = 0, (a, ~s) if g0 and not g1, and (~b', s)
    if g1 = g0 c_w, except V[0] = ~s_0 whatever g1 is. The carries of a + ~s and of b' + ~s equal
    those of a + b', so after the up-sweep the up-sweep states of (a, ~s) and (b', ~s) are a CX or
    two a bit away (B_i = ~(p_i ^ k_i); A_(i-1) = k_i ^ p_i and B_i = ~k_i), and the same sweep run
    backwards leaves them. No carry copies and no uncopy chain. The whole ring is the sweeps'
    temps, which they leave clean; copies of g0 and g1 sit on it only for the middle step and
    for the select."""
    w = len(U)
    lay = ring_layout_u6(len(ring))
    assert w >= 2 and lay is not None
    r, _ = lay
    if U6_OWN_COPIES:
        r = min(w, U6_OWN_COPIES)
        own = C.alloc(2 * r)
        G0, G1 = own[:r], own[r:]
    else:
        G0, G1 = list(ring[:r]), list(ring[r:2 * r])
    (cin,) = C.alloc(1)
    for wire in V:
        C.x(wire)
    up = _Gates()
    for i in range(w):
        up.cx(U[i], V[i])
    up.cx(U[0], cin)
    adjacent_xor(up, U)
    lookahead_sweep(up, [cin] + list(U), [None] + list(V), list(ring))
    if off is None:
        _emit(C, up.ops)
    else:
        # V[0] takes the held offset here rather than before: the XOR commutes with the X and the
        # CX on V[0] ahead of it, and no gate reads V[0] in between.
        _emit(C, up.ops[:w])
        C.cx(off, V[0])
        _emit(C, up.ops[w:])
    if peek is not None:
        # peek holds ~(b_1 ^ off), with b_1 V[1]'s input. Now V[1] = p_1 = a_1 ^ ~b_1 and
        # U[0] = c_1 ^ a_1 with c_1 = a_0 b'_0 = ~g0 (a_0 = 1), so these XORs clear it.
        C.cx(V[1], peek)
        C.cx(U[0], peek)
        C.cx(g0, peek)
        C.cx(off, peek)
        C.x(peek)
    _fan(C, g0, G0)
    C.ccx(g0, U[-1], g1)
    _fan(C, g1, G1)
    for i in range(w):
        k = cin if i == 0 else U[i - 1]
        C.ccx(G1[i % r], V[i], k)
        C.ccx(G0[i % r], k, V[i])
        C.cx(G0[i % r], V[i])
    _unfan(C, g1, G1)
    _unfan(C, g0, G0)
    _emit(C, reversed(up.ops))
    _fan(C, g0, G0)
    _fan(C, g1, G1)
    for i in range(w):
        C.cx(G1[i % r], U[i])
        C.x(V[i])
        C.cx(G0[i % r], V[i])
        if i:
            C.cx(G1[i % r], V[i])
    _unfan(C, g1, G1)
    _unfan(C, g0, G0)
    C.free([cin])
    if U6_OWN_COPIES:
        C.free(own)


# Widths at which the step uses the in-place ripple below instead of a ring variant: there it has
# fewer ticks and fewer CCX than u1 and u6 (measured in isolation, w = 2..12).
RIPPLE_MAX_W = 12


def shared_difference_ripple(C, U, V, g0, g1, ring):
    """The step's shared difference as one in-place majority chain, from 2elliot depth/v2 (9d34aef);
    the same map as the chain fallback on every input.

    Up: each MAJ's two CXs run ahead of every CCX (they commute with the CCXs between), so each
    CCX's first control is the previous CCX's target. Down: the selection moves into (c, v) before
    the carry arrives, by maps that keep c.~v, which is all the restore and the carry read:
    v ^= g0.~c, c ^= v while the up chain runs, then v ^= h1.~c with h1 = g1 ^ g0.carry_out.
    After the restore the tail is linear: u ^= v, v ^= u ^ g0 ^ h1, u ^= g0.h1. k and m copy g0
    and h1 so that the two selections alternate controls, then hold g0.h1 and g0 ^ h1.
    cin, k and m are the first clean ring wires, or fresh wires when the ring is shorter; all three
    end at zero."""
    w = len(U)
    lent = list(ring[:3])
    fresh = C.alloc(3 - len(lent))
    cin, k, m = lent + fresh
    carries = [cin] + U[:-1]
    for wire in V:
        C.x(wire)
    C.cx(g0, k)
    C.cx(U[0], V[0])
    C.cx(U[0], cin)
    for i in range(1, w):
        C.cx(U[i], V[i])
        C.cx(U[i], carries[i])
    C.cx(V[0], U[0])
    for i in range(1, w):
        C.ccx(carries[i], V[i], U[i])
        C.x(carries[i])
        C.ccx(k if (w - 1 - i) % 2 == 0 else g0, carries[i], V[i])
        C.cx(V[i], carries[i])
    C.ccx(g0, U[-1], g1)
    C.cx(g1, m)
    for i in reversed(range(1, w)):
        C.ccx(m if (w - 1 - i) % 2 == 0 else g1, carries[i], V[i])
        C.x(carries[i])
        C.ccx(carries[i], V[i], U[i])
    C.cx(g0, k)
    C.ccx(g0, g1, k)
    C.cx(g0, m)
    for i in reversed(range(1, w)):
        C.cx(U[i], carries[i])
        C.cx(V[i], U[i])
        C.cx(U[i], V[i])
        C.cx(m, V[i])
        C.cx(k, U[i])
    C.cx(V[0], U[0])
    C.cx(g0, V[0])
    C.cx(U[0], cin)
    C.cx(U[0], V[0])
    C.cx(g1, m)
    C.cx(g0, m)
    C.ccx(g0, g1, k)
    C.free(fresh)
    for wire in V:
        C.x(wire)


# Hand each Bezout tick's B1 mark to the next tick when nothing runs between them.
HAND_OFF_MARK = True

# Run the Bezout pass's block codec calls on the free wires freed longest ago.
CODEC_ASIDE = True
# The Euclid pass's back-to-back block encodes take their clean wires from this many sets in turn.
ENCODE_SETS = 3


def _aside(C, fn):
    """Run fn on the free wires freed longest ago, and leave the wires it frees under the rest.
    A codec call between two ticks then takes no wire the tick before it just freed, which it
    would wait on, and the tick after it takes the codec's wires last."""
    if not CODEC_ASIDE:
        return fn()
    rest = C.free_stack
    C.free_stack = rest[::-1]
    out = fn()
    st = C.free_stack
    m = 0
    while m < min(len(st), len(rest)) and st[m] == rest[-1 - m]:
        m += 1
    C.free_stack = st[m:] + rest[len(rest) - m:]
    return out


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for canonical 1 <= x < q; see model.py."""

    def __init__(self, C, mod):
        self.C, self.mod = C, mod
        q, n = mod.q, mod.n
        self.iters = tick_bound(q)
        self.widths = tick_widths(q, self.iters)
        self._euclid_steps = {}
        self._euclid_tick0_block = None
        self._encode_block = None
        self._encode_block0 = None
        self._bezout_step = None
        self._bezout_ticks = {}
        self.narrow = TwoBitMod(q)
        self.block = None

    @property
    def pass_tail_width(self):
        return self.widths[self.iters - 1]

    @property
    def tape_external(self):
        return 2 * self.iters - (self.mod.n - self.pass_tail_width) - self.u_lent

    @property
    def blocks(self):
        return self.iters // TICKS

    @property
    def u_lent(self):
        """U's wires that late tape cells take once they are zero, one per zero wire the pass
        ends with: V's tail and each coded block's free wires, which U takes back."""
        coded_free = sum(
            len(self._block_free(j, self.swept)) for j in range(self.blocks))
        return min(self.pass_tail_width + coded_free, self.mod.n - self.pass_tail_width)

    @property
    def lent_after_pass(self):
        coded_free = sum(
            len(self._block_free(j, self.swept)) for j in range(self.blocks))
        return self.pass_tail_width + coded_free - self.u_lent

    @property
    def swept(self):
        """Whether the tape at rest is the B1-swept one, which the code reads last cell first."""
        return self.narrow.merged

    def with_free_slots(self, wires):
        """A coded tape from its wires, with None where each block's free wires were."""
        G = list(wires)
        slots = sorted(slot for j in range(self.blocks) for slot in self._free_slots(j, self.swept))
        for slot in slots:
            G.insert(slot, None)
        return G

    def coded_tape(self, pool):
        """The tape after a pass."""
        return self.with_free_slots(pool[self.lent_after_pass:])

    def _block_free(self, j, swept):
        return block0_code.FREE if swept and j == 0 else FREE

    def _free_slots(self, j, swept):
        return [2 * ((j + 1) * TICKS - 1 - f // 2 if swept else j * TICKS + f // 2)
                + f % 2 for f in self._block_free(j, swept)]

    def _block_wires(self, G, j, swept):
        cells = range(j * TICKS, (j + 1) * TICKS)
        return [G[2 * k + h] for k in (reversed(cells) if swept else cells) for h in (0, 1)]

    def _encoder(self, j, swept):
        if swept and j == 0:
            if self._encode_block0 is None:
                self._encode_block0 = self.C.record_shape(
                    lambda W: block0_code.encode(self.C, W), [2 * TICKS])
            return self._encode_block0
        if self._encode_block is None:
            self._encode_block = self.C.record_shape(lambda W: tape_block.encode(self.C, W), [2 * TICKS])
        return self._encode_block

    def _encode(self, G, j, swept):
        """Code block j. Its free wires leave G, which keeps None in those places."""
        W = self._block_wires(G, j, swept)
        self.C.apply(self._encoder(j, swept), [W])
        wires = [W[f] for f in self._block_free(j, swept)]
        for slot, wire in zip(self._free_slots(j, swept), wires):
            assert G[slot] == wire
            G[slot] = None
        return wires

    def _decode(self, G, j, wires, swept):
        assert len(wires) == len(self._block_free(j, swept))
        for slot, wire in zip(self._free_slots(j, swept), wires):
            G[slot] = wire
        self.C.apply(self._encoder(j, swept),
                     [self._block_wires(G, j, swept)], inverse=True)

    def _shared_difference(self, U, V, g0, g1, ring_wires=(), peek=None, off=None):
        if peek is not None:
            shared_difference_u6(self.C, U, V, g0, g1, list(ring_wires), peek, off)
            return
        if len(U) <= RIPPLE_MAX_W:
            shared_difference_ripple(self.C, U, V, g0, g1, list(ring_wires))
            return
        if ring_layout_u6(len(ring_wires)) is not None:
            shared_difference_u6(self.C, U, V, g0, g1, list(ring_wires))
            return
        if len(ring_wires) >= 3:
            shared_difference_u1(self.C, U, V, g0, g1, list(ring_wires[:17]))
            return
        """Compare and form the even difference using one majority chain."""
        C = self.C
        (cin,) = C.alloc(1)
        carries = [cin] + U[:-1]
        for wire in V:
            C.x(wire)
        for i in range(len(U)):
            if i == 0:
                C.cx(U[i], V[i])
                C.cx(U[i], carries[i])
                C.cx(V[i], U[i])
            else:
                _maj(C, carries[i], V[i], U[i])
        C.ccx(g0, U[-1], g1)
        for i in reversed(range(len(U))):
            # Restore a_i and select S_i while unwinding the carry chain.
            if i == 0:
                C.cx(V[i], U[i])
                C.cx(g0, V[i])
            else:
                C.ccx(carries[i], V[i], U[i])
                C.ccx(g0, carries[i], V[i])
            C.cx(carries[i], V[i])
            C.cx(U[i], carries[i])
            if i != 0:
                C.ccx(g1, V[i], U[i])
                C.cx(g1, U[i])
            C.cx(carries[i], V[i])
            if i != 0:
                C.cx(g1, V[i])
        C.free([cin])
        for wire in V:
            C.x(wire)

    def _euclid_step(self, w: int, k: int = 0):
        if k == 0 and w in self._euclid_steps:
            return self._euclid_steps[w]
        if (w, k) not in self._euclid_steps:
            C = self.C

            def step_early(U, V, E, cell, cp, TT):
                # The step below with TWO_EARLY: TT is `two` and its copies, all clean.
                g0, b = cell
                offset = E[0]
                two, tw = TT[0], list(TT[1:])
                # g0 = V[0] ^ offset; the shared difference XORs the offset into V[0] itself. V[1] is
                # read through a copy holding ~(b_1 ^ offset), which the shared difference clears,
                # so neither waits for the CCX below or for the offset.
                (peek,) = C.alloc(1)
                C.cx(V[0], g0)
                C.cx(V[1], peek)
                C.cx(offset, g0)
                C.cx(offset, peek)
                C.x(peek)
                C.x(g0)
                C.ccx(g0, peek, two)
                C.x(g0)
                L = w - 1
                _fanout(C, two, tw)
                self._shared_difference(U[:w], V[:w], g0, b, list(cp), peek, offset)
                C.free([peek])
                V[:] = V[1:] + V[:1]
                cswap(C, g0, [offset], [b])
                ctl = [two] + tw
                Vr = V[:L]
                for p in range(L // 2):
                    cswap(C, ctl[p % len(ctl)], [Vr[p]], [Vr[L - 1 - p]])
                for p in range((L - 1) // 2):
                    cswap(C, ctl[p % len(ctl)], [Vr[p]], [Vr[L - 2 - p]])
                # The stripped bit, offset when two = 1, now sits at V[L - 1].
                C.ccx(two, offset, V[L - 1])
                C.cx(two, b)
                _unfanout(C, two, tw)
                C.x(g0)
                C.ccx(g0, b, two)
                C.x(g0)

            if k and isinstance(k, tuple):
                self._euclid_steps[w, k] = C.record_shape(step_early, [w, w, 1, 2, k[0], k[1]])
                return self._euclid_steps[w, k]

            def step(U, V, E, cell, cp=()):
                g0, b = cell
                offset = E[0]
                C.cx(offset, V[0])
                C.cx(V[0], g0)
                self._shared_difference(U[:w], V[:w], g0, b, list(cp))
                V[:] = V[1:] + V[:1]
                # A no-subtract tick strips a second bit when V + E is even.
                (two,) = C.alloc(1)
                C.cx(offset, V[0])
                C.x(g0)
                C.x(V[0])
                C.ccx(g0, V[0], two)
                C.x(V[0])
                C.x(g0)
                C.cx(offset, V[0])
                C.cx(two, b)
                C.ccx(two, offset, V[0])
                # A subtract keeps its swap on E and moves the old offset to the cell. It follows
                # the gates above so that `two` does not wait for it: on g0 they do nothing (two
                # is 0 and the two CXs from offset cancel), and off g0 the swap does nothing.
                cswap(C, g0, [offset], [b])
                if cp:
                    _par_rotate(C, two, V[:w - 1], list(cp))
                else:
                    for i in range(w - 2):
                        cswap(C, two, [V[i]], [V[i + 1]])
                C.x(g0)
                C.ccx(g0, b, two)
                C.x(g0)
                C.free([two])

            if k:
                self._euclid_steps[w, k] = C.record_shape(step, [w, w, 1, 2, k])
            else:
                self._euclid_steps[w] = C.record_shape(step, [w, w, 1, 2])
        return self._euclid_steps[w] if k == 0 else self._euclid_steps[w, k]

    def _euclid_tick0(self):
        """First Euclid step on U=q, E=0, using unwritten tape for carries."""
        if self._euclid_tick0_block is None:
            C, q, n = self.C, self.mod.q, self.mod.n
            f1 = (1 << n) - q + 1
            k = self._tick0_scratch()

            def step(U, V, E, cell, scr):
                C.free(scr)
                g0, b = cell
                C.cx(V[0], g0)
                cps = C.alloc(TICK0_COPIES) if TICK0_TREE else []
                fan = _Gates()
                _doubling(fan, g0, cps)
                _emit(C, fan.ops)
                src = [g0] + cps
                for i in range(1, n):
                    if q >> i & 1:
                        C.x(V[i])
                    C.ccx(src[i % len(src)], V[i], U[i])
                    if q >> i & 1:
                        C.x(V[i])
                _emit(C, reversed(fan.ops))
                if TICK0_TREE:
                    _add_selected_tree(C, V, [g0], [f1], max_temps=TICK0_SCRATCH - n - TICK0_COPIES - 1)
                    fan_into(C, g0, V, cps)
                else:
                    add_const(C, V, f1, ctrl=g0)
                    for wire in V:
                        C.cx(g0, wire)
                C.free(cps)
                C.cx(g0, E[0])
                V[:] = V[1:] + V[:1]
                C.x(g0)
                C.x(V[0])
                C.ccx(g0, V[0], b)
                C.x(V[0])
                C.x(g0)
                if TICK0_TREE:
                    cp = C.alloc((n - 1) // 2 - 1)
                    _par_rotate(C, b, V[:n - 1], cp)
                    C.free(cp)
                else:
                    for i in range(n - 2):
                        cswap(C, b, [V[i]], [V[i + 1]])

            self._euclid_tick0_block = C.record_shape(step, [n, n, 1, 2, k])
        return self._euclid_tick0_block

    def _tick0_scratch(self):
        """The unwritten tape the first step takes: its constant add's carries, or with TICK0_TREE
        the trees' wires."""
        n = self.mod.n
        f1 = (1 << n) - self.mod.q + 1
        return TICK0_SCRATCH if TICK0_TREE else n - 1 - ((f1 & -f1).bit_length() - 1)

    def _euclid_pass(self, U: list[int], pool: list[int]) -> None:
        # The pool starts as V plus clean external tape. It ends as the zero wires U did not
        # take back, then tape; U ends as its low wires and those zero wires. The boundary is a
        # permutation of U and the pool together.
        n = self.mod.n
        V, available, G = pool[:n], pool[n:], []
        lend = U[n - self.u_lent:]
        E = self.C.alloc(1)
        k = self._tick0_scratch()
        if TWO_EARLY:
            two_wires = self.C.alloc(2)
            tw_wires = self.C.alloc(min(TWO_COPIES, (n - 1) // 2 - 1))
        for i in range(self.iters):
            w = self.widths[i]
            assert len(V) == w and len(available) >= 2
            cell, available = available[:2], available[2:]
            tick0 = self._euclid_tick0() if i == 0 and w == n and len(available) >= k else None
            if tick0 is not None and tick0.wires <= set(tick0.flat_in) and set(tick0.flat_out) == set(tick0.flat_in):
                self.C.apply(tick0, [U, V, E, cell, available[:k]])
            else:
                # Alternate clean scratch so adjacent carry and strip sweeps can overlap.
                hold = self.C.alloc(1) if i % 2 else []
                k = _ring_size(w, len(available))
                full = tree_temps(w + 1, 1)
                borrow = (self.C.alloc(full - len(available))
                          if WIDE_RING and k < full <= len(available) + RING_BORROW else [])
                ring = available + borrow if borrow else available[:k]
                if TWO_EARLY and w > RIPPLE_MAX_W and ring_layout_u6(len(ring)) is not None:
                    ntw = min(TWO_COPIES, (w - 1) // 2 - 1)
                    TT = [two_wires[i % 2]] + tw_wires[:ntw]
                    self.C.apply(self._euclid_step(w, (len(ring), len(TT))), [U[:w], V, E, cell, ring, TT])
                else:
                    self.C.apply(self._euclid_step(w, len(ring)), [U[:w], V, E, cell, ring])
                self.C.free(borrow)
                self.C.free(hold)
            G.extend(cell)
            next_width = self.widths[i + 1] if i + 1 < self.iters else self.pass_tail_width
            assert w - 1 <= next_width <= w
            available.extend(V[next_width:])
            V[:] = V[:next_width]
            # By L6, u < 2**next_width from here on, and later steps touch only U[:w] for
            # w <= next_width, so U's wires from next_width up are zero for the rest of the pass.
            available.extend(reversed([u for u in U[next_width:w] if u in lend]))
        assert len(V) == self.pass_tail_width and not available
        if TWO_EARLY:
            self.C.free(two_wires + tw_wires)
        # Sweep down, coding each block as it finishes. The offset sweep leaves b on a
        # no-subtract tick alone. On the merged path the B1 sweep follows one cell behind it,
        # so a block is coded only after the B1 step that reads its bottom cell.
        T, swept, clean = self.iters, self.swept, []
        both = None

        def mark(k):
            self.C.x(G[2 * k])
            self.C.x(G[2 * k - 2])
            self.C.ccx(G[2 * k], G[2 * k - 2], both)
            self.C.x(G[2 * k - 2])
            self.C.x(G[2 * k])

        if swept and T >= 2:
            (both,) = self.C.alloc(1)
            mark(T - 1)
        for i in range(T - 1, -1, -1):
            cswap(self.C, G[2 * i], E, [G[2 * i + 1]])
            if swept and i + 1 < T:
                k = i + 1
                self._b1_step(G, k, both)
                if k >= 2:
                    # Turn mark(k,k-1) into mark(k-1,k-2) before any encode.
                    self.C.cx(G[2 * k - 4], G[2 * k])
                    self.C.x(G[2 * k - 2])
                    self.C.ccx(G[2 * k - 2], G[2 * k], both)
                    self.C.x(G[2 * k - 2])
                    self.C.cx(G[2 * k - 4], G[2 * k])
                else:
                    mark(1)
            done = i + 1 if swept else i
            if done % TICKS == 0 and done // TICKS < self.blocks:
                # Rotate the encoder's clean wires, so that a block's encode overlaps the next ones.
                hold = self.C.alloc(tape_block.ANCILLAS * (done // TICKS % ENCODE_SETS))
                clean.extend(self._encode(G, done // TICKS, swept))
                self.C.free(hold)
        if both is not None:
            self.C.free([both])
        if swept and self.blocks:
            clean.extend(self._encode(G, 0, swept))
        self.C.free(E)
        zeros = V + clean
        U[:] = U[:n - self.u_lent] + zeros[:self.u_lent]
        pool[:] = zeros[self.u_lent:] + [w for w in G if w is not None]

    def _bezout(self):
        if self._bezout_step is None:
            C, mod = self.C, self.mod

            def step(g, R, S):
                self.narrow.step(C, g, R, S)

            self._bezout_step = C.record_shape(step, [2, mod.n, mod.n])
        return self._bezout_step

    def _bezout_tick(self, first, last, pebble=False, consume=False, produce=False):
        """consume: the previous cell's B1 mark arrives on a wire from the tick before.
        produce: this tick leaves the next tick's mark, NOR of its own cell, on a given clean wire."""
        key = first, last, pebble, consume, produce
        if key not in self._bezout_ticks:
            C, n = self.C, self.mod.n

            def tick(prev, g, nxt, R, S, pin, pout):
                self.narrow.tick(C, prev or None, g, nxt[0] if nxt else None,
                                 R, S, pebble=pebble, pin=pin[0] if pin else None,
                                 pout=pout[0] if pout else None)

            sizes = [0 if first else 2, 2, 0 if last else 1, n, n, int(consume), int(produce)]
            self._bezout_ticks[key] = C.record_shape(tick, sizes)
        return self._bezout_ticks[key]

    def _b1_step(self, G, i, both):
        """Use the pass-wide mark while both cells are still raw."""
        cswap(self.C, both, [G[2 * i + 1]], [G[2 * i - 1]])

    def _bezout_pass(self, G, Yr, S):
        """G is the coded tape. Each block is decoded on a clean wire before a tick reads it and
        coded again once no tick reads it."""
        C, T, swept = self.C, self.iters, self.swept
        if not self.narrow.merged:
            bz = self._bezout()
            for i in reversed(range(T)):
                j = i // TICKS
                if j < self.blocks and i % TICKS == TICKS - 1:
                    self._decode(G, j, C.alloc(len(self._block_free(j, swept))), swept)
                C.apply(bz, [G[2 * i:2 * i + 2], Yr, S])
                if j < self.blocks and i % TICKS == 0:
                    C.free(self._encode(G, j, swept))
            return
        # A tick reads the previous cell and writes the next cell's g0, so a block is decoded
        # one tick before its first tick and coded one tick after its last. Across each edge
        # two blocks are decoded for two ticks.
        if self.blocks * TICKS == T:
            _aside(C, lambda: self._decode(G, self.blocks - 1,
                                           C.alloc(len(self._block_free(self.blocks - 1, swept))), swept))
            if self.blocks == 1:
                block0_code.hold(C, self._block_wires(G, 0, swept))
                C.free([G[1]])
                G[1] = None
        # The final swept cell is (1, 0) or (0, 1): the terminal subtract
        # is followed only by no-subtract cells. Lend its b only when uncoded.
        lend_tail = T > 1 and self.blocks * TICKS < T
        tail = (lend_tail and T >= tail_code.CELLS + 1
                and T - self.blocks * TICKS >= tail_code.CELLS
                and tick_widths(self.mod.q, T)[T - tail_code.CELLS] <= 6)
        lo = 2 * (T - tail_code.CELLS)
        tail_encoded = False

        def quiet(i):
            """Nothing is emitted between tick call i and tick call i - 1, which has a prev."""
            j = i - 1
            if not HAND_OFF_MARK or i == T - 1 or j < 0:
                return False
            enc = (i + 1) % TICKS == 0 and i + 1 < T and (i + 1) // TICKS < self.blocks
            tl = tail and i == T - tail_code.CELLS - 1
            dec = j % TICKS == 0 and 0 <= j // TICKS - 1 < self.blocks
            hold = j == 1 and self.blocks
            np = lend_tail and j == T - 2
            return not (enc or tl or dec or hold or np)

        hand, mark_at = None, None
        for i in reversed(range(T)):
            nxt_block = i // TICKS - 1
            if i % TICKS == 0 and 0 <= nxt_block < self.blocks:
                _aside(C, lambda: self._decode(G, nxt_block,
                                               C.alloc(len(self._block_free(nxt_block, swept))), swept))
                if nxt_block == 0:
                    block0_code.hold(C, self._block_wires(G, 0, swept))
                    C.free([G[1]])
                    G[1] = None
            if i == 1 and self.blocks:
                G[1] = C.alloc(1)[0]
                block0_code.hold(C, self._block_wires(G, 0, swept), inverse=True)
            no_prev = lend_tail and i == T - 2
            if no_prev:
                # The terminal cell is never B1; release b before its last reader.
                C.cx(G[-2], G[-1])
                C.x(G[-1])
                C.free([G[-1]])
                G[-1] = None
            first, last = i == T - 1 or no_prev, i == 0
            edge = any(i in (j * TICKS, j * TICKS - 1) for j in range(1, self.blocks))
            pebble = edge or (tail and self.blocks > 0 and i == self.blocks * TICKS)
            pebble = pebble or (self.blocks > 0 and i in (0, 1))
            prev = [] if first else G[2 * i + 2:2 * i + 4]
            nxt = [] if last else G[2 * i - 2:2 * i - 1]
            if i == T - 1:
                # S starts at zero and the terminal cell is A0, A1 or J, never B1.
                for r, s in zip(Yr, S):
                    C.ccx(G[2 * i], r, s)
            else:
                consume, produce = hand is not None, quiet(i)
                if consume:
                    # The mark is valid only if nothing ran since the tick that made it.
                    assert mark_at == (len(C.prims), len(C.calls or ())) and not first
                pin = [hand] if consume else []
                pout = C.alloc(1) if produce else []
                C.apply(self._bezout_tick(first, last, pebble, consume, produce),
                        [prev, G[2 * i:2 * i + 2], nxt, Yr, S, pin, pout])
                C.free(pin)
                hand = pout[0] if produce else None
                mark_at = (len(C.prims), len(C.calls or ()))
            prev_block = (i + 1) // TICKS
            if (i + 1) % TICKS == 0 and i + 1 < T and prev_block < self.blocks:
                _aside(C, lambda: C.free(self._encode(G, prev_block, swept)))
            if tail and i == T - tail_code.CELLS - 1:
                W = G[lo:lo + tail_code.LIVE]
                tail_code.apply(C, W)
                tail_encoded = True
                C.free([W[tail_code.FREE[0]]])
                G[lo + tail_code.FREE[0]] = None
        if self.blocks:
            _aside(C, lambda: C.free(self._encode(G, 0, swept)))
        if tail:
            assert tail_encoded, "tail inverse requires its forward apply"
            G[lo + tail_code.FREE[0]] = C.alloc(1)[0]
            tail_code.apply(C, G[lo:lo + tail_code.LIVE], inverse=True)
        if lend_tail:
            G[-1] = C.alloc(1)[0]
            C.x(G[-1])
            C.cx(G[-2], G[-1])

    def _body(self, Xr, Yr):
        C, q, n, T = self.C, self.mod.q, self.mod.n, self.iters
        U = C.alloc(n)
        load_const(C, U, q)
        pool = list(Xr) + C.alloc(self.tape_external)

        E = C.run(self._euclid_pass, [U, pool])
        lent = self.lent_after_pass
        G = self.coded_tape(pool)
        # Terminal V and the wires the code freed are zero; lend them during Bezout.
        C.free(pool[:lent])
        # U is (1,0,...,0). Reuse its entire word for the reverse coefficient;
        # the surviving V wire stays zero for the inverse Euclid boundary.
        C.x(U[0])
        S = U
        self._bezout_pass(G, Yr, S)
        # Restore the tape boundary using the newly allocated tail wire.
        pool[lent:] = [w for w in G if w is not None]
        C.x(Yr[0])
        pool[:lent] = C.alloc(lent)
        C.apply(E, [Yr, pool], inverse=True)
        load_const(C, Yr, q)
        Xr[:] = pool[:n]
        C.free(Yr + pool[n:])
        Yr[:] = S

    def apply(self, Xr, Yr, inverse=False):
        if self.block is None:
            self.block = self.C.record(self._body, [Xr, Yr])
        self.C.apply(self.block, [Xr, Yr], inverse)


# AND the flags' bits by balanced trees instead of chains: the same CCX, log depth.
EQ_TREE = True


def and_tree(C, wires):
    """The AND of wires (at least one) on fresh wires by a balanced tree. Returns the root, which is
    wires[0] itself for one wire, and the CCXs (a, b, t) in emission order, for clear_and."""
    level, history = list(wires), []
    while len(level) > 1:
        nxt = []
        for k in range(0, len(level) - 1, 2):
            (t,) = C.alloc(1)
            C.ccx(level[k], level[k + 1], t)
            history.append((level[k], level[k + 1], t))
            nxt.append(t)
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    return level[0], history


def clear_and(C, history):
    for a, b, t in reversed(history):
        C.ccx(a, b, t)
    C.free([t for _, _, t in history])


def eq_const_into(C, ctrl, bits, value, target):
    """XOR equality into target, optionally controlled by ctrl."""
    flips = [w for i, w in enumerate(bits) if not value >> i & 1]
    for w in flips:
        C.x(w)
    if EQ_TREE and (ctrl is not None or len(bits) > 1):
        # ctrl joins last: callers pass a ctrl that is itself a tree's root, ready late.
        if ctrl is not None:
            root, history = and_tree(C, bits)
            C.ccx(ctrl, root, target)
            clear_and(C, history)
        else:
            h = len(bits) // 2
            (lo, hl), (hi, hh) = and_tree(C, bits[:h]), and_tree(C, bits[h:])
            C.ccx(lo, hi, target)
            clear_and(C, hh)
            clear_and(C, hl)
        for w in flips:
            C.x(w)
        return
    chain = []
    acc = ctrl if ctrl is not None else bits[0]
    factors = bits[:-1] if ctrl is not None else bits[1:-1]
    for w in factors:
        (a,) = C.alloc(1)
        C.ccx(acc, w, a)
        chain.append((acc, w, a))
        acc = a
    C.ccx(acc, bits[-1], target)
    for acc, w, a in reversed(chain):
        C.ccx(acc, w, a)
    C.free([a for _, _, a in chain])
    for w in flips:
        C.x(w)


def encode(P, n):
    x, y = (0, 0) if P is None else P
    return x | y << n


class ControlledPointAdd:
    def __init__(self, curve: Curve, A):
        assert A is not None and curve.on_curve(A)
        self.curve, self.A = curve, A
        self.mod = Mod(curve.q)
        self.n = self.mod.n
        neg = curve.neg
        twoA = curve.add(A, A)
        # (exceptional input, its sum with A); the generic kernel divides by zero on these.
        self.exceptions = [(None, A), (A, twoA), (neg(A), None), (neg(twoA), neg(A))]
        enc = [encode(e, self.n) for e, _ in self.exceptions]
        assert len(set(enc)) == 4, "addend order must exceed 3"
        self._square_blk = None

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * Y^2 mod q, Y untouched."""
        C, mod, n = self.C, self.mod, self.n
        if SquareFold.applies(mod.q):
            self._square_mulsub(gen, Yr, Xr)
            return
        T = C.alloc(n)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        if self._modadd_blk is None:
            self._modadd_blk = C.record_shape(lambda k, a, b: mod.modadd(C, a, b, ctrl=k[0]), [1, n, n])
            self._dbl = C.record_shape(lambda s: mod.moddbl(C, s), [n])
        for i in range(n):
            (k,) = C.alloc(1)
            C.ccx(gen, Yr[i], k)
            C.apply(self._modadd_blk, [[k], T, Xr], inverse=True)
            C.ccx(gen, Yr[i], k)
            C.free([k])
            C.apply(self._dbl, [T])
        for _ in range(n):
            C.apply(self._dbl, [T], inverse=True)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        C.free(T)

    def _square_mulsub(self, gen, Yr, Xr):
        C, mod, n = self.C, self.mod, self.n
        if self._square_blk is None:
            fold = SquareFold(mod)
            self._square_blk = C.record_shape(lambda y, j: fold.compute(C, y, j), [n, fold.junk])
            self._modadd_blk = C.record_shape(lambda k, a, b: fold.modadd(C, a, b, k[0]), [1, n, n])
            self._square_junk = fold.junk
        J = C.alloc(self._square_junk)
        C.apply(self._square_blk, [Yr, J])
        C.apply(self._modadd_blk, [[gen], J[:n], Xr], inverse=True)
        C.apply(self._square_blk, [Yr, J], inverse=True)
        C.free(J)

    def _kernel(self, c, gen, Xr, Yr, f0, f1):
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A

        def recode(inverse=False):
            if not inverse:
                C.cx(gen, c)
                C.cx(gen, f0)
            C.x(c)
            C.x(f1)
            (t,) = C.alloc(1)
            C.ccx(c, f0, t)
            C.ccx(t, f1, gen)
            C.ccx(c, f0, t)
            C.free([t])
            C.x(f1)
            C.x(c)
            if inverse:
                C.cx(gen, f0)
                C.cx(gen, c)

        mod.cmodadd_const(C, gen, q - a, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        C.ccx(f0, f1, Xr[0])
        recode()
        C.free([gen])
        self.M.apply(Xr, Yr, inverse=True)
        (gen,) = C.alloc(1)
        recode(inverse=True)
        mod.cmodadd_const(C, gen, 3 * a, Xr)
        self._mulsub(gen, Yr, Xr)
        recode()
        C.free([gen])
        self.M.apply(Xr, Yr)
        (gen,) = C.alloc(1)
        recode(inverse=True)
        C.ccx(f0, f1, Xr[0])
        mod.creflect(C, gen, a, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        return gen

    def _eq_flags(self, c, Xr, Yr, flags, points):
        C, n = self.C, self.n
        groups = {}
        for f, p in zip(flags, points):
            x, y = (0, 0) if p is None else p
            groups.setdefault(x, []).append((f, y))
        assert len(groups) == 3 and len(groups.get(0, [])) == 1, "0/a-pair/2A groups"
        pair_x = next(x for x, leaves in groups.items() if len(leaves) == 2)
        single_x = next(x for x in groups if x not in (0, pair_x))
        shared_zero = [i for i in range(n) if not ((pair_x | single_x) >> i & 1)]
        shared_one = [i for i in range(n) if (pair_x & single_x) >> i & 1]
        differing = [i for i in range(n) if (pair_x ^ single_x) >> i & 1]

        def chain(acc, wires):
            if EQ_TREE and wires:
                root, history = and_tree(C, wires)
                if acc is None:
                    return root, history
                (target,) = C.alloc(1)
                C.ccx(acc, root, target)
                return target, history + [(acc, root, target)]
            history = []
            for wire in wires:
                if acc is None:
                    acc = wire
                else:
                    (target,) = C.alloc(1)
                    C.ccx(acc, wire, target)
                    history.append((acc, wire, target))
                    acc = target
            return acc, history

        def clear(history):
            for left, right, target in reversed(history):
                C.ccx(left, right, target)
            C.free([target for _, _, target in history])

        # The zero leaf is uncontrolled; c joins only on the nonzero branches.
        for i in shared_zero:
            C.x(Xr[i])
        root, root_chain = chain(None, [Xr[i] for i in shared_zero])
        zero_bits = [Xr[i] for i in shared_one + differing]
        eq_const_into(C, root, zero_bits, 0, groups[0][0][0])
        common, common_chain = chain(root, [c] + [Xr[i] for i in shared_one])
        flips = [Xr[i] for i in differing if not (pair_x >> i & 1)]
        for wire in flips:
            C.x(wire)
        leaf, leaf_chain = chain(common, [Xr[i] for i in differing])
        for f, y in groups[pair_x]:
            eq_const_into(C, leaf, Yr[:1], y & 1, f)
        clear(leaf_chain)
        for wire in flips:
            C.x(wire)
        f, y = groups[single_x][0]
        bits = [Xr[i] for i in differing] + Yr[:1]
        value = sum(((single_x >> i) & 1) << j for j, i in enumerate(differing))
        value |= (y & 1) << len(differing)
        eq_const_into(C, common, bits, value, f)
        clear(common_chain)
        clear(root_chain)
        for i in shared_zero:
            C.x(Xr[i])

    def build(self, hierarchical=False):
        curve, A = self.curve, self.A
        twoA = curve.add(A, A)
        if A[0] == 0 or twoA is None or twoA[0] == 0:
            raise ValueError("kernel zero flag requires nonzero x(A) and x(2A)")
        if pow(curve.b, (curve.q - 1) // 2, curve.q) != curve.q - 1:
            raise ValueError("kernel zero flag refuses curves with a point at x = 0")
        self.C = C = Circuit(hierarchical)
        n = self.n
        (c,) = C.alloc(1)
        Xr, Yr = C.alloc(n), C.alloc(n)
        declared = Xr + Yr
        self.M = DialogMultiplier(C, self.mod)
        self._modadd_blk = self._square_blk = None

        flags = C.alloc(4)
        self._eq_flags(c, Xr, Yr, flags, [e for e, _ in self.exceptions])
        f0, f1, f2, gen = flags
        for f, (e_in, e_out) in zip(flags[1:], self.exceptions[1:]):
            load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f)
        # Move the uncontrolled x = 0 marker to f2 when c = 0.
        C.x(c)
        C.ccx(c, f0, f2)
        C.ccx(c, f2, f0)
        C.x(c)
        e_in, e_out = self.exceptions[0]
        load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f0)
        # Only the c-controlled part of f2 contributes to the generic flag.
        for w in (c, f0, f1):
            C.cx(w, gen)
        C.ccx(c, f2, gen)
        # Two wires name the one set flag among three; the wire left over is zero.
        C.cx(f2, f0)
        C.cx(f2, f1)
        C.ccx(f0, f1, f2)
        C.free([f2])

        gen = self._kernel(c, gen, Xr, Yr, f0, f1)

        (f2,) = C.alloc(1)
        C.ccx(f0, f1, f2)
        C.cx(f2, f1)
        C.cx(f2, f0)
        for w in (c, f0, f1):
            C.cx(w, gen)
        C.ccx(c, f2, gen)
        flags = [f0, f1, f2, gen]
        self._eq_flags(c, Xr, Yr, flags, [e for _, e in self.exceptions])
        C.free(flags)

        _restore_layout(C, Xr + Yr, declared)
        self.layout = {"control": [c], "x": declared[:n], "y": declared[n:]}
        return C


def _restore_layout(C, cur, declared):
    """SWAP logical bit i from wire cur[i] onto declared[i]. A displaced wire
    holds either another live bit, which is tracked, or a zero ancilla."""
    cur = list(cur)
    where = {w: i for i, w in enumerate(cur)}
    for i, d in enumerate(declared):
        s = cur[i]
        if s == d:
            continue
        C.swap(s, d)
        j = where.pop(d, None)
        del where[s]
        if j is not None:
            cur[j] = s
            where[s] = j
        cur[i] = d
        where[d] = i
