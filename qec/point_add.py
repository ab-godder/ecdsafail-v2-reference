"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, _maj, cswap, load_const
from .builder import Circuit
from .curve import Curve
from .square import HalfSquares, SquareFold
from .two_bit import TwoBitMod
from .model import tick_bound, tick_joint, tick_widths
from . import block0_code, tail_code, tape_block
from .tape_block import FREE, TICKS
from .builder import CCX, CX
from typing import NamedTuple

# K1: the CCX pairs, by op index within one Mod.cmodadd_const(ctrl, value, b) call, that equal CX
# from ctrl on every canonical b < q and either ctrl, keyed by (q, value mod q). scripts/k1_check.py
# proves each entry with its own BDD. Only the kernel's add of 3a uses it.
K1_CUTS = {
    (0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
     3 * 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
     % 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f): (2929, 5373),
}


def _k1_cut(C, start, q, value, ctrl):
    """Turn the K1 pair of the cmodadd_const call whose ops began at op start into CX from ctrl."""
    for pos in K1_CUTS.get((q, value % q), ()):
        j = 4 * (start + pos)
        op, x, y, t = C.prims[j:j + 4]
        assert op == CCX and ctrl in (x, y), (pos, op, x, y)
        C.prims[j:j + 4] = type(C.prims)("i", (CX, ctrl, t, t))


def _fan(C, src, copies):
    """CX tree that XORs src onto each copy, doubling the wires that hold it each layer. Returns
    the tree's (source, target) pairs in order; replaying them reversed undoes it."""
    tree, have, todo = [], [src], list(copies)
    while todo:
        layer = [(h, todo.pop(0)) for h in list(have)[:len(todo)]]
        for h, t in layer:
            C.cx(h, t)
        tree += layer
        have += [t for _, t in layer]
    return tree


def _rotate_left(C, ctrl, V, copies):
    """Under ctrl, V[:] = V[1:] + V[:1] on the wires: the serial strip's permutation, as two
    layers of disjoint swaps, all of V reversed and then all but its last wire. copies are clean
    wires that take ctrl for the layers' swaps and are clean again on return."""
    L = len(V)
    tree = _fan(C, ctrl, copies)
    ctl = [ctrl] + list(copies)
    for p in range(L // 2):
        cswap(C, ctl[p % len(ctl)], [V[p]], [V[L - 1 - p]])
    for p in range((L - 1) // 2):
        cswap(C, ctl[p % len(ctl)], [V[p]], [V[L - 2 - p]])
    for h, t in reversed(tree):
        C.cx(h, t)


class Stretch(NamedTuple):
    """Ticks start .. stop - 1 run their bits K and up as one high word on one label wire."""
    start: int
    stop: int
    K: int


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for canonical 1 <= x < q; see model.py."""

    # The tape stays decoded: no block, block-zero or tail code, and so no pebbled ticks. Every
    # code call below is gated on blocks or tail, and each encode was paired with a decode that
    # nothing between reads past, so the passes are the same maps on the same states.
    RAW_TAPE = False
    # Carries the forward pass keeps for the inverse pass to clear, each saving one CCX per pass.
    KEEP = 0
    # Give each cell 1 to T - 2 a wire for its swap bit; see DESIGN.md, "Swap bits on the tape".
    SWAP_WIRES = True
    # Give each cell a third wire for its old offset; see DESIGN.md, "Held offsets beside the tape".
    # Without it a tick strips on a fresh wire and the pass ends with the offset sweep.
    HELD = False
    # Above a tick's joint threshold K, U and V are never both 2^K or more, so a tick may run those
    # bits as one high word, 3 CCX a bit, on a label naming the register that holds it. Ticks go in
    # stretches of one K; each stretch leaves its label wire for the inverse pass. See high_plan.
    HIGH = True
    HIGH_MIN = 2
    # high_plan's gain model: CCX a high tick costs beyond its bits, and CCX a pass is charged for
    # each label wire, which is held through the Bezout pass at the multiplier's wire peak.
    HIGH_OVH = 5
    HIGH_PRICE = 285
    # Clean copies of g0*label for the high word's per-bit gates.
    HIGH_COPIES = 16
    # If set, the stretches to run instead of high_plan's; the checks use it.
    HIGH_STRETCHES = None

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
    def held_count(self):
        """Wires the Euclid pass takes for held offsets, one a cell with HELD."""
        return self.iters if self.HELD else 0

    @property
    def cell_wires(self):
        """Wires a Euclid step's cell register has without swap wires."""
        return 3 if self.HELD else 2

    @property
    def pass_tail_width(self):
        return self.widths[self.iters - 1]

    @property
    def tape_external(self):
        return 2 * self.iters - (self.mod.n - self.pass_tail_width) - self.u_lent

    @property
    def blocks(self):
        return 0 if self.RAW_TAPE else self.iters // TICKS

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

    def _block0_places(self, U, pool):
        """Where a merged pass leaves block 0's free wires: the last of its zero words, which fill
        U's top u_lent wires and then the pool's first lent_after_pass."""
        n = self.mod.n
        places = [(U, i) for i in range(n - self.u_lent, n)] + [(pool, i) for i in range(self.lent_after_pass)]
        return places[len(places) - len(self._block_free(0, self.swept)):]

    def _code_block0(self, U, pool):
        """Code block 0 of a pass that left it raw, with its free wires at _block0_places."""
        G = self.coded_tape(pool)
        for slot, (reg, i) in zip(self._free_slots(0, self.swept), self._block0_places(U, pool)):
            G[slot] = reg[i]
        self._encode(G, 0, self.swept)

    def _decode(self, G, j, wires, swept):
        assert len(wires) == len(self._block_free(j, swept))
        for slot, wire in zip(self._free_slots(j, swept), wires):
            G[slot] = wire
        self.C.apply(self._encoder(j, swept),
                     [self._block_wires(G, j, swept)], inverse=True)


    SELECT_COPIES = 16

    def _shared_difference(self, U, V, g0, g1, top=None, kept=None, copies=(), high=None):
        """Compare and form the even difference using one majority chain. With top set, bit w - 1
        is done by _top_bit on the domain that mode names. With kept, the chain's carries go onto
        those zero wires instead and stay there, which the inverse pass clears. copies are clean
        wires that take g1 and g0 for the selects and are clean again on return. With high =
        (m, z), the top m bits are done by _high_word on label z, with copies as its scratch."""
        if kept is not None:
            assert high is None
            return self._kept_difference(U, V, g0, g1, top, kept, copies)
        C = self.C
        (cin,) = C.alloc(1)
        carries = [cin] + U[:-1]
        for wire in V:
            C.x(wire)
        chain = len(U) if top is None else len(U) - 1
        if high is not None:
            assert top is None and 1 <= high[0] < len(U)
            chain = len(U) - high[0]
        # The majority chain, with each bit's two CX first: they read U[i] before the chain writes
        # it and only XOR into wires the chain also only XORs into, so the CCX ripple runs alone.
        C.cx(U[0], V[0])
        C.cx(U[0], carries[0])
        for i in range(1, chain):
            C.cx(U[i], V[i])
            C.cx(U[i], carries[i])
        C.cx(V[0], U[0])
        for i in range(1, chain):
            C.ccx(carries[i], V[i], U[i])
        if high is not None:
            self._high_word(U[chain:], V[chain:], carries[chain], g0, g1, high[1], copies)
        elif top is None:
            C.ccx(g0, U[-1], g1)
        else:
            self._top_bit(C, top, g0, U[-1], V[-1], U[-2], g1)
        # The selects go in before the unwind reaches each bit, so the unwind's serial chain
        # holds one CCX and one CX a bit. With S = carries[i] = c ^ a, B = V[i] = b ^ a: S ^= g1*B,
        # then B ^= g0*NOT(S). The restore A ^= S*B then leaves A = a ^ g1*(b ^ a), since
        # S*NOT(S) = 0, and S ^= A leaves S = c. Three CX after the chain give the old outputs.
        # This needs g1 => g0, which both g1 forms above give.
        # g1 and g0 are copied onto the clean copies, which take the selects round-robin, so the
        # selects are not one serial path on each; neither is written until the copies are cleared.
        k = min(self.SELECT_COPIES, len(copies) // 2)
        trees = [_fan(C, g1, copies[:k]), _fan(C, g0, copies[k:2 * k])]
        c1, c0 = [g1] + list(copies[:k]), [g0] + list(copies[k:2 * k])
        for j, i in enumerate(reversed(range(1, chain))):
            C.ccx(c1[j % len(c1)], V[i], carries[i])
            C.x(carries[i])
            C.ccx(c0[j % len(c0)], carries[i], V[i])
            C.x(carries[i])
        # The unwind, each bit's CX after the whole ripple: CX(U[i], carries[i]) only XORs into a
        # wire the next lower CCX also only XORs into, and comes after the CCX that reads it. Run
        # backward, as the inverse pass runs it, the CX then go first and the ripple runs alone.
        for i in reversed(range(1, chain)):
            C.ccx(carries[i], V[i], U[i])
        for i in reversed(range(1, chain)):
            C.cx(U[i], carries[i])
        # Bit 0 as before: c ^= a; v ^= a; v ^= c with the empty X pair gone is v ^= c; c ^= a.
        C.cx(V[0], U[0])
        C.cx(g0, V[0])
        C.cx(carries[0], V[0])
        C.cx(U[0], carries[0])
        C.cx(carries[0], V[0])
        for j, i in enumerate(range(1, chain)):
            C.cx(U[i], V[i])
            C.cx(c1[j % len(c1)], U[i])
            C.cx(c0[j % len(c0)], V[i])
        C.free([cin])
        # _high_word leaves its bits un-complemented already.
        for wire in V[:chain] if high is not None else V:
            C.x(wire)
        for j, wire in enumerate(V[1:]):
            C.cx(c1[j % len(c1)], wire)
        for tree in trees:
            for h, t in reversed(tree):
                C.cx(h, t)

    def _kept_difference(self, U, V, g0, g1, top, kept, copies=()):
        """_shared_difference with c_(i+1) built on kept[i - 1] by one CCX and never cleared, so
        U[i] keeps a_i and each bit's descent is the select and the exchange alone. Bit 0 is as
        there: U[0] holds c_1 = NOT g0 until its descent restores a_0 = 1. The returned carries
        are those of a + NOT(vhat) on the step's input, a function of that input.
        With fewer kept wires than chain bits, bits m + 1 and up, m = len(kept), use the in-place
        majority of _shared_difference on top of the kept c_(m+1), which its descent restores."""
        C = self.C
        chain = len(U) if top is None else len(U) - 1
        m = len(kept)
        assert 1 <= m <= chain - 1
        carry = [None, U[0]] + list(kept) + U[m + 1:chain]
        (cin,) = C.alloc(1)
        for wire in V:
            C.x(wire)
        C.cx(U[0], V[0])
        C.cx(U[0], cin)
        C.cx(V[0], U[0])
        for i in range(1, m + 1):
            # c_(i+1) = a ^ (a ^ b)(a ^ c), with a, b and c restored.
            C.cx(U[i], V[i])
            C.cx(U[i], carry[i])
            C.ccx(carry[i], V[i], carry[i + 1])
            C.cx(U[i], carry[i + 1])
            C.cx(U[i], carry[i])
            C.cx(U[i], V[i])
        for i in range(m + 1, chain):
            _maj(C, carry[i], V[i], U[i])
        if top is None:
            C.ccx(g0, carry[chain], g1)
        else:
            self._top_bit(C, top, g0, U[-1], V[-1], carry[chain], g1)
        # As in _shared_difference, the selects take g0 and g1 round-robin from clean copies.
        k = min(self.SELECT_COPIES, len(copies) // 2)
        trees = [_fan(C, g1, copies[:k]), _fan(C, g0, copies[k:2 * k])]
        c1, c0 = [g1] + list(copies[:k]), [g0] + list(copies[k:2 * k])
        for i in reversed(range(m + 1, chain)):
            # The in-place descent of _shared_difference, which restores c_i on carry[i].
            C.ccx(carry[i], V[i], U[i])
            C.ccx(c0[i % len(c0)], carry[i], V[i])
            C.cx(U[i], carry[i])
            C.cx(U[i], V[i])
            C.cx(carry[i], V[i])
            C.x(V[i])
            C.ccx(c1[i % len(c1)], V[i], U[i])
            C.x(V[i])
            C.cx(carry[i], V[i])
        for i in reversed(range(1, m + 1)):
            # B[i] = b ^ g0*(a ^ c), S_i when g0 = 1; then a_i ^= g1 * NOT(S_i ^ c_i).
            C.cx(U[i], carry[i])
            C.ccx(c0[i % len(c0)], carry[i], V[i])
            C.cx(U[i], carry[i])
            C.cx(carry[i], V[i])
            C.x(V[i])
            C.ccx(c1[i % len(c1)], V[i], U[i])
            C.x(V[i])
            C.cx(carry[i], V[i])
        C.cx(V[0], U[0])
        C.cx(g0, V[0])
        C.cx(U[0], cin)
        C.cx(U[0], V[0])
        C.cx(cin, V[0])
        C.cx(cin, V[0])
        C.free([cin])
        for wire in V:
            C.x(wire)
        for j, wire in enumerate(V[1:]):
            C.cx(c1[j % len(c1)], wire)
        for tree in trees:
            for h, t in reversed(tree):
                C.cx(h, t)

    def _high_copies(self, m):
        return min(self.HIGH_COPIES, m // 2)

    def _high_scratch(self, m):
        """Clean wires _high_word borrows: the borrows E_0 .. E_m, g0*z and its copies."""
        return m + 2 + self._high_copies(m)

    def _high_word(self, U, V, cw, g0, g1, z, scratch):
        """The chain's top m bits, where U or V is zero (the joint threshold). On entry V holds
        NOT v, cw the carry c into these bits, and g1 is 0. z names the register holding the word
        X = U ^ V: 1 for U, 0 for V, either when X = 0. A subtract takes one borrow
        e = g0*(c ^ z) off X, the word of the larger minus the smaller, so V's bits end as
        X - e and U's are cleared when z: the borrow chain E_(i+1) = E_i AND NOT X_i, up and down,
        and one CCX a bit for U, 3 in place of the majority cell's 4. The swap is g0*z, or the
        borrow leaving the top, g1 = g0*z ^ E_m. V's bits are left so that the caller's CX from
        g1, with no NOT, gives X - e, the top borrow's all-ones word back to 0. scratch goes back
        clean."""
        C, m = self.C, len(U)
        E, zg = scratch[:m + 1], scratch[m + 1]
        zc = [zg] + list(scratch[m + 2:m + 2 + self._high_copies(m)])
        for u, v in zip(U, V):
            C.cx(u, v)
        C.ccx(g0, z, zg)
        tree = _fan(C, zg, zc[1:])
        # U ^= zg * X, with V = NOT X: zg ^ zg * V.
        for i, (u, v) in enumerate(zip(U, V)):
            C.cx(zc[i % len(zc)], u)
            C.ccx(zc[i % len(zc)], v, u)
        C.cx(z, cw)
        C.ccx(g0, cw, E[0])
        C.cx(z, cw)
        for i, v in enumerate(V):
            C.ccx(E[i], v, E[i + 1])
        C.cx(zg, g1)
        C.cx(E[m], g1)
        # The borrows go in before the chain comes down, so neither ripple waits on a CX: with
        # V = X ^ E, the borrowed word, E_(i+1) = E_i AND NOT X_i is E_i AND V_i. The NOT sits
        # between CX from E_i, so in both directions E_i and V_i reach each chain's CCX at the same
        # tick and the depth pass orients those CCX by the chain alone.
        for i, v in enumerate(V):
            C.cx(E[i], v)
            C.x(v)
            C.cx(E[i], v)
            C.cx(E[i], v)
        for i, v in reversed(list(enumerate(V))):
            C.ccx(E[i], v, E[i + 1])
        C.cx(z, cw)
        C.ccx(g0, cw, E[0])
        C.cx(z, cw)
        for i, (u, v) in enumerate(zip(U, V)):
            C.cx(u, v)
            C.cx(zc[i % len(zc)], v)
        for h, t in reversed(tree):
            C.cx(h, t)
        C.ccx(g0, z, zg)

    def _or_into(self, bits, target, scratch):
        """target ^= OR of bits, as a tree of OR gates on scratch, which goes back clean."""
        C = self.C

        def gate(a, b, t):
            C.x(a)
            C.x(b)
            C.ccx(a, b, t)
            C.x(a)
            C.x(b)
            C.x(t)

        if len(bits) == 1:
            C.cx(bits[0], target)
            return
        free, steps, level = list(scratch), [], list(bits)
        while len(level) > 2:
            nxt = []
            for k in range(0, len(level) - 1, 2):
                t = free.pop(0)
                steps.append((level[k], level[k + 1], t))
                nxt.append(t)
            if len(level) % 2:
                nxt.append(level[-1])
            level = nxt
        for a, b, t in steps:
            gate(a, b, t)
        gate(level[0], level[1], target)
        for a, b, t in reversed(steps):
            gate(a, b, t)

    @staticmethod
    def _top_bit(C, mode, g0, A, B, Cw, H):
        """Bit w - 1 of the shared difference in place of its MAJ, g1 and unwind (5 CCX). On entry
        A = U[w-1] holds a, B = V[w-1] holds b = NOT vhat_top, Cw = U[w-2] holds the carry c into
        bit w - 1, H = g1 is 0. On exit A = U'_top, B = b ^ g0*(a ^ c), Cw = c, H = g0*MAJ(a, b, c),
        as the five CCX leave them, on the domain the mode names:
        top_clear: b = 1 (V's top wire is 0 on entry), so g1 = g0*(a | c) and U'_top = a & ~g0.
        narrow: the step leaves U's top wire 0, so a implies g0 & b, g1 = a ^ g0*c*(a ^ b) and
            B = a ^ b ^ g0*c.
        both: both hold."""
        if mode == "top_clear":
            C.cx(Cw, A)
            C.ccx(g0, A, B)
            C.cx(Cw, A)
            C.ccx(g0, A, H)
            C.cx(H, A)
            C.x(B)
            C.ccx(Cw, B, H)
            C.x(B)
        elif mode == "narrow":
            C.ccx(g0, Cw, H)
            C.cx(A, B)
            C.ccx(H, B, A)
            C.cx(H, B)
            C.ccx(g0, Cw, H)
            C.cx(A, H)
            C.cx(H, A)
        elif mode == "both":
            # P1 = ~a & ~c on B, P2 = ~g0 & ~c on A; then g1 = g0 ^ P1 ^ P2, B = a ^ g0 ^ c ^ P2.
            C.x(B)
            C.x(A)
            C.x(Cw)
            C.ccx(A, Cw, B)
            C.x(A)
            C.x(Cw)
            C.cx(A, B)
            C.x(g0)
            C.x(Cw)
            C.ccx(g0, Cw, A)
            C.x(g0)
            C.x(Cw)
            C.cx(A, H)
            C.cx(B, H)
            C.cx(g0, H)
            C.cx(H, B)
            C.cx(Cw, B)
            C.cx(B, A)
            C.cx(g0, A)
            C.cx(Cw, A)
        else:
            raise ValueError(mode)

    @staticmethod
    def _kept_count(w, top_clear, narrow, last):
        """Carries a step builds on kept wires: one per chain bit above bit 0."""
        if last:
            return 0
        top = w >= 2 and (top_clear or narrow)
        return (w - 1 if top else w) - 1

    def kept_plan(self):
        """Per tick, how many carries it keeps: narrow ticks from the last back, within KEEP. The
        first tick that does not fit keeps its lowest KEEP - total carries, if any."""
        plan, total = [0] * self.iters, 0
        for i in reversed(range(1, self.iters)):
            w = self.widths[i]
            if w == self.mod.n:
                break
            top_clear = w >= 3 and self.widths[i - 1] == w
            last = i == self.iters - 1
            narrow = not last and self.widths[i + 1] == w - 1
            k = self._kept_count(w, top_clear, narrow, last)
            if total + k > self.KEEP:
                plan[i] = self.KEEP - total
                break
            plan[i], total = k, total + k
        return plan

    def _spare_counts(self):
        """Per tick, the clean wires the pass has for the step once the tick has its cell."""
        n, W = self.mod.n, self.widths
        free, out = self.tape_external, []
        for i in range(self.iters):
            free -= 2
            out.append(free)
            nxt = W[i + 1] if i + 1 < self.iters else self.pass_tail_width
            free += W[i] - nxt + sum(1 for j in range(nxt, W[i]) if j >= n - self.u_lent)
        return out

    def _high_fits(self, t, K, spare):
        m = self.widths[t] - K
        return (1 <= t < self.iters - 1 and m >= max(1, self.HIGH_MIN)
                and self._high_scratch(m) <= spare[t])

    def high_plan(self):
        """The stretches of ticks that run their bits K and up as a high word. A stretch's K is
        its first tick's joint threshold; later ticks' thresholds are no higher, so it holds for
        them. The first tick sets the label wire to [U >= 2^K] by an OR, and every tick XORs its
        g1 into it: a swap leaves U below 2^K, so the word is V's after it, or there is none.
        No other tick changes which register holds it, and K never moves inside a stretch, so the
        label stays right. It is left as the stretch's junk wire for the inverse pass. A dynamic
        program picks the stretches that maximise the gain model: per tick, its high bits less
        HIGH_OVH; per stretch, less the OR and HIGH_PRICE."""
        if self.HIGH_STRETCHES is not None:
            return [Stretch(*s) for s in self.HIGH_STRETCHES]
        if not self.HIGH or self.HELD or self.KEEP:
            return []
        T, W = self.iters, self.widths
        joint, spare = tick_joint(self.mod.q, T), self._spare_counts()
        best, back = [0] * (T + 1), [None] * (T + 1)
        for a in range(T):
            if best[a] > best[a + 1]:
                best[a + 1], back[a + 1] = best[a], None
            K = joint[a]
            if not self._high_fits(a, K, spare):
                continue
            acc = -(2 * (W[a] - K) - 3) - self.HIGH_PRICE
            b = a
            while b < T and self._high_fits(b, K, spare):
                acc += W[b] - K - self.HIGH_OVH
                if best[a] + acc > best[b + 1]:
                    best[b + 1], back[b + 1] = best[a] + acc, a
                b += 1
        out, t = [], T
        while t > 0:
            if back[t] is None:
                t -= 1
            else:
                out.append(Stretch(back[t], t, joint[back[t]]))
                t = back[t]
        return out[::-1]

    @staticmethod
    def _strip_copies(w, top_clear):
        """Clean wires the strip's rotation of _strip_len wires takes for its control."""
        return max(0, DialogMultiplier._strip_len(w, top_clear) // 2 - 1)

    @staticmethod
    def _strip_len(w, top_clear):
        return w - 2 if top_clear else w - 1

    def _euclid_step(self, w: int, top_clear: bool = False, narrow: bool = False, last: bool = False,
                     keep: int = 0, sw: bool = False, cp: int = 0, hi=None):
        """top_clear: the tick before ran at the same width, so V is below 2^(w - 1) on entry
        and the second strip's last swap, of two zero bits, never fires.
        narrow: the next width is w - 1, so by E7 this step leaves U below 2^(w - 1).
        last: the pass's last tick, whose input reaches (U, V, E) = (1, 0, 0) by E3.
        sw: the cell has a fourth wire, zero on entry, which the step leaves holding s = g0 AND b.
        cp: the step takes that many clean wires, as copies of the strip's control and scratch,
        and returns them clean.
        hi: (m, init) for a tick of a stretch; see high_plan. The step takes the stretch's label
        wire last. init: the stretch's first tick, which sets the label to [U >= 2^(w - m)]."""
        if last:
            top_clear = narrow = False
            keep = 0
            sw = False
            cp = 0
            hi = None
        assert hi is None or (not self.HELD and not keep and cp >= self._high_scratch(hi[0]))
        key = (w, top_clear, narrow, last, keep, sw, cp, hi)
        if key not in self._euclid_steps:
            C = self.C

            def step(U, V, E, cell, *extra):
                rest = list(extra)
                kept = (rest.pop(0),) if keep else ()
                copies = rest.pop(0) if cp else []
                label = rest.pop(0)[0] if hi else None
                g0, b = cell[:2]
                offset = E[0]
                if not self.HELD:
                    self._unheld_step(U, V, E, cell, top_clear, narrow, last, kept, copies,
                                      hi and (hi[0], hi[1], label))
                    return
                held = cell[2]
                if last:
                    # The input is (1, 0, 0), (1, 1, 0) or (1, 0, 1). The cell becomes
                    # (V0 ^ E, NOT g0) with held = E, and V and E end zero: all by CX.
                    C.cx(V[0], g0)
                    C.cx(offset, g0)
                    C.cx(g0, V[0])
                    C.cx(offset, V[0])
                    C.cx(g0, b)
                    C.x(b)
                    C.cx(offset, held)
                    C.cx(held, offset)
                    V[:] = V[1:] + V[:1]
                    return
                top = None
                if w >= 2 and (top_clear or narrow):
                    top = "both" if top_clear and narrow else "top_clear" if top_clear else "narrow"
                C.cx(offset, V[0])
                C.cx(V[0], g0)
                self._shared_difference(U[:w], V[:w], g0, b, top, kept[0] if keep else None, copies)
                V[:] = V[1:] + V[:1]
                # A no-subtract tick strips a second bit when V + E is even. held, zero on
                # entry, is the strip's control: two = NOT g0 AND NOT (V0 ^ E), E still e.
                C.cx(offset, V[0])
                C.x(g0)
                C.x(V[0])
                C.ccx(g0, V[0], held)
                C.x(V[0])
                C.x(g0)
                C.cx(offset, V[0])
                C.ccx(held, offset, V[0])
                # Under held, rotate V[:L] left by one: L - 1 swaps, as _rotate_left's two layers.
                _rotate_left(C, held, V[:self._strip_len(w, top_clear)],
                             copies[:self._strip_copies(w, top_clear)])
                # b = s ^ two, the cell's swept bit; then held = s, as two ^ (s ^ two).
                C.cx(held, b)
                C.cx(b, held)
                if sw:
                    C.cx(held, cell[3])
                # A subtract keeps its swap s on E and holds the old offset e beside the cell.
                cswap(C, g0, [offset], [held])

            sizes = [w, w, 1, (4 if sw else 3) if self.HELD else 2] + ([keep] if keep else []) + \
                ([cp] if cp else []) + ([1] if hi else [])
            self._euclid_steps[key] = C.record_shape(step, sizes)
        return self._euclid_steps[key]

    def _unheld_step(self, U, V, E, cell, top_clear, narrow, last, kept, copies, hi=None):
        """The Euclid step on a two-wire cell, as before the held offsets: a subtract keeps its swap
        on E and moves the old offset onto b, and the strip runs on a fresh wire that the step
        clears. The pass's offset sweep then leaves the same swept tape as the held step.
        hi: (m, init, z) runs the top m bits as a high word on label z, which then takes g1."""
        C, w = self.C, len(V)
        g0, b = cell
        offset = E[0]
        if last:
            # The input is (1, 0, 0), (1, 1, 0) or (1, 0, 1), and the cell becomes
            # (V0 ^ E, NOT V0): only V[0], E and the cell change, all by CX.
            C.cx(V[0], g0)
            C.cx(offset, g0)
            C.x(b)
            C.cx(V[0], b)
            C.cx(b, V[0])
            C.x(V[0])
            C.cx(g0, offset)
            C.cx(b, offset)
            C.x(offset)
            V[:] = V[1:] + V[:1]
            return
        top = None
        if w >= 2 and (top_clear or narrow):
            top = "both" if top_clear and narrow else "top_clear" if top_clear else "narrow"
        high = None
        if hi:
            m, init, z = hi
            top, high = None, (m, z)
            if init:
                self._or_into(U[w - m:w], z, copies[:m - 2])
        C.cx(offset, V[0])
        C.cx(V[0], g0)
        self._shared_difference(U[:w], V[:w], g0, b, top, kept[0] if kept else None, copies, high)
        if hi:
            # The label's update: a swap moves the word to V or leaves both registers below 2^K.
            C.cx(b, z)
        V[:] = V[1:] + V[:1]
        # A no-subtract tick strips a second bit when V + E is even. two implies NOT g0, and
        # the subtract's swap below is the identity when g0 = 0, so the strip and the clear of
        # two read the same E and b before that swap as after it, and go first.
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
        _rotate_left(C, two, V[:self._strip_len(w, top_clear)],
                     copies[:self._strip_copies(w, top_clear)])
        C.x(g0)
        C.ccx(g0, b, two)
        C.x(g0)
        C.free([two])
        # A subtract keeps its swap on E and moves the old offset to the cell.
        cswap(C, g0, [offset], [b])

    def _euclid_tick0(self):
        """First Euclid step on U = q, E = 0. U is zero once q's bits flip, so it takes the carries
        of x + g0*f1, which then become U = x without being cleared (Joel's 4748bb3)."""
        if self._euclid_tick0_block is None:
            C, q, n = self.C, self.mod.q, self.mod.n
            f1 = (1 << n) - q + 1
            j = (f1 & -f1).bit_length() - 1

            def step(U, V, E, cell):
                g0, b = cell
                C.cx(V[0], g0)
                for i in range(n):
                    if q >> i & 1:
                        C.x(U[i])
                # U[i] = g0*c_i, c_i the carry into bit i of x + f1, which has no bit below j.
                # Where f1_i = 1, g0*(x_i OR c_i) = d ^ x_i*(g0 ^ d) with d = g0*c_i.
                for i in range(j, n - 1):
                    if i == j:
                        C.ccx(g0, V[i], U[i + 1])
                    elif f1 >> i & 1:
                        C.cx(U[i], U[i + 1])
                        C.cx(g0, U[i])
                        C.ccx(U[i], V[i], U[i + 1])
                        C.cx(g0, U[i])
                    else:
                        C.ccx(U[i], V[i], U[i + 1])
                # V ^= g0*(f1 ^ c) makes V = x + g0*f1; then U ^= V ^ g0*f1 leaves U = x.
                for i in range(j, n):
                    if i > j:
                        C.cx(U[i], V[i])
                    if f1 >> i & 1:
                        C.cx(g0, V[i])
                for i in range(n):
                    C.cx(V[i], U[i])
                    if f1 >> i & 1:
                        C.cx(g0, U[i])
                # The select: U = q > x, so a subtract swaps and U keeps x; off g0, V still holds x
                # and U takes q back, U ^= NOT g0*(x ^ q). x_0 = g0 and q_0 = 1.
                C.x(g0)
                C.cx(g0, U[0])
                for i in range(1, n):
                    if q >> i & 1:
                        C.x(V[i])
                    C.ccx(g0, V[i], U[i])
                    if q >> i & 1:
                        C.x(V[i])
                C.x(g0)
                for wire in V:
                    C.cx(g0, wire)
                C.cx(g0, E[0])
                V[:] = V[1:] + V[:1]
                C.x(g0)
                C.x(V[0])
                C.ccx(g0, V[0], b)
                C.x(V[0])
                C.x(g0)
                for i in range(n - 2):
                    cswap(C, b, [V[i]], [V[i + 1]])
                # The cell's swept bit s ^ two, with s = g0 since U = q > v; e = 0 needs no held.
                # Without held the offset sweep does this at cell 0.
                if self.HELD:
                    C.cx(g0, b)

            self._euclid_tick0_block = C.record_shape(step, [n, n, 1, 2])
        return self._euclid_tick0_block

    def _euclid_pass(self, U: list[int], pool: list[int], held: list[int],
                     block0_raw: bool = False, kept: list[int] | None = None,
                     sweep: bool = True, swr: list[int] | None = None,
                     labels: list[int] | None = None) -> None:
        # The pool starts as V plus clean external tape. It ends as the zero wires U did not
        # take back, then tape; U ends as its low wires and those zero wires. The boundary is a
        # permutation of U and the pool together. held, iters wires, starts zero and ends with
        # each subtract cell's old offset, zero elsewhere; only the inverse pass clears it.
        # swr, if given, starts zero and ends with cells 1 to T - 2's swap bits s = g0 AND b.
        # labels, one wire a stretch of high_plan, starts zero and ends with each stretch's label;
        # only the inverse pass clears it. A caller that gives none gets them as fresh wires.
        n = self.mod.n
        assert swr is None or len(swr) == self.iters - 2
        assert len(held) == self.held_count
        assert self.HELD or swr is None
        plan = self.kept_plan() if kept is not None else [0] * self.iters
        assert kept is None or len(kept) == sum(plan)
        stretches = self.high_plan()
        if labels is None:
            labels = self.C.alloc(len(stretches))
        assert len(labels) == len(stretches)
        joint, spare_at, high = tick_joint(self.mod.q, self.iters), self._spare_counts(), {}
        for j, s in enumerate(stretches):
            for t in range(s.start, s.stop):
                assert t not in high and s.K >= joint[t] and self._high_fits(t, s.K, spare_at)
                assert not plan[t] and not self.HELD
                high[t] = (self.widths[t] - s.K, t == s.start, j)
        at = 0
        V, available, G = pool[:n], pool[n:], []
        lend = U[n - self.u_lent:]
        E = self.C.alloc(1)
        for i in range(self.iters):
            w = self.widths[i]
            assert len(V) == w and len(available) >= 2
            cell, available = available[:2], available[2:]
            tick0 = self._euclid_tick0() if i == 0 and w == n else None
            if tick0 is not None and tick0.wires <= set(tick0.flat_in) and set(tick0.flat_out) == set(tick0.flat_in):
                self.C.apply(tick0, [U, V, E, cell])
            else:
                # A tick leaves V below 2^(w' - 1), w' its own width.
                top_clear = i >= 1 and w >= 3 and self.widths[i - 1] == w
                # E7: the next tick's U fits its width. E3: the last tick ends at (1, 0, 0).
                last = i == self.iters - 1
                narrow = not last and self.widths[i + 1] == w - 1
                cell = cell + held[i:i + 1]
                assert len(cell) == self.cell_wires
                on_sw = swr is not None and 1 <= i <= self.iters - 2
                if on_sw:
                    cell = cell + [swr[i - 1]]
                regs = [U[:w], V, E, cell]
                if plan[i]:
                    assert plan[i] <= self._kept_count(w, top_clear, narrow, last)
                    regs.append(kept[at:at + plan[i]])
                    at += plan[i]
                # The strip's control copies go on later cells' held wires, which are zero until
                # their own steps. Without held they go on the pass's unwritten wires, the last of
                # available, which are zero until later steps take them.
                spare = self.iters - 1 - i if self.HELD else len(available)
                cp = 0 if last else min(max(self._strip_copies(w, top_clear), 2 * self.SELECT_COPIES),
                                        spare)
                hi = high.get(i)
                if hi:
                    # The high word's scratch is clean again before the copies are taken.
                    assert spare == spare_at[i]
                    cp = max(cp, self._high_scratch(hi[0]))
                if cp:
                    regs.append(held[self.iters - cp:] if self.HELD else available[len(available) - cp:])
                if hi:
                    regs.append(labels[hi[2]:hi[2] + 1])
                self.C.apply(self._euclid_step(w, top_clear, narrow, last, plan[i], on_sw, cp,
                                               hi and hi[:2]), regs)
            G.extend(cell[:2])
            next_width = self.widths[i + 1] if i + 1 < self.iters else self.pass_tail_width
            assert w - 1 <= next_width <= w
            available.extend(V[next_width:])
            V[:] = V[:next_width]
            # By L6, u < 2**next_width from here on, and later steps touch only U[:w] for
            # w <= next_width, so U's wires from next_width up are zero for the rest of the pass.
            available.extend(reversed([u for u in U[next_width:w] if u in lend]))
        assert len(V) == self.pass_tail_width and not available
        # Each tick leaves its own swap on a subtract cell's b, so only the B1 sweep is left.
        # Sweep down, coding each block as it finishes, once the B1 step that reads its
        # bottom cell has run.
        T, swept, clean = self.iters, self.swept, []
        # Without sweep the pass leaves the tape unswept, for a Bezout pass that sweeps it inline.
        assert sweep or not self.blocks

        def done(i):
            k = i + 1 if swept else i
            if k % TICKS == 0 and k // TICKS < self.blocks:
                clean.extend(self._encode(G, k // TICKS, swept))

        self._b1_sweep(G, swept and sweep, done, None if self.HELD else E)
        self.C.free(E)
        # Block 0's free wires take their places among the zero words raw, and block 0 is coded
        # last, unless block0_raw leaves that to the caller.
        raw = []
        if swept and self.blocks:
            raw = [G[slot] for slot in self._free_slots(0, swept)]
            for slot in self._free_slots(0, swept):
                G[slot] = None
        zeros = V + clean + raw
        U[:] = U[:n - self.u_lent] + zeros[:self.u_lent]
        pool[:] = zeros[self.u_lent:] + [w for w in G if w is not None]
        if raw and not block0_raw:
            self._code_block0(U, pool)

    def _b1_sweep(self, G, b1=True, done=None, E=None):
        """Rotate each no-subtract run of the tape G so B1 is last in Bezout order, cell T - 1 down
        to 0, calling done(i) once cell i's steps are over; with b1 False only done runs.
        With E, the offset wire of a pass without held offsets, each cell first takes the offset
        sweep, which moves each subtract cell's swap from E onto b and leaves E zero."""
        T = self.iters
        both = None

        def mark(k):
            self.C.x(G[2 * k])
            self.C.x(G[2 * k - 2])
            self.C.ccx(G[2 * k], G[2 * k - 2], both)
            self.C.x(G[2 * k - 2])
            self.C.x(G[2 * k])

        # The pass reaches (1, 0) within T ticks: its last subtract leaves E = 0, and each
        # no-subtract cell after it has b = 1. So when cells T - 1 and T - 2 are both
        # no-subtract their b are equal: the B1 step between them never swaps, and mark(T - 2)
        # is made at once.
        late = b1 and T >= 3
        if b1 and T >= 2:
            (both,) = self.C.alloc(1)
            mark(T - 2 if late else T - 1)
        for i in range(T - 1, -1, -1):
            if E is not None:
                # At cell T - 1, g0 = 0 implies E ^ b = 1; at cell 0, tick 0 on U = q > x, a
                # subtract always swaps, so there g0 = 1 implies E ^ b = 1.
                g0, b = G[2 * i], G[2 * i + 1]
                if i == 0 or i == T - 1:
                    self.C.cx(b, E[0])
                    self.C.cx(g0, b)
                    if i:
                        self.C.cx(E[0], b)
                        self.C.x(b)
                    self.C.cx(b, E[0])
                else:
                    cswap(self.C, g0, E, [b])
            if b1 and i + 1 < T and not (late and i == T - 2):
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
            if done is not None:
                done(i)
        if both is not None:
            self.C.free([both])

    def _b1_sweep_sw(self, G, swr):
        """The B1 sweep that a Bezout pass with swap wires applies to the tape, cell T - 2 down to
        0, as TwoBitMod.tick applies it: one CCX a cell, since s = g0 AND x on swr makes B1 linear.
        Its inverse unsweeps the tape that pass leaves."""
        C, T, tw = self.C, self.iters, self.narrow
        (tau,) = C.alloc(1)
        for i in range(T - 2, -1, -1):
            g0, b = G[2 * i], G[2 * i + 1]
            up = None if i == T - 2 else G[2 * i + 2]
            C.cx(tau, b)
            if i == 0:
                # t_1 = NOT g0_1 AND B1_0, and B1_0 = NOT x_0. t_0 = 0.
                C.x(b)
                C.x(up)
                C.ccx(b, up, tau)
                C.x(up)
                C.x(b)
            else:
                tw._sweep_in(C, g0, b, swr[i - 1], up, G[2 * i - 2], tau)
        C.free([tau])

    def _bezout(self):
        if self._bezout_step is None:
            C, mod = self.C, self.mod

            def step(g, R, S):
                self.narrow.step(C, g, R, S)

            self._bezout_step = C.record_shape(step, [2, mod.n, mod.n])
        return self._bezout_step

    def _bezout_tick(self, first, last, pebble=False, inline=False, tape=False, split=None,
                     after=False):
        """after marks the tick after the terminal cell, which takes that cell's g0 before flag."""
        key = (first, last, pebble, inline, tape, split, after)
        if key not in self._bezout_ticks:
            C, n = self.C, self.mod.n
            if tape:
                def tick(prev, g, nxt, R, S, tau, psw, gsw, hand):
                    # gsw is this cell's swap wire, psw the previous cell's, which holds its
                    # B1 mark on entry and its swap bit again on return.
                    self.narrow.tick(C, prev or None, g, nxt[0] if nxt else None, R, S,
                                     pebble=pebble, held=psw[0] if psw else None, cell0=last,
                                     hand=hand, tau=tau[0], sw=gsw[0] if gsw else None, tape=True,
                                     split=split)

                sizes = [0 if first else 2, 2, 0 if last else 1, n, n, 1,
                         0 if first else 1, 0 if last else 1, 2]
                self._bezout_ticks[key] = C.record_shape(tick, sizes)
                return self._bezout_ticks[key]

            def tick(prev, g, nxt, R, S, *rest):
                *tau, flag, hand = rest
                term = tau.pop() if after else None
                # flag comes in holding [previous cell is B1], or 0 on a first tick, which has no
                # previous cell, and goes out holding [this cell is B1], or empty after the last.
                # hand carries a B1's top two bits to the next tick the same way, zero elsewhere.
                # tau, with inline, carries the B1 sweep's toggle from tick to tick.
                out = self.narrow.tick(C, prev or None, g, nxt[0] if nxt else None, R, S,
                                       pebble=pebble, held=None if first else flag[0],
                                       cell0=last, hand=hand, tau=tau[0][0] if tau else None,
                                       split=split, term=term[0] if after else None)
                flag[:] = [] if out is None else [out]

            # flag and hand come last: a last tick's output leaves them empty.
            sizes = ([0 if first else 2, 2, 0 if last else 1, n, n] + ([1] if inline else [])
                     + ([1] if after else []) + [1, 2])
            self._bezout_ticks[key] = C.record_shape(tick, sizes)
        return self._bezout_ticks[key]

    def _split(self, i):
        """The add-or-shift of Bezout tick i: None for the default, False for one of SPLIT_SKIP
        ticks spread evenly over the pass."""
        skip, T = self.narrow.SPLIT_SKIP, self.iters
        if not (self.narrow.SPLIT_ADD and skip) or (i * skip) // T == ((i + 1) * skip) // T:
            return None
        return False

    @property
    def inline_sweep(self):
        """The Bezout pass sweeps B1 itself, on a tape the Euclid pass leaves unswept."""
        return self.swept and not self.blocks

    def _b1_step(self, G, i, both):
        """Use the pass-wide mark while both cells are still raw."""
        cswap(self.C, both, [G[2 * i + 1]], [G[2 * i - 1]])

    def _bezout_pass(self, G, Yr, S, block0_raw=False, swr=None):
        """G is the coded tape. Each block is decoded on a clean wire before a tick reads it and
        coded again once no tick reads it, except block 0 on the merged path with block0_raw."""
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
            self._decode(G, self.blocks - 1,
                         C.alloc(len(self._block_free(self.blocks - 1, swept))), swept)
            if self.blocks == 1:
                block0_code.hold(C, self._block_wires(G, 0, swept))
                C.free([G[1]])
                G[1] = None
        # The final swept cell is (1, 0) or (0, 1): the terminal subtract
        # is followed only by no-subtract cells. Lend its b only when uncoded.
        lend_tail = T > 1 and self.blocks * TICKS < T
        tail = (not self.RAW_TAPE and lend_tail and T >= tail_code.CELLS + 1
                and T - self.blocks * TICKS >= tail_code.CELLS
                and tick_widths(self.mod.q, T)[T - tail_code.CELLS] <= 6)
        lo = 2 * (T - tail_code.CELLS)
        tail_encoded = False
        flag, hand = [], []
        inline = self.inline_sweep
        tau = C.alloc(1) if inline else []
        # With swr, cell i's swap bit s = g0 AND b sits on swr[i - 1] for 1 <= i <= T - 2; cells
        # T - 1 and 0 need none, since s is 0 on the terminal cell and g0 on cell 0.
        assert swr is None or (inline and lend_tail and not tail and T >= 3 and len(swr) == T - 2)
        for i in reversed(range(T)):
            nxt_block = i // TICKS - 1
            if i % TICKS == 0 and 0 <= nxt_block < self.blocks:
                self._decode(G, nxt_block,
                             C.alloc(len(self._block_free(nxt_block, swept))), swept)
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
            elif swr is not None:
                if not hand:
                    hand = C.alloc(2)
                psw = [] if first else [swr[i]]
                gsw = [] if last else [swr[i - 1]]
                C.apply(self._bezout_tick(first, last, pebble, inline, tape=True,
                                          split=self._split(i)),
                        [prev, G[2 * i:2 * i + 2], nxt, Yr, S, tau, psw, gsw, hand])
            else:
                if not flag:
                    # The pass's first tick call: no tick before it handed a flag over.
                    flag = C.alloc(1)
                    if not first:
                        self.narrow._mark_b1(C, prev, flag[0])
                if not hand:
                    # No tick before this one handed a pair over, so its cell is not B1.
                    hand = C.alloc(2)
                # The tick after the terminal cell enters with S = g0*R, that cell's g0.
                after = no_prev and not last
                C.apply(self._bezout_tick(first, last, pebble, inline, split=self._split(i),
                                          after=after),
                        [prev, G[2 * i:2 * i + 2], nxt, Yr, S] + ([tau] if inline else [])
                        + ([[G[2 * T - 2]]] if after else []) + [flag, hand])
            prev_block = (i + 1) // TICKS
            if (i + 1) % TICKS == 0 and i + 1 < T and prev_block < self.blocks:
                C.free(self._encode(G, prev_block, swept))
            if tail and i == T - tail_code.CELLS - 1:
                W = G[lo:lo + tail_code.LIVE]
                tail_code.apply(C, W)
                tail_encoded = True
                C.free([W[tail_code.FREE[0]]])
                G[lo + tail_code.FREE[0]] = None
        assert not flag and not hand
        C.free(tau)
        if self.blocks and not block0_raw:
            C.free(self._encode(G, 0, swept))
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

        # On the merged path block 0 stays decoded from Bezout's end into the inverse pass, so
        # the pass is recorded without its last encode, which runs here instead.
        raw0 = self.swept and self.blocks > 0
        held = C.alloc(self.held_count)
        inline = self.inline_sweep
        # Each cell's swap bit, copied by CX in the Euclid step, saves the Bezout tick its own
        # CCX for it and the CCX that cleared the B1 mark made from it, and halves the unsweep.
        swr = C.alloc(T - 2) if self.SWAP_WIRES and self.HELD and inline and not self.blocks and T >= 3 else None
        kept = C.alloc(sum(self.kept_plan()))
        labels = C.alloc(len(self.high_plan()))
        regs = [U, pool, held, kept, labels] + ([swr] if swr is not None else [])
        E = C.run(lambda U, pool, held, kept, labels, *sw: self._euclid_pass(
            U, pool, held, raw0, kept, sweep=not inline, swr=sw[0] if sw else None,
            labels=labels), regs)
        if raw0:
            self._code_block0(U, pool)
        lent = self.lent_after_pass
        G = self.coded_tape(pool)
        # Terminal V and the wires the code freed are zero; lend them during Bezout.
        C.free(pool[:lent])
        # U is (1,0,...,0). Reuse its entire word for the reverse coefficient;
        # the surviving V wire stays zero for the inverse Euclid boundary.
        C.x(U[0])
        S = U
        self._bezout_pass(G, Yr, S, raw0, swr)
        if swr is not None:
            # Bezout swept the tape it read; the inverse pass needs it unswept.
            C.apply(C.record(lambda G, W: self._b1_sweep_sw(G, W), [G, swr]), [G, swr], inverse=True)
        elif inline:
            # Bezout swept the tape it read; the inverse pass needs it unswept.
            C.apply(C.record(lambda G: self._b1_sweep(G), [G]), [G], inverse=True)
        if raw0:
            raw = [G[slot] for slot in self._free_slots(0, True)]
            for slot in self._free_slots(0, True):
                G[slot] = None
        # Restore the tape boundary using the newly allocated tail wire.
        pool[lent:] = [w for w in G if w is not None]
        C.x(Yr[0])
        pool[:lent] = C.alloc(lent)
        if raw0:
            # Block 0's free wires take the places of the last zero words, which are released.
            places = self._block0_places(Yr, pool)
            C.free([reg[i] for reg, i in places])
            for (reg, i), w in zip(places, raw):
                reg[i] = w
        C.apply(E, [Yr, pool, held, kept, labels] + ([swr] if swr is not None else []), inverse=True)
        C.free(held + kept + labels + (swr or []))
        load_const(C, Yr, q)
        Xr[:] = pool[:n]
        C.free(Yr + pool[n:])
        Yr[:] = S

    def apply(self, Xr, Yr, inverse=False):
        if self.block is None:
            self.block = self.C.record(self._body, [Xr, Yr])
        self.C.apply(self.block, [Xr, Yr], inverse)


def eq_const_into(C, ctrl, bits, value, target):
    """XOR equality into target, optionally controlled by ctrl."""
    flips = [w for i, w in enumerate(bits) if not value >> i & 1]
    for w in flips:
        C.x(w)
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
    # Square with HalfSquares where it applies: two groups of half squares, each built, taken
    # into X and unbuilt, in place of SquareFold's whole product and its Karatsuba junk.
    HALVES = True

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

    def _halves(self):
        return HalfSquares(self.mod) if self.HALVES and HalfSquares.applies(self.mod.q) else None

    def _square_offset(self):
        """K, where _mulsub leaves X - gen*(Y^2 - K): _kernel's constant add just before it
        takes gen*K off."""
        hs = self._halves()
        return hs.K if hs is not None else 0

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * (Y^2 - _square_offset()) mod q, Y untouched."""
        C, mod, n = self.C, self.mod, self.n
        hs = self._halves()
        if hs is not None:
            hs.mulsub(C, gen, Yr, Xr)
            return
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
            self._square_zeros = fold.zero_slots
            self._modadd_blk = C.record_shape(lambda k, a, b, z: fold.modadd(C, a, b, k[0], z),
                                              [1, n, n, len(fold.zero_slots)])
            self._square_junk = fold.junk
        J = C.alloc(self._square_junk)
        C.apply(self._square_blk, [Yr, J])
        # The square leaves these junk wires 0, so the modadd's ancillas can borrow them.
        C.apply(self._modadd_blk, [[gen], J[:n], Xr, [J[i] for i in self._square_zeros]],
                inverse=True)
        C.apply(self._square_blk, [Yr, J], inverse=True)
        C.free(J)

    def _kernel(self, c, gen, Xr, Yr, f0, f1):
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A

        def recode(inverse=False):
            if not inverse:
                C.cx(gen, c)
                C.cx(gen, f0)
            # gen ^= NOT c AND (f0 XOR f1): no c = 0 word has f0 = 0 and f1 = 1.
            C.cx(f1, f0)
            C.x(c)
            C.ccx(c, f0, gen)
            C.x(c)
            C.cx(f1, f0)
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
        start = len(C.prims) // 4
        # _mulsub leaves gen*_square_offset() too many on X; this add takes it off first, since
        # nothing between them reads X.
        mod.cmodadd_const(C, gen, 3 * a - self._square_offset(), Xr)
        _k1_cut(C, start, q, 3 * a - self._square_offset(), gen)
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
