"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, _maj, add_const, cswap, load_const
from .builder import Circuit
from .curve import Curve
from .square import SquareFold
from .two_bit import TwoBitMod
from .model import tick_bound, tick_widths
from . import block0_code, tail_code, tape_block
from .tape_block import FREE, TICKS


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for canonical 1 <= x < q; see model.py."""

    # The tape stays decoded: no block, block-zero or tail code, and so no pebbled ticks. Every
    # code call below is gated on blocks or tail, and each encode was paired with a decode that
    # nothing between reads past, so the passes are the same maps on the same states.
    RAW_TAPE = True
    # Carries the forward pass keeps for the inverse pass to clear, each saving one CCX per pass.
    KEEP = 1480

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


    def _shared_difference(self, U, V, g0, g1, top=None, kept=None):
        """Compare and form the even difference using one majority chain. With top set, bit w - 1
        is done by _top_bit on the domain that mode names. With kept, the chain's carries go onto
        those zero wires instead and stay there, which the inverse pass clears."""
        if kept is not None:
            return self._kept_difference(U, V, g0, g1, top, kept)
        C = self.C
        (cin,) = C.alloc(1)
        carries = [cin] + U[:-1]
        for wire in V:
            C.x(wire)
        chain = len(U) if top is None else len(U) - 1
        for i in range(chain):
            if i == 0:
                C.cx(U[i], V[i])
                C.cx(U[i], carries[i])
                C.cx(V[i], U[i])
            else:
                _maj(C, carries[i], V[i], U[i])
        if top is None:
            C.ccx(g0, U[-1], g1)
        else:
            self._top_bit(C, top, g0, U[-1], V[-1], U[-2], g1)
        for i in reversed(range(chain)):
            # Restore a_i and select S_i while unwinding the carry chain.
            if i == 0:
                C.cx(V[i], U[i])
                C.cx(g0, V[i])
            else:
                C.ccx(carries[i], V[i], U[i])
                C.ccx(g0, carries[i], V[i])
            C.cx(U[i], carries[i])
            C.cx(U[i], V[i])
            # S_i XOR c_i XOR 1 is a_i XOR v_i on the swap branch.
            C.cx(carries[i], V[i])
            C.x(V[i])
            if i != 0:
                C.ccx(g1, V[i], U[i])
            C.x(V[i])
            C.cx(carries[i], V[i])
        C.free([cin])
        for wire in V:
            C.x(wire)
        for wire in V[1:]:
            C.cx(g1, wire)

    def _kept_difference(self, U, V, g0, g1, top, kept):
        """_shared_difference with c_(i+1) built on kept[i - 1] by one CCX and never cleared, so
        U[i] keeps a_i and each bit's descent is the select and the exchange alone. Bit 0 is as
        there: U[0] holds c_1 = NOT g0 until its descent restores a_0 = 1. The returned carries
        are those of a + NOT(vhat) on the step's input, a function of that input."""
        C = self.C
        chain = len(U) if top is None else len(U) - 1
        assert len(kept) == chain - 1
        carry = [None, U[0]] + list(kept)
        (cin,) = C.alloc(1)
        for wire in V:
            C.x(wire)
        C.cx(U[0], V[0])
        C.cx(U[0], cin)
        C.cx(V[0], U[0])
        for i in range(1, chain):
            # c_(i+1) = a ^ (a ^ b)(a ^ c), with a, b and c restored.
            C.cx(U[i], V[i])
            C.cx(U[i], carry[i])
            C.ccx(carry[i], V[i], carry[i + 1])
            C.cx(U[i], carry[i + 1])
            C.cx(U[i], carry[i])
            C.cx(U[i], V[i])
        if top is None:
            C.ccx(g0, carry[chain], g1)
        else:
            self._top_bit(C, top, g0, U[-1], V[-1], carry[chain], g1)
        for i in reversed(range(1, chain)):
            # B[i] = b ^ g0*(a ^ c), S_i when g0 = 1; then a_i ^= g1 * NOT(S_i ^ c_i).
            C.cx(U[i], carry[i])
            C.ccx(g0, carry[i], V[i])
            C.cx(U[i], carry[i])
            C.cx(carry[i], V[i])
            C.x(V[i])
            C.ccx(g1, V[i], U[i])
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
        for wire in V[1:]:
            C.cx(g1, wire)

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
        """Per tick, how many carries it keeps: narrow ticks from the last back, within KEEP."""
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
                break
            plan[i], total = k, total + k
        return plan

    def _euclid_step(self, w: int, top_clear: bool = False, narrow: bool = False, last: bool = False,
                     keep: bool = False):
        """top_clear: the tick before ran at the same width, so V is below 2^(w - 1) on entry
        and the second strip's last swap, of two zero bits, never fires.
        narrow: the next width is w - 1, so by E7 this step leaves U below 2^(w - 1).
        last: the pass's last tick, whose input reaches (U, V, E) = (1, 0, 0) by E3."""
        if last:
            top_clear = narrow = keep = False
        key = (w, top_clear, narrow, last, keep)
        if key not in self._euclid_steps:
            C = self.C

            def step(U, V, E, cell, *kept):
                g0, b, held = cell
                offset = E[0]
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
                self._shared_difference(U[:w], V[:w], g0, b, top, kept[0] if keep else None)
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
                for i in range(w - 3 if top_clear else w - 2):
                    cswap(C, held, [V[i]], [V[i + 1]])
                # b = s ^ two, the cell's swept bit; then held = s, as two ^ (s ^ two).
                C.cx(held, b)
                C.cx(b, held)
                # A subtract keeps its swap s on E and holds the old offset e beside the cell.
                cswap(C, g0, [offset], [held])

            sizes = [w, w, 1, 3] + ([self._kept_count(w, top_clear, narrow, last)] if keep else [])
            self._euclid_steps[key] = C.record_shape(step, sizes)
        return self._euclid_steps[key]

    def _euclid_tick0(self):
        """First Euclid step on U=q, E=0, using unwritten tape for carries."""
        if self._euclid_tick0_block is None:
            C, q, n = self.C, self.mod.q, self.mod.n
            f1 = (1 << n) - q + 1
            k = n - 1 - ((f1 & -f1).bit_length() - 1)

            def step(U, V, E, cell, scr):
                C.free(scr)
                g0, b = cell
                C.cx(V[0], g0)
                for i in range(1, n):
                    if q >> i & 1:
                        C.x(V[i])
                    C.ccx(g0, V[i], U[i])
                    if q >> i & 1:
                        C.x(V[i])
                # The select left U[j] = q_j ^ g0*(v_j ^ q_j), so the first carry g0*v_j,
                # j the lowest set bit of f1, is U[j], or U[j] ^ g0 ^ 1 when q_j = 1.
                j = (f1 & -f1).bit_length() - 1

                def low(t):
                    C.cx(U[j], t)
                    if q >> j & 1:
                        C.cx(g0, t)
                        C.x(t)

                add_const(C, V, f1, ctrl=g0, low=low)
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
                C.cx(g0, b)

            self._euclid_tick0_block = C.record_shape(step, [n, n, 1, 2, k])
        return self._euclid_tick0_block

    def _euclid_pass(self, U: list[int], pool: list[int], held: list[int],
                     block0_raw: bool = False, kept: list[int] | None = None) -> None:
        # The pool starts as V plus clean external tape. It ends as the zero wires U did not
        # take back, then tape; U ends as its low wires and those zero wires. The boundary is a
        # permutation of U and the pool together. held, iters wires, starts zero and ends with
        # each subtract cell's old offset, zero elsewhere; only the inverse pass clears it.
        n = self.mod.n
        assert len(held) == self.iters
        plan = self.kept_plan() if kept is not None else [0] * self.iters
        assert kept is None or len(kept) == sum(plan)
        at = 0
        V, available, G = pool[:n], pool[n:], []
        lend = U[n - self.u_lent:]
        E = self.C.alloc(1)
        f1 = (1 << n) - self.mod.q + 1
        k = n - 1 - ((f1 & -f1).bit_length() - 1)
        for i in range(self.iters):
            w = self.widths[i]
            assert len(V) == w and len(available) >= 2
            cell, available = available[:2], available[2:]
            tick0 = self._euclid_tick0() if i == 0 and w == n and len(available) >= k else None
            if tick0 is not None and tick0.wires <= set(tick0.flat_in) and set(tick0.flat_out) == set(tick0.flat_in):
                self.C.apply(tick0, [U, V, E, cell, available[:k]])
            else:
                # A tick leaves V below 2^(w' - 1), w' its own width.
                top_clear = i >= 1 and w >= 3 and self.widths[i - 1] == w
                # E7: the next tick's U fits its width. E3: the last tick ends at (1, 0, 0).
                last = i == self.iters - 1
                narrow = not last and self.widths[i + 1] == w - 1
                cell = cell + held[i:i + 1]
                regs = [U[:w], V, E, cell]
                if plan[i]:
                    assert plan[i] == self._kept_count(w, top_clear, narrow, last)
                    regs.append(kept[at:at + plan[i]])
                    at += plan[i]
                self.C.apply(self._euclid_step(w, top_clear, narrow, last, plan[i] > 0), regs)
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
        late = swept and T >= 3
        if swept and T >= 2:
            (both,) = self.C.alloc(1)
            mark(T - 2 if late else T - 1)
        for i in range(T - 1, -1, -1):
            if swept and i + 1 < T and not (late and i == T - 2):
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
                clean.extend(self._encode(G, done // TICKS, swept))
        if both is not None:
            self.C.free([both])
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

    def _bezout(self):
        if self._bezout_step is None:
            C, mod = self.C, self.mod

            def step(g, R, S):
                self.narrow.step(C, g, R, S)

            self._bezout_step = C.record_shape(step, [2, mod.n, mod.n])
        return self._bezout_step

    def _bezout_tick(self, first, last, pebble=False):
        if (first, last, pebble) not in self._bezout_ticks:
            C, n = self.C, self.mod.n

            def tick(prev, g, nxt, R, S, flag, hand):
                # flag comes in holding [previous cell is B1], or 0 on a first tick, which has no
                # previous cell, and goes out holding [this cell is B1], or empty after the last.
                # hand carries a B1's top two bits to the next tick the same way, zero elsewhere.
                out = self.narrow.tick(C, prev or None, g, nxt[0] if nxt else None, R, S,
                                       pebble=pebble, held=None if first else flag[0],
                                       cell0=last, hand=hand)
                flag[:] = [] if out is None else [out]

            sizes = [0 if first else 2, 2, 0 if last else 1, n, n, 1, 2]
            self._bezout_ticks[first, last, pebble] = C.record_shape(tick, sizes)
        return self._bezout_ticks[first, last, pebble]

    def _b1_step(self, G, i, both):
        """Use the pass-wide mark while both cells are still raw."""
        cswap(self.C, both, [G[2 * i + 1]], [G[2 * i - 1]])

    def _bezout_pass(self, G, Yr, S, block0_raw=False):
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
            else:
                if not flag:
                    # The pass's first tick call: no tick before it handed a flag over.
                    flag = C.alloc(1)
                    if not first:
                        self.narrow._mark_b1(C, prev, flag[0])
                if not hand:
                    # No tick before this one handed a pair over, so its cell is not B1.
                    hand = C.alloc(2)
                C.apply(self._bezout_tick(first, last, pebble),
                        [prev, G[2 * i:2 * i + 2], nxt, Yr, S, flag, hand])
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
        held = C.alloc(T)
        kept = C.alloc(sum(self.kept_plan()))
        E = C.run(lambda U, pool, held, kept: self._euclid_pass(U, pool, held, raw0, kept),
                  [U, pool, held, kept])
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
        self._bezout_pass(G, Yr, S, raw0)
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
        C.apply(E, [Yr, pool, held, kept], inverse=True)
        C.free(held + kept)
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
