"""The four-symbol tick on a held offset, g = [g0, b].

    g0 = 1:         S = 2S + R mod q, then swap R and S on b
    g0 = 0, b = 0:  S = 2S mod q
    g0 = 0, b = 1:  S = 4S mod q

docs/mux-proof.md states and proves what follows. Credit: the offset
T = S + R between a doubling and the following subtract is the representation of section VI of
arXiv:2609.05625. The tick below keeps that offset as an unreduced integer, so nothing in it is
negative.

After the sweep in point_add.py, every B1 is followed by a subtract or ends the pass. A B1 then
adds R and does no reduction, leaving H = 2T + R, and the subtract after it makes 2H - R = 4T + R.

Add and shift now have complementary controls, so one rolling XOR wire D, the sum of the x lines
so far, lets one CCX carry both: a pair of neighbouring word bits takes D, or on J its own XOR.

H's bit 0 is R's, so it is not stored. On a B1 a pair of wires that the tick hands to the next
one holds H's two top bits between the ticks, and is 0 elsewhere. The only other wire held across them is
[this cell is B1], which the tick makes from its own sw with an X and two CX and the next tick
reads as held, in place of computing it with a CCX.
"""
from .arith import Mod, _maj, _selected_carry, cswap, gt_into


class Prefix:
    """AND(bits), and bits += carry, at 3 CCX per bit on about log2(len(bits)) clean wires.

    tops[j] holds a_j = AND(bits[:starts[j]]). When a_j is 1 every wire below starts[j] and
    every lower top reads 1, whether it holds its own value or a lower block's chain, so
    block j's chain can sit on those wires under an X and block lengths double. A chain
    value is only meaningful when a_j is 1, and reaches anything else only through a_j.
    Nothing may read bits between this and write, and write runs top down.
    """

    def __init__(self, C, bits, cap=None):
        self.C, self.bits = C, bits
        self.blocks = []
        self.released = None
        h = len(bits)
        self.top = bits[0] if h == 1 else None
        if h < 2:
            return
        (self.first,) = C.alloc(1)
        C.ccx(bits[0], bits[1], self.first)
        tops, start = [self.first], 2
        while start < h:
            carriers = bits[:start] + tops[:-1]
            m = min(len(carriers) + 1, h - start, cap or h)
            x, u = bits[start:start + m], carriers[:m - 1]
            for w in u:
                C.x(w)
            v = x[0]
            for i in range(1, m):
                C.ccx(v, x[i], u[i - 1])
                v = u[i - 1]
            (p,) = C.alloc(1)
            C.ccx(tops[-1], v, p)
            self.blocks.append((tops[-1], p, x, u, v))
            tops.append(p)
            start += m
        self.top = tops[-1]

    def release(self):
        """Clear the top below the last when no wire of the last block's chain holds it or the two
        wires that rebuild it. Only write reads it, and restore rebuilds it first."""
        if len(self.blocks) < 2:
            return
        a, p, _, _, v = self.blocks[-2]
        used = set(self.blocks[-1][3])
        if not {a, p, v} & used:
            self.C.ccx(a, v, p)
            self.C.free([p])
            self.released = len(self.blocks) - 2

    def restore(self):
        if self.released is None:
            return
        a, p, x, u, v = self.blocks[self.released]
        (q,) = self.C.alloc(1)
        self.C.ccx(a, v, q)
        self.blocks[self.released] = a, q, x, u, v
        last = self.blocks[-1]
        self.blocks[-1] = (q,) + last[1:]
        self.released = None

    # A block's carry is copied onto up to COPIES freed tops, one per COPY_EVERY of its bits past
    # the first, so its increments fan out.
    COPIES, COPY_EVERY = 3, 3

    def write(self, t):
        """bits += t, every wire restored or cleared."""
        C, bits = self.C, self.bits
        # Tops of the blocks already written, freed and so zero, take copies of later carries.
        pool = []
        assert self.released is None
        if not bits:
            return
        if len(bits) == 1:
            C.cx(t, bits[0])
            return
        blocks = self.blocks

        def carry_in(j):
            # p = a*v becomes a*t, the carry into block j.
            a, p, x, u, v = blocks[j]
            C.cx(t, v)
            C.ccx(a, v, p)
            C.cx(t, v)

        if blocks:
            carry_in(len(blocks) - 1)
        for j in reversed(range(len(blocks))):
            a, p, x, u, v = blocks[j]
            # The carry's copies are cleared by the same CX once its last increment is done; p
            # itself is not written meanwhile. The top freed last is left out: its clear is still
            # settling when this block starts.
            ctl = [p] + pool[:-1][:min(self.COPIES, (len(x) - 1) // self.COPY_EVERY)]
            for c in ctl[1:]:
                C.cx(p, c)
            # p takes the top positions, the copies the rest down to bit 1, top down.
            per = -(-(len(x) - 1) // len(ctl))

            def own(i):
                return ctl[(len(x) - 1 - i) // per]

            # x ^= p*below as x ^= p*NOT(below) before the uncompute, which reads below*x and
            # so is unchanged, then x ^= p once the block's chain is down: the chain holds only
            # the uncompute's CCX, and no later gate of the block reads x[i] or writes p.
            for i in reversed(range(1, len(x))):
                below = x[0] if i == 1 else u[i - 2]
                C.x(below)
                C.ccx(own(i), below, x[i])
                C.x(below)
                C.ccx(below, x[i], u[i - 1])
            for i in reversed(range(1, len(x))):
                C.cx(own(i), x[i])
            C.cx(p, x[0])
            for c in ctl[1:]:
                C.cx(p, c)
            for w in reversed(u):
                C.x(w)
            if j:
                # The block below takes its carry now, while its chain is still up, and p = a*t
                # clears as (a'*t)*v' with a = a'*v': the lower block's carry and chain top. So
                # the next block waits on this clear alone, not on a clear then a carry.
                carry_in(j - 1)
                _, p_below, _, _, v_below = blocks[j - 1]
                C.ccx(p_below, v_below, p)
            else:
                C.ccx(t, a, p)
            C.free([p])
            pool.append(p)
        C.x(bits[0])
        C.ccx(t, bits[0], bits[1])
        C.x(bits[0])
        C.ccx(bits[0], bits[1], self.first)
        C.cx(t, bits[1])
        C.cx(t, bits[0])
        C.free([self.first])


class TwoBitMod(Mod):
    # The tick's add-or-shift: the split add then a rotation, on SPLIT_COPIES copies of J, or
    # the one-pass add with the rolling XOR D.
    SPLIT_ADD = True
    SPLIT_COPIES = 4
    # With SPLIT_ADD, this many Bezout ticks of each pass, spread evenly, keep the one-pass add.
    SPLIT_SKIP = 0
    # The longest block of the tick's high-word Prefix P, or None for doubling blocks.
    P_CAP = 37

    def __init__(self, q):
        super().__init__(q)
        n = self.n
        self.f = f = (1 << n) - q
        self.corner = f ^ 2 * f ^ 3 * f
        self.window = (4 * f).bit_length()
        self.merged = 4 * f < 1 << n and self.corner <= 4 * f
        self.runs = self._runs() if self.merged else []

    def _runs(self):
        """Runs of two or more window bits above bit 0 on which every reduction adds the same bit:
        0 on a quiet run, 1 for every k >= 1 on a ones run (Joel's ca60e8a). As (start, end, ones),
        low to high."""
        f, w = self.f, self.window
        used = f | 2 * f | 3 * f | 4 * f | self.corner
        ones = f & 2 * f & 3 * f & 4 * f
        kind = [None] + ['ones' if ones >> i & 1 else None if used >> i & 1 else 'quiet'
                         for i in range(1, w)]
        runs, start = [], 1
        for i in range(2, w + 1):
            if i == w or kind[i] != kind[start]:
                if kind[start] and i - start >= 2:
                    runs.append((start, i, kind[start] == 'ones'))
                start = i
        return runs

    def _run_prefixes(self, C, L):
        """One Prefix per run. A ones run takes 2^h - 1 + c = 2^h - (1 - c), a decrement by 1 - c,
        which is an increment of the complement, so its bits sit under an X from here to the write."""
        prefixes = []
        for a, b, ones in self.runs:
            if ones:
                for x in L[a:b]:
                    C.x(x)
            prefixes.append(Prefix(C, L[a:b]))
        return prefixes

    def _bridge(self, C, selectors, masks, a, carry, top, out):
        """out ^= s ^ (carry ^ s)*top, the carry out of the run at a, whose addend bits all equal s.

        s is the XOR of the selectors set at a. On a quiet run s = 0. On a ones run the complement
        is incremented by 1 - carry, and the run carries out unless that increment wraps. When s = 0
        on a ones run the whole addend is 0, so carry = 0 and top is not read."""
        active = [sel for sel, mask in zip(selectors, masks) if mask >> a & 1]
        for sel in active:
            C.cx(sel, carry)
        C.ccx(carry, top, out)
        for sel in active:
            C.cx(sel, carry)
            C.cx(sel, out)

    def step(self, C, g, R, S):
        if self.merged:
            self.tick(C, None, g, None, R, S)
            return
        g0, b = g
        sw = self._begin(C, g0, b)
        self._cdbl(C, b, S)
        C.cx(sw, b)
        self.dbl_add(C, g0, R, S)
        cswap(C, sw, R, S)
        self._end(C, g0, b, sw)

    def _begin(self, C, g0, b):
        """b holds J and a fresh wire holds sw."""
        (sw,) = C.alloc(1)
        C.ccx(g0, b, sw)
        C.cx(sw, b)
        return sw

    def _end(self, C, g0, b, sw):
        """b holds b again; release sw."""
        C.ccx(g0, b, sw)
        C.free([sw])

    def _cdbl(self, C, bit, s):
        """s = (1 + bit)*s mod q. moddbl fixes zero, so park s when bit is 0."""
        parked = C.alloc(self.n)
        C.x(bit)
        cswap(C, bit, s, parked)
        C.x(bit)
        self.moddbl(C, s)
        C.x(bit)
        cswap(C, bit, s, parked)
        C.x(bit)
        C.free(parked)

    def _perches(self, selectors, masks):
        """Bits whose carry out can sit on their selector, as {bit: selector}."""
        perch, read = {}, set()
        for i in reversed(range(1, self.window)):
            on = [s for s, mask in zip(selectors, masks) if mask >> i & 1]
            free = [s for s in on if s not in read]
            if free and not any(a <= i < b for a, b, _ in self.runs):
                perch[i] = free[0], [s for s in on if s != free[0]]
            read.update(on)
        return perch

    def _chain_up(self, C, L, selectors, masks, M, perch, take, bottom=None):
        """M holds one Prefix per run. bottom, when given, is a wire that already holds the bottom
        carry chain[1]."""
        chain = {0: None}

        def rise(i):
            if i in perch:
                (s, rest), c = perch[i], chain[i]
                for other in rest:
                    C.cx(other, s)
                C.cx(s, L[i])
                C.cx(s, c)
                C.ccx(c, L[i], s)
                chain[i + 1] = s
            else:
                chain[i + 1] = take()
                _selected_carry(C, L[i], selectors, masks, i, chain[i], chain[i + 1])

        i = 0
        if bottom is not None:
            assert 0 not in perch
            chain[1], i = bottom, 1
        for (a, b, _), Mr in zip(self.runs, M):
            for j in range(i, a):
                rise(j)
            chain[b] = take()
            self._bridge(C, selectors, masks, a, chain[a], Mr.top, chain[b])
            i = b
        for j in range(i, self.window):
            rise(j)
        return chain

    def _chain_down(self, C, L, selectors, masks, M, perch, chain, drop, write, keep=0):
        """keep = 1 leaves the bottom carry chain[1] built, for the caller to clear."""

        def fall(i):
            if i in perch:
                (s, rest), c = perch[i], chain[i]
                C.ccx(c, L[i], s)
                C.cx(s, c)
                C.cx(s, L[i])
                for other in reversed(rest):
                    C.cx(other, s)
            else:
                _selected_carry(C, L[i], selectors, masks, i, chain[i], chain[i + 1])
                drop(chain[i + 1])
            if write and i:
                for selector, mask in zip(selectors, masks):
                    if mask >> i & 1:
                        C.cx(selector, L[i])
                C.cx(chain[i], L[i])

        i = self.window
        for (a, b, ones), Mr in reversed(list(zip(self.runs, M))):
            for j in reversed(range(b, i)):
                fall(j)
            self._bridge(C, selectors, masks, a, chain[a], Mr.top, chain[b])
            drop(chain[b])
            if write:
                # The run's sum is its bits + s*(2^h - 1) + carry, which is the Prefix's write of
                # carry, or for a ones run under s, of carry ^ s into the complement.
                active = [sel for sel, mask in zip(selectors, masks) if mask >> a & 1]
                for sel in active:
                    C.cx(sel, chain[a])
                Mr.restore()
                Mr.write(chain[a])
                for sel in active:
                    C.cx(sel, chain[a])
                if ones:
                    for x in L[a:b]:
                        C.x(x)
            i = a
        for j in reversed(range(keep, i)):
            fall(j)

    def _after_b1(self, C, prev):
        """Fresh wire = [the previous cell is B1], or None on a pass's first tick."""
        if prev is None:
            return None
        (p,) = C.alloc(1)
        self._mark_b1(C, prev, p)
        return p

    def _mark_b1(self, C, prev, p):
        for w in prev:
            C.x(w)
        C.ccx(prev[0], prev[1], p)
        for w in prev:
            C.x(w)

    def _drop_after_b1(self, C, prev, p):
        if prev is not None:
            self._mark_b1(C, prev, p)
            C.free([p])

    def _add_or_shift(self, C, J, a, word, held_lane=True, split=None, g0=None):
        """With g0, a B1 lane's sum bits n - 1 and n leave the word for two clean wires, returned
        as (hn, hb) and 0 on every other lane (Joel's 5f16404). A B1 lane, g0 = J = 0, is never
        held, so its word[n] is 0 on entry."""
        if self.SPLIT_ADD if split is None else split:
            return self._add_then_rotate(C, J, a, word, held_lane, g0)
        return self._add_or_shift_d(C, J, a, word, held_lane, g0)

    def _add_then_rotate(self, C, J, a, word, held_lane=True, g0=None):
        """The same map as _add_or_shift_d: off J an add of a >> 1 on its own, on J a rotation.

        The add is the ripple add with each sum bit made conditional in the up pass, while its carry
        wire is idle: with x the bit's carry line, u ^= NOT(J) AND NOT(x) leaves x*u, the only thing
        the down pass reads of u, unchanged, so the down pass is the plain undo, one CCX a bit, and
        two CX after it leave u = w ^ NOT(J)*x, the sum off J and w on J. On J the rotation of word's
        low n + 1 bits left is word = 2*word, as word[n] = 0 there; it is two layers of disjoint
        swaps, all n + 1 reversed and then all but word[0], under copies of J, each swap emitted
        once both its bits are final. Only one of the two acts on a lane, so their order is free.

        Gate order keeps the carry off CX in both directions: in the up pass each carry line's
        c ^= a is emitted before the CCX that sets the carry (they share only that target), so the
        carry is a CCX control and nothing else; in the down pass the c ^= a undo is emitted after
        the next lower CCX (again sharing only that target), so the carry line is a CCX target and
        nothing else. Run backwards, each pass keeps the same property.
        """
        n, k = self.n, self.SPLIT_COPIES
        (cin, carry, t) = C.alloc(3)
        copies = C.alloc(k)
        tree, sources = [], [J]
        while len(tree) < k:
            for s in list(sources):
                if len(tree) == k:
                    break
                c = copies[len(tree)]
                C.cx(s, c)
                tree.append((s, c))
                sources.append(c)
        C.x(J)
        ops, tgt = a[1:], word[:n - 1]
        m = n - 1
        carries = [cin] + ops[:-1]
        # Up pass. u = w ^ a, then every carry line's c ^= a, bottom up, each read before written.
        for i in range(m):
            C.cx(ops[i], tgt[i])
        for i in range(m):
            C.cx(ops[i], carries[i])
        for i in range(m):
            c, u = carries[i], tgt[i]
            C.ccx(c, u, ops[i])
            # x*(u ^ NOT J ^ NOT J*x) = x*u: u ^= NOT J AND NOT x.
            C.cx(J, u)
            C.cx(c, u)
            C.ccx(copies[i % k], c, u)
        # Off J, word[n - 1:] += the carry into bit n - 1, held lanes through hi.
        top, mid, hi = word[n - 1], word[n], word[n + 1]
        C.ccx(J, ops[-1], carry)
        C.ccx(carry, top, t)
        if held_lane:
            C.ccx(t, mid, hi)
        C.cx(t, mid)
        C.cx(carry, top)
        if g0 is None:
            C.x(top)
            C.ccx(carry, top, t)
            C.x(top)
        else:
            # J's wire holds NOT J, so g0 ^ NOT J is B1. hn = B1*top takes the final bit n - 1, and
            # with top cleared there, t's clear reads NOT top ^ B1, which is 0 on B1: t keeps the
            # carry into bit n, which is bit n itself, and mid gives it back.
            (hn,) = C.alloc(1)
            C.cx(J, g0)
            C.ccx(g0, top, hn)
            C.cx(hn, top)
            C.x(top)
            C.cx(g0, top)
            C.ccx(carry, top, t)
            C.cx(g0, top)
            C.x(top)
            C.cx(J, g0)
            C.cx(t, mid)
        C.ccx(J, ops[-1], carry)
        # On J, word[0..n] rotates left: layer one reverses all n + 1 bits, pairs (j, n - j), and
        # layer two all but word[0], pairs (j + 1, n - j). Each swap goes out once both its bits are
        # final, as the down pass finishes them top down.
        one, two, used = set(), set(), [0]

        def swap(x, y):
            cswap(C, copies[used[0] % k], [word[x]], [word[y]])
            used[0] += 1

        def ready(p):
            """Position p's layer-one swap is out, or p is the middle bit, which has none."""
            return 2 * p == n or min(p, n - p) in one

        # Down pass. After CCX i, ops[i] holds a_i ^ a_(i+1) (a_i at the top) until the line above
        # is final and lends its a_(i+1).
        def finish(i):
            o, u = ops[i], tgt[i]
            if i + 1 < m:
                C.cx(ops[i + 1], o)
            C.cx(J, u)
            C.cx(o, u)
            if i < (n + 1) // 2:
                swap(i, n - i)
                one.add(i)
            for j in (i, i - 1):
                if 0 <= j < n // 2 and j not in two and ready(j + 1) and ready(n - j):
                    swap(j + 1, n - j)
                    two.add(j)

        for i in reversed(range(m)):
            C.ccx(carries[i], tgt[i], ops[i])
            if i + 1 < m:
                finish(i + 1)
        finish(0)
        C.cx(ops[0], cin)
        assert len(one) == (n + 1) // 2 and len(two) == n // 2
        C.x(J)
        for s, c in reversed(tree):
            C.cx(s, c)
        C.free(copies)
        if g0 is None:
            C.free([cin, carry, t])
            return None
        C.free([cin, carry])
        return hn, t

    def _add_or_shift_d(self, C, J, a, word, held_lane=True, g0=None):
        """On J, word = 2*word. Off J, word += a >> 1, through word's top two bits.

        word[n] must be 0 on J: a J cell never follows a B1 (D1), so a J lane is never held.
        held_lane=False promises that no lane is held, so word[n] = 0 on every lane."""
        n = self.n
        (cin, D, t) = C.alloc(3)
        ops, tgt = a[1:], word[:n - 1]
        carries = [cin] + ops[:-1]
        for i in range(n - 1):
            _maj(C, carries[i], tgt[i], ops[i])
            C.cx(carries[i], D)
        c, top, mid, hi = ops[-1], word[n - 1], word[n], word[n + 1]
        # t = c * top is the carry out of bit n - 1 off J. On J it is junk, but there mid = 0, so
        # neither hi's update nor mid's needs a J gate on t.
        C.ccx(c, top, t)
        if held_lane:
            C.ccx(t, mid, hi)
        # mid ^= t off J; on J, mid = 0 takes top, J's top move.
        C.cx(t, top)
        C.ccx(J, top, mid)
        C.cx(t, top)
        C.cx(t, mid)
        if g0 is None:
            C.ccx(c, top, t)
        else:
            # hn = B1*(top ^ c), with B1 = 1 ^ g0 ^ J on the g0 wire, and top ^= hn leaves c there.
            # t's clear then reads top ^ B1: c*top off B1, and c*(c ^ 1) = 0 on B1, where t keeps
            # the carry c*top, bit n.
            (hn,) = C.alloc(1)
            C.x(g0)
            C.cx(J, g0)
            C.cx(c, top)
            C.ccx(g0, top, hn)
            C.cx(c, top)
            C.cx(hn, top)
            C.cx(g0, top)
            C.ccx(c, top, t)
            C.cx(g0, top)
            C.cx(J, g0)
            C.x(g0)
        # top ^= c ^ D off J: its carry, and the cancel of the D that pair n - 2 hands it. On J top
        # becomes 0, read through mid, which holds top's old value there.
        C.cx(c, D)
        C.cx(D, top)
        C.cx(mid, D)
        C.ccx(J, D, top)
        C.cx(mid, D)
        C.cx(c, D)
        if g0 is not None:
            C.cx(t, mid)
        for i in reversed(range(n - 1)):
            c, u, v, o = carries[i], word[i], word[i + 1], ops[i]
            C.ccx(c, u, o)
            C.cx(o, u)
            # u and v both take D, or on J their own XOR. Then c still holds x_i, for D.
            C.cx(v, u)
            C.cx(D, v)
            C.cx(u, D)
            C.ccx(J, D, v)
            C.cx(u, D)
            C.cx(v, u)
            C.cx(c, D)
            C.cx(o, c)
        if g0 is None:
            C.free([cin, D, t])
            return None
        C.free([cin, D])
        return hn, t

    def _sweep_in(self, C, g0, b, sw, up, nxt, tau):
        """The B1 sweep's step into this cell, on an unswept tape read top down.

        Sweep step k swaps b_k and b_(k-1) when g0_k = g0_(k-1) = 0. On the tape no B1 has a
        no-subtract cell above it in Euclid order, so that b_(k-1) is 1, and the step toggles both
        bits by t_k = NOT g0_k AND NOT g0_(k-1) AND NOT x_k, with x_k = b_k XOR t_(k+1) this cell's
        bit after the step above. tau holds t_(k+1) on entry, and b holds x_k, sw = g0_k AND x_k.
        t_(k+1) = NOT g0_(k+1) AND B1_k with B1_k = NOT g0_k AND NOT x_k, since x_k = 0 when t_(k+1)
        is 1 and b_k = 1 below a no-subtract cell, so t_k = t_(k+1) XOR B1_k AND (g0_(k-1) XOR
        g0_(k+1)), one CCX. up is g0_(k+1), or None on a pass's first tick, where tau is 0 and
        g0_(k+1) is taken as 1. On return tau holds t_k and b the swept bit x_k XOR t_k; sw is
        unchanged, as t_k implies NOT g0_k."""
        C.x(sw)
        C.cx(g0, sw)
        C.cx(b, sw)
        if up is None:
            C.x(nxt)
        else:
            C.cx(up, nxt)
        C.ccx(sw, nxt, tau)
        if up is None:
            C.x(nxt)
        else:
            C.cx(up, nxt)
        C.cx(b, sw)
        C.cx(g0, sw)
        C.x(sw)
        C.cx(tau, b)

    def _last_tick(self, C, prev, g, a, s, held, hand, tau):
        """Cell 0's tick, where g0 implies b, with no add and no swap (Joel's d21aeb3). The pass
        ends with R = 0, so on entry a J or B1 lane has R = 0 and a subtract lane has
        2T + R = 0 mod q, or 2H - R = 0 mod q when held. The result on a subtract lane is then R
        itself, made from S alone:

            plain:  S = -2T mod q = ~(V + (k + 1)f - 1) for V = 2T, k = floor(V / q)
            held:   S = 2H mod q
            J:      S = 4T mod q        B1:  S = 2T mod q

        after which R ^= g0*S clears R in place of the swap. V + kf = r + e*2^n for r = V mod q,
        so the plain lane's sum overflows exactly when e is set, as the other lanes' does, and
        its complement is q - r. That is R, except on the lane T = R = 0, where it is q."""
        n, f, w = self.n, self.f, self.window
        g0, b = g
        if hand is None:
            assert prev is None
            (spill, high) = C.alloc(2)
        else:
            (spill, high) = hand
        word = s + [spill, high]
        if prev is None:
            p = None
        elif held is not None:
            p = held
        else:
            p = self._after_b1(C, prev)
        up = prev[0] if prev else None
        if tau is not None:
            C.cx(tau, b)
            if up is not None:
                # t_1 = NOT g0_1 AND B1_0, and B1_0 = NOT x_0 because g0_0 implies x_0. t_0 = 0.
                C.x(b)
                C.x(up)
                C.ccx(b, up, tau)
                C.x(up)
                C.x(b)
        # b holds J = b ^ g0 from here to the end.
        C.cx(g0, b)
        for i in reversed(range(n)):
            cswap(C, b, [word[i]], [word[i + 1]])
        if p is not None:
            # A held word is H - R0, with bit 0 clear.
            C.ccx(p, a[0], word[0])
        (l0,) = C.alloc(1)
        low = [l0] + word[:n - 1]
        L, H = low[:w], low[w:]
        k0, k1, k2 = word[n - 1:n + 2]
        (one,) = C.alloc(1)
        C.x(one)
        P, M = Prefix(C, H, cap=self.P_CAP), self._run_prefixes(C, L)
        spare = H[:] if len(H) > 1 else []

        def alloc():
            (c,) = C.alloc(1)
            return c

        def take():
            if not spare:
                return alloc()
            c = spare.pop()
            C.x(c)
            return c

        def drop(c):
            if c in H:
                C.x(c)
                spare.append(c)
            else:
                C.free([c])

        def plain():
            """Toggle the g0 wire between g0 and [plain subtract lane] = g0 ^ p."""
            if p is not None:
                C.cx(p, g0)

        (three,) = C.alloc(1)
        C.ccx(k0, k1, three)
        selectors = [one, k0, k1, three]
        masks = [f, f ^ 2 * f, f ^ 3 * f, self.corner ^ 4 * f]
        chain = self._chain_up(C, L, selectors, masks, M, {}, take)
        (e,) = C.alloc(1)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        self._chain_down(C, L, selectors, masks, M, {}, chain, drop, False)
        C.x(one)
        C.ccx(three, e, k2)
        C.ccx(e, k1, three)
        C.free([one])
        one = three
        C.ccx(k0, e, k1)
        C.cx(e, k0)
        # A plain lane has k <= 1 and adds (k + 1)f - 1, which is f ^ 1 or 2f ^ 3.
        plain()
        (ck0,) = C.alloc(1)
        C.ccx(g0, k0, ck0)
        selectors = [one, k0, k1, k2, g0, ck0]
        masks = [self.corner, f, 2 * f, 4 * f, f ^ 1, 2 * f ^ 2]
        perch = self._perches(selectors, masks)
        chain = self._chain_up(C, L, selectors, masks, M, perch, alloc)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        C.free([e])
        P.write(chain[w])
        self._chain_down(C, L, selectors, masks, M, perch, chain, drop, True)
        C.ccx(g0, k0, ck0)
        C.free([ck0])
        for x in word[:n]:
            C.cx(g0, x)
        plain()

        if p is not None:
            cswap(C, p, [k1], [k2])
        # one = k0*k1 becomes z = (p ^ wire k2)*~k0, which is [k is 0 or 4] on a held lane and 0
        # elsewhere, as in tick with p for g0.
        z = one
        if p is not None:
            C.cx(p, k1)
        C.ccx(k0, k1, z)
        if p is not None:
            C.cx(p, k1)
            C.cx(p, z)
        C.cx(k2, z)

        for r, x in zip(a, [k0] + word[:n - 1]):
            C.ccx(g0, x, r)
        # R is now 0, or q on the lane T = R = 0, where S is q too. R's top bit flags that lane.
        q = self.q
        for i, (r, x) in enumerate(zip(a, [k0] + word[:n - 1])):
            C.cx(r, x)
            if i < n - 1 and q >> i & 1:
                C.cx(a[n - 1], r)
        plain()
        C.ccx(g0, a[n - 1], k1)
        C.ccx(g0, k1, a[n - 1])
        # Wire k1 holds the flag on a plain lane and k's bit 2 on a held lane. With R = 0,
        # [S >> 1 > R >> 1] is [S >= 2]: it clears bit 2 under z, and the flag is
        # [S is 0] = ~S0 * ~[S >= 2] on a plain lane.
        (zc,) = C.alloc(1)
        C.x(k0)
        C.ccx(g0, k0, zc)
        C.cx(zc, k1)
        C.cx(zc, z)
        gt_into(C, low[1:], a[1:], k1, ctrl=z)
        C.cx(zc, z)
        C.ccx(g0, k0, zc)
        C.x(k0)
        C.free([zc])
        plain()
        if p is not None:
            C.cx(k2, p)
            C.x(k0)
            C.ccx(p, k0, z)
            C.x(k0)
            C.cx(k2, p)
            # 2H = S + kq and H's bit 0 is S's, so k = f*(S - 2*S0) mod 4.
            if not f >> 1 & 1:
                C.cx(k0, low[1])
            C.ccx(p, low[1], k2)
            if not f >> 1 & 1:
                C.cx(k0, low[1])
        C.free([z])
        self._drop_after_b1(C, prev, p)
        C.ccx(b, low[1], k1)
        if f >> 1 & 1:
            C.ccx(b, k0, k1)
        C.cx(g0, b)
        C.free([l0, spill, high])
        if hand is not None:
            hand[:] = []
        s[:] = s[-1:] + s[:-1]
        return None

    def tick(self, C, prev, g, nxt, a, s, pebble=False, held=None, cell0=False, hand=None,
             tau=None, sw=None, tape=False, split=None, term=None):
        """One tick of a swept tape. prev is the cell before, or None on the first tick.
        split picks the add-or-shift (see _add_or_shift); None takes SPLIT_ADD.
        nxt is the next cell's g0 wire, or None on the last tick, whose B1 reduces normally.

        held is the wire the tick before handed over, holding [prev is B1], or None to compute
        it here. cell0 says this is cell 0, where g0 implies b (Euclid's tick 0 runs on U = q > x,
        so its subtract always swaps). Off the tape it is _last_tick; on it, sw is g0 itself and J
        is b ^ g0.
        Returns the wire that holds [this cell is B1] for the next tick, or None when nxt is None.

        hand is the pair of wires the tick before handed over, holding H's bits n and n + 1 when
        the previous cell is B1 and 0 otherwise; this tick uses them as word's top two wires. On
        return, hand holds this tick's pair for the next tick, or is empty when nxt is None.
        hand may be None only when there is neither a previous nor a next cell.
        sw and [prev is B1] each stay on one wire for the whole tick.

        tau, when given, carries the B1 sweep through an unswept tape: it holds the toggle the
        sweep step above leaves on this cell, and on return the one this cell's step leaves on the
        next, or 0 after cell 0. See _sweep_in.

        With tape, each cell's swap bit s = g0 AND b sits on the tape's swap wire: sw is this
        cell's, read in place of a computed sw, and the returned [this cell is B1] is formed on it
        in place. held is then the previous cell's swap wire holding its mark, which this tick
        turns back into that cell's s with X and CX instead of clearing it with a CCX.

        term, when given, is the terminal cell's g0 on the tick after it, which has no previous
        cell and enters with S = term*R."""
        if cell0 and not tape:
            return self._last_tick(C, prev, g, a, s, held, hand, tau)
        n, f, w = self.n, self.f, self.window
        tick_start = len(C.prims)
        g0, b = g
        if hand is None:
            assert prev is None and nxt is None
            (spill, high) = C.alloc(2)
        else:
            (spill, high) = hand
        word = s + [spill, high]
        if term is not None:
            assert prev is None
            (rtop,) = C.alloc(1)
            C.cx(s[n - 1], rtop)
        if prev is None:
            p = None
        elif held is not None:
            p = held
        else:
            p = self._after_b1(C, prev)
        up = prev[0] if prev else None
        if tau is not None:
            C.cx(tau, b)
        if cell0:
            if tau is not None and up is not None:
                # t_1 = NOT g0_1 AND B1_0, and B1_0 = NOT x_0 because g0_0 implies x_0. t_0 = 0.
                C.x(b)
                C.x(up)
                C.ccx(b, up, tau)
                C.x(up)
                C.x(b)
            sw = g0
            C.cx(g0, b)
        elif tau is not None:
            if not tape:
                (sw,) = C.alloc(1)
                C.ccx(g0, b, sw)
            self._sweep_in(C, g0, b, sw, up, nxt, tau)
            C.cx(sw, b)
        else:
            sw = self._begin(C, g0, b)
        if p is not None:
            for x in word:
                C.cx(p, x)
        ahead = self._add_or_shift(C, b, a, word, held_lane=prev is not None, split=split,
                                   g0=None if nxt is None else g0)
        if p is not None:
            for x in word:
                C.cx(p, x)
        (l0,) = C.alloc(1)
        C.ccx(g0, a[0], l0)
        low = [l0] + word[:n - 1]
        L, H = low[:w], low[w:]
        k0, k1, k2 = word[n - 1:n + 2]

        # one = 0 only on a B1, which is the tick's inert lane. b holds J here, so this needs no CCX.
        (one,) = C.alloc(1)
        if nxt is None:
            C.x(one)
        else:
            C.cx(g0, one)
            C.cx(b, one)
            # On a B1, the add moved H's bits n and n + 1 to ahead, the next tick's word top.
        P, M = Prefix(C, H, cap=self.P_CAP), self._run_prefixes(C, L)

        # With at most one H bit, P reads a data wire or nothing, and there is nothing to borrow.
        spare = H[:] if len(H) > 1 else []

        def alloc():
            (c,) = C.alloc(1)
            return c

        def take():
            if not spare:
                return alloc()
            c = spare.pop()
            C.x(c)
            return c

        def drop(c):
            if c in H:
                C.x(c)
                spare.append(c)
            else:
                C.free([c])

        (three,) = C.alloc(1)
        C.ccx(k0, k1, three)
        selectors = [one, k0, k1, three]
        masks = [f, f ^ 2 * f, f ^ 3 * f, self.corner ^ 4 * f]
        # The first chain's bottom carry, l0*(one ^ K0), sits on a clean wire, not a borrowed H
        # wire, so that it is a true carry on every lane and can become the second chain's.
        (c1,) = C.alloc(1)
        _selected_carry(C, L[0], selectors, masks, 0, None, c1)
        chain = self._chain_up(C, L, selectors, masks, M, {}, take, bottom=c1)
        (e,) = C.alloc(1)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        self._chain_down(C, L, selectors, masks, M, {}, chain, drop, False, keep=1)
        if nxt is None:
            C.x(one)
        else:
            C.cx(b, one)
            C.cx(g0, one)
        C.cx(sw, b)
        C.ccx(three, e, k2)
        # three = K0*K1 becomes k0*k1 = K0*K1 ^ e*K1 for k = K + e, so it reads k1 before k1 changes.
        C.ccx(e, k1, three)
        C.free([one])
        one = three
        C.ccx(k0, e, k1)
        C.cx(e, k0)
        # c1 = l0*(one ^ K0) becomes the second chain's l0*k0 = l0*(K0 ^ e). l0 = 1 only off B1,
        # where one = 1, so the change is l0*(1 ^ e).
        C.cx(l0, c1)
        C.ccx(l0, e, c1)
        selectors, masks = [one, k0, k1, k2], [self.corner, f, 2 * f, 4 * f]
        perch = self._perches(selectors, masks)
        # Every low carry stays on its own wire until its fall, and e stays on its own wire.
        chain = self._chain_up(C, L, selectors, masks, M, perch, alloc, bottom=c1)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        C.free([e])
        P.write(chain[w])
        # The bottom carry c1 = l0*k0 stays: on a held lane l0 = R0, so it is the k0*R0 that bit one reads.
        self._chain_down(C, L, selectors, masks, M, perch, chain, drop, True, keep=1)
        c1 = chain[1]

        # On a held lane k1 and k2 change wires, so the compare's target is k's top bit there.
        if p is not None:
            cswap(C, p, [k1], [k2])
        def mark_even(z):
            C.cx(k2, g0)
            C.x(k0)
            C.ccx(g0, k0, z)
            C.x(k0)
            C.cx(k2, g0)

        compare_start = len(C.prims)
        # one = k0*k1 becomes z = (g0 ^ wire k2)*~k0 with one CCX (Joel's a43c00c). Wire k2 holds
        # k's bit 1 on a held lane and k's bit 2 elsewhere, and k <= 4 gives k0*k2 = 0, so
        # k0*k1 = k0*(wire k1 ^ wire k2) and z ^ k0*k1 = k0*(wire k1 ^ g0) ^ g0 ^ wire k2.
        z = one
        C.cx(g0, k1)
        C.ccx(k0, k1, z)
        C.cx(g0, k1)
        C.cx(g0, z)
        C.cx(k2, z)
        if term is None:
            gt_into(C, a[1:], low[1:], k1, ctrl=z)
        else:
            # Joel's bb0fe4e. V = (2*term + 1)*R here, so an even k is 2 only when V = 3R >= 2q.
            # Then R >= 2q/3, which is at least 2^(n - 1) because 4f < 2^n, and k = 0 has
            # R < q/3 < 2^(n - 1): k's bit 1 is S's top bit, copied before the add.
            C.ccx(z, rtop, k1)
            C.ccx(term, a[n - 1], rtop)
            C.free([rtop])
        mark_even(z)
        C.free([z])
        if p is not None:
            # S = R + k*f mod 4 names k's bit 1, which now sits on k2. c1 holds its k0*R0 term.
            def bit_one(m):
                C.cx(low[1], m)
                C.cx(a[1], m)
                if f >> 1 & 1:
                    C.cx(k0, m)

            bit_one(c1)
            C.ccx(p, c1, k2)
            bit_one(c1)
        if tape and p is not None:
            # p = 1 ^ g0 ^ b ^ s, formed by the tick before from this same cell's g0 and swept b.
            C.x(p)
            C.cx(prev[0], p)
            C.cx(prev[1], p)
        else:
            self._drop_after_b1(C, prev, p)
        C.ccx(l0, k0, c1)
        C.free([c1])
        C.cx(sw, b)
        C.ccx(b, low[1], k1)
        if f >> 1 & 1:
            C.ccx(b, k0, k1)
        C.cx(l0, k0)
        C.cx(sw, b)
        C.ccx(g0, a[0], l0)
        C.free([l0, spill, high])
        if hand is not None:
            hand[:] = [] if nxt is None else ahead
        s[:] = s[-1:] + s[:-1]
        # The swap's control goes by a star of CX onto free wires: each is zero, as every wire is
        # once freed, and the reversed star clears it. Only wires this tick has touched, so its
        # arity and peak stay, and none touched since the compare began, which would hold the
        # copies back until the compare is done, and in the inverse hold the compare back.
        touched = set()
        for k in (1, 2, 3):
            touched.update(C.prims[tick_start + k::4])
        late = set()
        for k in (1, 2, 3):
            late.update(C.prims[compare_start + k::4])
        copies = sorted(touched.intersection(C.free_stack) - late)
        for c in copies:
            C.cx(sw, c)
        ctl = copies or [sw]
        for j, (u, v) in enumerate(zip(a, s)):
            cswap(C, ctl[j * len(ctl) // n], [u], [v])
        for c in reversed(copies):
            C.cx(sw, c)
        if cell0:
            return None
        if nxt is None:
            self._end(C, g0, b, sw)
            return None
        # [this cell is B1] = (1 ^ g0)(1 ^ b) = 1 ^ g0 ^ b ^ sw, the next tick's held.
        C.x(sw)
        C.cx(g0, sw)
        C.cx(b, sw)
        return sw
