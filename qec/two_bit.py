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
from functools import lru_cache
import heapq

from .arith import Mod, _maj, _selected_carry, cswap, gt_into


@lru_cache(maxsize=None)
def _cap4_schedule(n, top):
    """The handover's reversible low-carry pebble search, per chain length."""
    if not 1 <= n <= 15:
        return None
    dist, back, heap, best = {0: 0}, {0: None}, [(0, 0)], None
    while heap:
        cost, state = heapq.heappop(heap)
        if cost != dist[state]:
            continue
        live = state.bit_count()
        if state >> (n - 1) & 1 and live + top <= 4:
            candidate = (cost - n + n - live, state)
            if best is None or candidate < best:
                best = candidate
        for i in range(1, n + 1):
            if i == 1 or state >> (i - 2) & 1:
                nxt = state ^ (1 << (i - 1))
                if nxt.bit_count() <= 4 and cost + 1 < dist.get(nxt, 1 << 30):
                    dist[nxt], back[nxt] = cost + 1, (state, i)
                    heapq.heappush(heap, (cost + 1, nxt))
    if best is None:
        return None
    state = best[1]
    kept = frozenset(i + 1 for i in range(n) if state >> i & 1)
    moves = []
    while back[state] is not None:
        state, i = back[state]
        moves.append(i)
    return tuple(reversed(moves)), kept


class Prefix:
    """AND(bits), and bits += carry, at 3 CCX per bit on about log2(len(bits)) clean wires.

    tops[j] holds a_j = AND(bits[:starts[j]]). When a_j is 1 every wire below starts[j] and
    every lower top reads 1, whether it holds its own value or a lower block's chain, so
    block j's chain can sit on those wires under an X and block lengths double. A chain
    value is only meaningful when a_j is 1, and reaches anything else only through a_j.
    Nothing may read bits between this and write, and write runs top down.
    """

    def __init__(self, C, bits):
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
            m = min(len(carriers) + 1, h - start)
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

    def write(self, t):
        """bits += t, every wire restored or cleared."""
        C, bits = self.C, self.bits
        assert self.released is None
        if not bits:
            return
        if len(bits) == 1:
            C.cx(t, bits[0])
            return
        for a, p, x, u, v in reversed(self.blocks):
            # p = a*v becomes a*t, the carry into this block.
            C.cx(t, v)
            C.ccx(a, v, p)
            C.cx(t, v)
            for i in reversed(range(1, len(x))):
                below = x[0] if i == 1 else u[i - 2]
                C.ccx(below, x[i], u[i - 1])
                C.ccx(p, below, x[i])
            C.cx(p, x[0])
            for w in reversed(u):
                C.x(w)
            C.ccx(t, a, p)
            C.free([p])
        C.ccx(bits[0], bits[1], self.first)
        C.ccx(t, bits[0], bits[1])
        C.cx(t, bits[0])
        C.free([self.first])


class TwoBitMod(Mod):
    def __init__(self, q):
        super().__init__(q)
        n = self.n
        self.f = f = (1 << n) - q
        self.corner = f ^ 2 * f ^ 3 * f
        self.window = (4 * f).bit_length()
        self.merged = 4 * f < 1 << n and self.corner <= 4 * f
        self.run = self._quiet_run() if self.merged else (self.window, self.window)

    def _quiet_run(self):
        """Longest run of window bits above bit 0 with no mask set, or an empty one."""
        f, w = self.f, self.window
        used = f | 2 * f | 3 * f | 4 * f | self.corner
        best, start = (w, w), None
        for i in range(1, w + 1):
            quiet = i < w and not used >> i & 1
            if quiet and start is None:
                start = i
            if not quiet and start is not None:
                if i - start >= 2 and i - start > best[1] - best[0]:
                    best = (start, i)
                start = None
        return best

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
        a, b = self.run
        perch, read = {}, set()
        for i in reversed(range(1, self.window)):
            if a <= i < b:
                continue
            on = [s for s, mask in zip(selectors, masks) if mask >> i & 1]
            free = [s for s in on if s not in read]
            if free:
                perch[i] = free[0], [s for s in on if s != free[0]]
            read.update(on)
        return perch

    def _chain_up(self, C, L, selectors, masks, M, perch, take, kept=None, moves=None,
                  bottom=None):
        """kept, when given, names the low carries that stay. Each other one is cleared once the
        next kept carry is built, so only kept carries are live above it; _chain_down rebuilds them.
        bottom, when given, is a wire that already holds the bottom carry chain[1]."""
        a, b, w = *self.run, self.window
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

        start = 0
        if bottom is not None:
            assert kept is None and moves is None and a >= 1 and 0 not in perch
            chain[1], start = bottom, 1
        if moves is not None:
            for j in moves:
                if j in chain:
                    _selected_carry(C, L[j - 1], selectors, masks,
                                    j - 1, chain[j - 1], chain[j])
                    C.free([chain.pop(j)])
                else:
                    rise(j - 1)
            start = a - 1
        last = start
        for i in range(start, a):
            rise(i)
            if kept is not None and (i + 1 in kept or i + 1 == a):
                for j in range(i, last, -1):
                    if j not in kept:
                        _selected_carry(C, L[j - 1], selectors, masks,
                                    j - 1, chain[j - 1], chain[j])
                        C.free([chain.pop(j)])
                last = i + 1
        if b > a:
            chain[b] = take()
            C.ccx(chain[a], M.top, chain[b])
        for i in range(b, w):
            rise(i)
        return chain

    def _chain_down(self, C, L, selectors, masks, M, perch, chain, drop, write, take=None, keep=0):
        """keep = 1 leaves the bottom carry chain[1] built, for the caller to clear."""
        a, b, w = *self.run, self.window

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

        for i in reversed(range(b, w)):
            fall(i)
        if b > a:
            C.ccx(chain[a], M.top, chain[b])
            drop(chain[b])
            if write:
                M.restore()
                M.write(chain[a])
        for i in reversed(range(keep, a)):
            if i not in chain:
                start = max(j for j in chain if j < i)
                for j in range(start + 1, i + 1):
                    chain[j] = take()
                    _selected_carry(C, L[j - 1], selectors, masks,
                                    j - 1, chain[j - 1], chain[j])
            fall(i)

    def _pebble_schedule(self, perch):
        a, b = self.run
        # The checked construction pebbles the unperched prefix, then uses
        # the selector at a-1 as carry a. Other shapes retain today's schedule.
        if sorted(i for i in perch if i < a) != [a - 1]:
            return None
        return _cap4_schedule(a - 1, int(b > a))

    def _kept_carries(self, perch):
        """The low carries 1..a that the second chain keeps, chosen so that the most clean low wires
        live at once, while building a segment or at the chain's top, is least; then fewest dropped.

        A carry out of a perched bit is on a selector, so it costs no wire and is always kept, and so
        is the carry into a perched bit, which the perch changes until its fall. Each dropped carry
        costs 2 CCX: one to clear it going up, one to rebuild it going down.
        """
        a, b = self.run
        top = int(b > a)
        fixed = {i + 1 for i in perch if i < a} | {i for i in perch if 0 < i < a} | {a}
        free = [i for i in range(1, a) if i not in fixed]
        if len(free) > 16:
            return set(range(a + 1))
        on_wire = lambda i: i - 1 not in perch
        best = None
        for mask in range(1 << len(free)):
            kept = fixed | {i for k, i in enumerate(free) if mask >> k & 1}
            live = peak = last = 0
            for i in sorted(kept):
                peak = max(peak, live + sum(on_wire(j) for j in range(last + 1, i + 1)))
                live += on_wire(i)
                last = i
            key = (max(peak, live + top), a - len(kept))
            if best is None or key < best[0]:
                best = key, kept
        return best[1]

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

    def _add_or_shift(self, C, J, a, word, held_lane=True):
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
        C.ccx(c, top, t)
        # top ^= c ^ D off J: its carry, and the cancel of the D that pair n - 2 hands it. On J top
        # becomes 0, read through mid, which holds top's old value there.
        C.cx(c, D)
        C.cx(D, top)
        C.cx(mid, D)
        C.ccx(J, D, top)
        C.cx(mid, D)
        C.cx(c, D)
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
        C.free([cin, D, t])

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

    def tick(self, C, prev, g, nxt, a, s, pebble=False, held=None, cell0=False, hand=None,
             tau=None, sw=None, tape=False):
        """One tick of a swept tape. prev is the cell before, or None on the first tick.
        nxt is the next cell's g0 wire, or None on the last tick, whose B1 reduces normally.

        held is the wire the tick before handed over, holding [prev is B1], or None to compute
        it here. cell0 says this is cell 0, where g0 implies b (Euclid's tick 0 runs on U = q > x,
        so its subtract always swaps): there sw is g0 itself and J is b ^ g0.
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
        turns back into that cell's s with X and CX instead of clearing it with a CCX."""
        n, f, w = self.n, self.f, self.window
        g0, b = g
        if hand is None:
            assert prev is None and nxt is None
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
        self._add_or_shift(C, b, a, word, held_lane=prev is not None)
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
            C.x(one)
            # On a B1, H's bits n and n + 1 move to two fresh wires, the next tick's word top.
            ahead = C.alloc(2)
            for k, y in zip((k0, k1), ahead):
                C.ccx(one, k, y)
                C.cx(y, k)
            C.x(one)
        ra, rb = self.run
        P, M = Prefix(C, H), Prefix(C, L[ra:rb])

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
        self._chain_down(C, L, selectors, masks, M, perch, chain, drop, True, alloc, keep=1)
        c1 = chain[1]
        C.ccx(k0, k1, one)
        C.free([one])

        # On a held lane k1 and k2 change wires, so the compare's target is k's top bit there.
        if p is not None:
            cswap(C, p, [k1], [k2])
        def mark_even(z):
            C.cx(k2, g0)
            C.x(k0)
            C.ccx(g0, k0, z)
            C.x(k0)
            C.cx(k2, g0)

        (z,) = C.alloc(1)
        mark_even(z)
        gt_into(C, a[1:], low[1:], k1, ctrl=z)
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
        cswap(C, sw, a, s)
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
