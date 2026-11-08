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

H's bit 0 is R's, so it is not stored. On a B1 its wire and the next cell's g0, which must be 1,
hold H's two top bits between the ticks, and nothing clean is held across them.
"""
from functools import lru_cache
import heapq

from .arith import Mod, _maj, _selected_carry, add, carries_blocked, crotate, cswap, gt_into_blocked


# Copies of the A1 bit for the tick's last swap: the clean wires the tick's peak leaves there.
FAN = 18
# Blocks of the clear's compare.
COMPARE_BLOCKS = 10
# Blocks of the add in _add_or_shift, or 0 for the single ripple with the rolling XOR.
ADD_BLOCKS = 16
# Copies of J for the shift's rotation.
ROTATE_COPIES = 8


def _compare_blocks(w, m):
    """Block sizes for gt_into_blocked over w bits, low first. Block 0 ripples once and the
    others three times, so block 0 takes two shares; above it sizes rise half a bit per block,
    since the unwind reaches higher blocks first."""
    m = max(1, min(m, w // 2))
    first = 2 * w // (m + 1) if m > 1 else w
    rest = w - first
    cuts = [round(rest * j / (m - 1) + j * (j - m + 1) / 4) for j in range(m)] if m > 1 else [0]
    sizes = [first] + [y - x for x, y in zip(cuts, cuts[1:])]
    assert sum(sizes) == w and min(sizes) > 0
    return sizes


def _fanned_cswap(C, ctrl, a, b, k):
    """cswap(C, ctrl, a, b) with ctrl copied onto k - 1 clean wires, so no control wire is
    shared by more than ceil(len(a) / k) CCX. Pair i reads copy i mod k."""
    k = max(1, min(k, len(a)))
    copies = [ctrl] + C.alloc(k - 1)
    h = 1
    while h < k:
        for j in range(min(h, k - h)):
            C.cx(copies[j], copies[j + h])
        h *= 2
    for i, (x, y) in enumerate(zip(a, b)):
        cswap(C, copies[i % k], [x], [y])
    while h > 1:
        h //= 2
        for j in reversed(range(min(h, k - h))):
            C.cx(copies[j], copies[j + h])
    C.free(copies[1:])


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

    def _chain_up(self, C, L, selectors, masks, M, perch, take, kept=None, moves=None):
        """kept, when given, names the low carries that stay. Each other one is cleared once the
        next kept carry is built, so only kept carries are live above it; _chain_down rebuilds them."""
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

    def _chain_down(self, C, L, selectors, masks, M, perch, chain, drop, write, take=None):
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
        for i in reversed(range(a)):
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

    def _add_or_shift(self, C, J, a, word):
        """On J, word = 2*word. Off J, word += a >> 1, through word's top two bits."""
        if ADD_BLOCKS and self.n > 2:
            return self._add_or_shift_blocked(C, J, a, word)
        n = self.n
        (cin, D, carry, t) = C.alloc(4)
        ops, tgt = a[1:], word[:n - 1]
        carries = [cin] + ops[:-1]
        for i in range(n - 1):
            _maj(C, carries[i], tgt[i], ops[i])
            C.cx(carries[i], D)
        C.x(J)
        C.ccx(J, ops[-1], carry)
        C.ccx(carry, word[n - 1], t)
        C.ccx(t, word[n], word[n + 1])
        C.cx(t, word[n])
        C.cx(carry, word[n - 1])
        C.x(word[n - 1])
        C.ccx(carry, word[n - 1], t)
        C.x(word[n - 1])
        C.ccx(J, ops[-1], carry)
        # The pair below hands word[n - 1] the rolling XOR, which has no pair above to cancel it.
        C.ccx(J, D, word[n - 1])
        C.x(J)
        cswap(C, J, [word[n - 1]], [word[n]])
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
        C.free([cin, D, carry, t])

    def _add_or_shift_blocked(self, C, J, a, word):
        """_add_or_shift's map with the add cut into blocks that add side by side.

        Off J: the carry into every block of word[:n-1] + a[1:] is found at once by
        carries_blocked under not J, each block then adds with its carry in, the carry out of the
        top block increments word's top three bits, and the block carries are cleared by
        carries_blocked again on (~sum + a[1:]), whose boundary carries are [a_low > sum_low], the
        same carries. On J: word[0..n] rotates up one place, word[0] taking word[n], so a word
        whose bit n is 0 doubles; word[n + 1] is left alone, as in the ripple."""
        n = self.n
        x, w = a[1:], word[:n - 1]
        sizes = _compare_blocks(n - 1, ADD_BLOCKS)
        m = len(sizes)
        outs = C.alloc(m)
        C.x(J)
        carries_blocked(C, x, w, outs, J, sizes)
        # Each block's add reads its own copy of not J, so the blocks do not queue on one control.
        copies = [J] + C.alloc(m - 1)
        fan = [(copies[j - (1 << (j.bit_length() - 1))], copies[j]) for j in range(1, m)]
        for src, dst in fan:
            C.cx(src, dst)
        lo = 0
        for j, size in enumerate(sizes):
            add(C, x[lo:lo + size], w[lo:lo + size], ctrl=copies[j], cin=outs[j - 1] if j else None)
            lo += size
        for src, dst in reversed(fan):
            C.cx(src, dst)
        C.free(copies[1:])
        c = outs[-1]
        (t,) = C.alloc(1)
        C.ccx(c, word[n - 1], t)
        C.ccx(t, word[n], word[n + 1])
        C.cx(t, word[n])
        C.cx(c, word[n - 1])
        C.x(word[n - 1])
        C.ccx(c, word[n - 1], t)
        C.x(word[n - 1])
        C.free([t])
        for y in w:
            C.x(y)
        carries_blocked(C, x, w, outs, J, sizes)
        for y in w:
            C.x(y)
        C.x(J)
        C.free(outs)
        copies = C.alloc(ROTATE_COPIES - 1)
        crotate(C, [J] + copies, word[n::-1])
        C.free(copies)

    def tick(self, C, prev, g, nxt, a, s, pebble=False):
        """One tick of a swept tape. prev is the cell before, or None on the first tick.
        nxt is the next cell's g0 wire, or None on the last tick, whose B1 reduces normally."""
        n, f, w = self.n, self.f, self.window
        g0, b = g
        (spill, high) = C.alloc(2)
        word = s + [spill, high]
        p = self._after_b1(C, prev)
        if p is not None:
            cswap(C, p, [s[0]], [spill])
            cswap(C, p, [g0], [high])
            C.cx(p, g0)
        sw = self._begin(C, g0, b)
        if p is not None:
            for x in word:
                C.cx(p, x)
        self._add_or_shift(C, b, a, word)
        if p is not None:
            for x in word:
                C.cx(p, x)
        self._drop_after_b1(C, prev, p)
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
            cswap(C, one, [k0], [l0])
            C.cx(one, nxt)
            cswap(C, one, [k1], [nxt])
            C.x(one)
        ra, rb = self.run
        P, M = Prefix(C, H), Prefix(C, L[ra:rb])
        # Only M's skip reads M before its write, so the top below its last is cleared until then.
        M.release()

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
        chain = self._chain_up(C, L, selectors, masks, M, {}, take)
        (e,) = C.alloc(1)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        self._chain_down(C, L, selectors, masks, M, {}, chain, drop, False)
        if nxt is None:
            C.x(one)
        else:
            C.cx(b, one)
            C.cx(g0, one)
        C.cx(sw, b)
        self._end(C, g0, b, sw)
        C.ccx(three, e, k2)
        C.ccx(k0, k1, three)
        C.free([three])
        C.ccx(k0, e, k1)
        C.cx(e, k0)
        if spare:
            u = take()
            C.cx(e, u)
            C.ccx(P.top, u, e)
            C.free([e])
            e = u
        C.ccx(k0, k1, one)
        selectors, masks = [one, k0, k1, k2], [self.corner, f, 2 * f, 4 * f]
        perch = self._perches(selectors, masks)
        kept = self._kept_carries(perch)
        schedule = self._pebble_schedule(perch) if pebble else None
        moves = None
        if schedule is not None:
            moves, low_kept = schedule
            kept = set(low_kept) | {self.run[0]}
        chain = self._chain_up(C, L, selectors, masks, M, perch, alloc, kept, moves)
        if H:
            C.ccx(P.top, chain[w], e)
        else:
            C.cx(chain[w], e)
        drop(e)
        P.write(chain[w])
        self._chain_down(C, L, selectors, masks, M, perch, chain, drop, True, alloc)
        C.ccx(k0, k1, one)
        C.free([one])

        # On a held lane k1 and k2 change wires, so the compare's target is k's top bit there.
        p = self._after_b1(C, prev)
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
        # z = 1 forces k0 = 0, and wire k2 = 0 on the domain by S7, so both stand in for clean
        # wires. Three idle wires of any value serve as carry ins: the next g0, R's bit 0 and b
        # have been free longest, prev and p next. A pebbled tick has one clean wire fewer here.
        dirty = [x for x in (nxt, a[0], b, *(prev or []), p) if x is not None][:3]
        sizes = _compare_blocks(n - 1, COMPARE_BLOCKS - pebble)
        gt_into_blocked(C, a[1:], low[1:], k1, z, sizes, [k0, k2], dirty)
        mark_even(z)
        C.free([z])
        if p is not None:
            # S = R + k*f mod 4 names k's bit 1, which now sits on k2.
            def bit_one(m):
                C.ccx(k0, a[0], m)
                C.cx(low[1], m)
                C.cx(a[1], m)
                if f >> 1 & 1:
                    C.cx(k0, m)

            (m,) = C.alloc(1)
            bit_one(m)
            C.ccx(p, m, k2)
            bit_one(m)
            C.free([m])
        self._drop_after_b1(C, prev, p)
        sw = self._begin(C, g0, b)
        C.ccx(b, low[1], k1)
        if f >> 1 & 1:
            C.ccx(b, k0, k1)
        C.cx(l0, k0)
        if nxt is not None:
            C.cx(g0, b)
            C.x(b)
            C.ccx(b, k0, l0)
            C.x(b)
            C.cx(g0, b)
        C.cx(sw, b)
        C.ccx(g0, a[0], l0)
        C.free([l0, spill, high])
        s[:] = s[-1:] + s[:-1]
        _fanned_cswap(C, sw, a, s, FAN)
        self._end(C, g0, b, sw)
