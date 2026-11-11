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

from . import wide
from .arith import Mod, _Gates, _emit, _maj, _selected_carry, adjacent_xor, cswap, gt_into


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

    alloc_restore = None

    def restore(self):
        if self.released is None:
            return
        a, p, x, u, v = self.blocks[self.released]
        (q,) = [self.alloc_restore()] if self.alloc_restore else self.C.alloc(1)
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
            # Each bit's increment runs before the uncompute that ends its chain: x ^= p*(not below)
            # leaves below*x as it was for the uncompute, and the x ^= p that completes x ^= p*below
            # waits for the block's last bit, since no later bit reads x[i] or writes p.
            for i in reversed(range(1, len(x))):
                below = x[0] if i == 1 else u[i - 2]
                C.x(below)
                C.ccx(p, below, x[i])
                C.x(below)
                C.ccx(below, x[i], u[i - 1])
            for i in reversed(range(1, len(x))):
                C.cx(p, x[i])
            C.cx(p, x[0])
            for w in reversed(u):
                C.x(w)
            C.ccx(t, a, p)
            C.free([p])
        C.x(bits[0])
        C.ccx(t, bits[0], bits[1])
        C.x(bits[0])
        C.ccx(bits[0], bits[1], self.first)
        C.cx(t, bits[1])
        C.cx(t, bits[0])
        C.free([self.first])


# Cuts for the tick's 255-bit compare by block count, tuned on the emitted tick definitions 260 and 261
# (pass off, weighted 546:216) with the wires this function hands out. With 11 blocks: a 38-bit
# single-ripple block 0, two short blocks, then blocks of about 24 bits.
CUTS_255 = {
    7: (0, 55, 87, 120, 153, 186, 220, 255),
    8: (0, 58, 76, 105, 134, 164, 194, 224, 255),
    9: (0, 44, 58, 78, 107, 136, 166, 196, 226, 255),
    10: (0, 41, 54, 74, 99, 125, 151, 177, 203, 231, 255),
    11: (0, 38, 48, 61, 84, 108, 132, 155, 180, 204, 231, 255),
    12: (0, 44, 54, 65, 77, 90, 113, 134, 157, 182, 206, 232, 255),
}
MAX_BLOCKS = 11
# The tick's add as a blocked-lookahead chain (add_or_shift_la) instead of the ripple.
LOOKAHEAD_ADD = True
# The add's temps and J copies: the wires P's tenth block adds at the tick's peak are idle here.
LA_NTEMPS, LA_NCOPIES = 8, 6
# (size, count): that many wider blocks among F's / B's blocks of 2, or None for blocks of 2 only. An
# r-bit block costs 4r - 3 CCX, so each chain step saved costs 3 CCX for any r (DESIGN.md, "Blocks of
# 3 in the add's forward sweep"). Split between F and B, the sweep that runs at 7 ticks a step and
# starves for temps differs by direction ("Blocks of 3 split between the sweeps").
LA_F_WIDE = (3, 34)
LA_B_WIDE = (3, 35)
# The add's carry and t serve as temps in its sweeps; temps go out free-first under the tick model.
LA_TOP_TEMPS = True
LA_IDLE_TEMPS = True
# The tick's P as PrefixR, whose write runs its blocks in parallel. At 221 bits: ten blocks, the top
# one small, since the write starts with the top block's ladder (DESIGN.md, "Ten P blocks" and "P's
# blocks after the faster add").
PARALLEL_PREFIX = True
P_SIZES_221 = (3, 7, 14, 28, 28, 28, 28, 28, 28, 27)
HIGH_CARRIERS = True
# Put the final swap's sw and control copies on scratch the next tick reaches only after its add.
LATE_FINAL_WIRES = True
LATE_MIN_COPIES = 6
# Deal the final swap's bits to sw's copies in turn instead of in slices: the compare restores a
# block's bits together, and in a slice they would wait for one copy in turn.
FINAL_INTERLEAVE = True
# The final swap's first this many bits run on sw before its fan.
FINAL_EARLY = 1
# Copies of p for un-copying it from the word after the add (1 = none).
P_FAN = 32
# Run the tick's add as 2elliot depth/v2's CDKM ripple instead of either add above: 579,096
# fewer CCX on the whole circuit for 1,535,700 more ticks (measured after the post-pass).
CDKM_ADD = False
# Hand P's and M's tops, and the H wires the first chain borrows, out in the order this tick's ops
# so far free them, so neither waits on a wire the add holds to its end.
IDLE_WIRES = True
# M's restored top too, instead of the carrier the second chain has just dropped.
RESTORE_IDLE = True
# The tick's compare as wide.gt_tree on clean scratch, in blocks of this many bits (0 = compare_blocks).
WIDE_CMP = (2, 3, 4, 5)
# The tick's add in blocks of this many bits, each with its own temps (0 = LA_F_WIDE / LA_B_WIDE).
WIDE_ADD = 9
# The add's blocks ramp from 1 up to WIDE_ADD at both ends, each size this many times.
WIDE_ADD_RAMP = 1
WIDE_ADD_COPIES = 24
# The add's backward sweep in its own blocks, (top, ramp) as above (None = the forward sweep's).
# Its chain costs 11 ticks a block against the forward chain's 7, so its blocks are larger.
WIDE_ADD_BACK = (10, (2, 1))
# Explicit block sizes, bottom up, for the forward and backward sweeps at the width they sum to.
WIDE_ADD_FSIZES = (1, 1, 2, 2, 3, 4, 5, 6, 6, 8, 9, 9, 11, 10, 12, 9, 6, 14, 10, 8, 9, 10, 15, 9, 13, 9, 11, 11, 5, 8, 6, 6, 2, 3, 2)
WIDE_ADD_BSIZES = (1, 1, 1, 1, 3, 4, 5, 3, 6, 6, 5, 8, 7, 7, 9, 5, 8, 10, 10, 10, 9, 9, 7, 10, 11, 11, 10, 11, 11, 10, 9, 8, 6, 7, 6, 4, 3, 2, 1)
# Run the add's first CX wave through U as a log-depth network, and fan p into the word through
# copies. The forward sweep's blocks then start together instead of one bit a tick.
WIDE_WAVE = True
# The first wave's network in segments of this many bits (0 = one network).
WAVE_SEG = 64
# Hand out the add's copies of J by when each comes free, not in turn.
COPY_IDLE = True
# The second wave too. Left serial, it follows the backward sweep down the word, and the rotation
# starts on the middle bits meanwhile; as a network it waits for the sweep's last bit.
WIDE_WAVE_BACK = False
# The tick's P and M as wide.PrefixTree in blocks of this many bits (0 = PrefixR and Prefix). With
# WIDE_P the first chain's carriers are clean wires, not borrowed H wires.
WIDE_P = 4
WIDE_M = 3
# The tick's two low chains as wide.LowAdd in blocks of this many bits (0 = _chain_up and _chain_down).
WIDE_LOW = 3
# P's and M's writes in PrefixTree._write_pebbled's order (None: write's own). P clears its whole
# tree going down and may take fresh wires; M only rebuilds its tree node by node coming up.
P_PEBBLE = 99
M_PEBBLE = -1
P_FRESH = 39
M_FRESH = 0
# Fewer fresh wires for P on a pebble tick: the pass holds a second decoded block across a block
# edge, one wire more than elsewhere, so this keeps those ticks at the others' peak.
P_FRESH_EDGE = 1
# The compare's tree as two balanced trees, over its bits below this and the rest, joined at the
# root (0 = one balanced tree).
CMP_SPLIT = 14
# Give every scratch wire from the prefixes to the compare by when it comes free under the tick model,
# earliest first, instead of off the stack: a wire one tree has just freed is busy until that tree
# ends, and handing it to a tree that could run beside it makes the second wait for the first.
WIDE_IDLE = True
# Hand out the tick's wires by when they come free from its first gate, not only from its prefixes.
IDLE_FROM_START = True
# Take the wire for one before p's mark is dropped (see tick).
ONE_BEFORE_DROP = True
# With WIDE_ASIDE, emit the trial's undo before the second chain.
ASIDE_UNDO_FIRST = True
# Turn three into the corner selector k0' k1' instead of clearing it (see tick).
KUPD_FAST = True
# With WIDE_ASIDE, the second chain on wires the tick has not touched yet.
SECOND_FRESH = True
# Run the trial chain on copies, beside the second chain.
WIDE_ASIDE = True


# alloc_after counts a wire as free in time when it is free within this many ticks of the gate's
# other wires, and may take up to IDLE_FRESH wires the tick has not touched when none is.
IDLE_SLACK = -8
IDLE_FRESH = 0
# Scratch wires a tick may touch before _untouched hands out touched wires instead (0 = no cap).
TOUCH_CAP = 318


class _Idle:
    """When each wire this tick has touched comes free, under the board's tick model with the
    tick alone (every wire free at its start)."""

    def __init__(self, C, start, io=()):
        self.C, self.pos, self.r, self.fresh = C, start, {}, 0
        # The tick's own registers, and every wire handed out: what counts toward TOUCH_CAP.
        self.io, self.handed = set(io), set()

    def ready(self):
        r, ops = self.r, self.C.prims
        for pos in range(self.pos, len(ops), 4):
            k, a, b, c = ops[pos:pos + 4]
            if k == 0:
                r[a] = r.get(a, 0) + 1
            elif k == 2:
                s2 = max(r.get(b, 0), r.get(c, 0) + 1) + 1
                s4 = max(r.get(a, 0), s2 + 1) + 1
                r[c], r[a], r[b] = s4 + 6, s4 + 7, s4 + 7
            else:
                t = max(r.get(a, 0), r.get(b, 0)) + (1 if k == 1 else 3)
                r[a] = r[b] = t
        self.pos = len(ops)
        return r

    def first(self, wires):
        """Index of the wire in wires that comes free first, the later one on a tie."""
        r = self.ready()
        return min(range(len(wires)), key=lambda i: (r.get(wires[i], 0), -i))

    def alloc_after(self, *inputs):
        """One clean wire for a gate whose other wires are `inputs`: the touched free wire that
        comes free last among those free by the time the inputs are, so that the gate does not wait
        for its target and earlier wires stay for earlier gates. With none, as alloc."""
        C, r = self.C, self.ready()
        t = max((r.get(w, 0) for w in inputs), default=0) + IDLE_SLACK
        fit = [w for w in C.free_stack if w in r and r[w] <= t]
        if not fit:
            if self.fresh < IDLE_FRESH:
                self.fresh += 1
                return self._untouched()
            return self.alloc(1)[0]
        w = max(fit, key=lambda x: r[x])
        C.free_stack.remove(w)
        return w

    def _untouched(self):
        """A clean wire this tick has not touched: off the stack when it holds one, else a new wire.
        With TOUCH_CAP scratch wires touched already, the touched free wire that comes free first."""
        C, r = self.C, self.ready()
        if TOUCH_CAP and len((set(r) | self.handed) - self.io) >= TOUCH_CAP:
            touched = [w for w in C.free_stack if w in r]
            if touched:
                w = touched[self.first(touched)]
                C.free_stack.remove(w)
                return w
        w = self._new()
        self.handed.add(w)
        return w

    def _new(self):
        C, r = self.C, self.r
        for w in reversed(C.free_stack):
            if w not in r:
                C.free_stack.remove(w)
                return w
        C.num_wires += 1
        return C.num_wires - 1

    def alloc(self, k=1):
        """Clean wires this tick has touched, earliest free first; a fresh wire only when none is free,
        as the stack would, so the tick touches no more wires than before."""
        C, out = self.C, []
        for _ in range(k):
            r = self.ready()
            touched = [w for w in C.free_stack if w in r]
            if not touched:
                out += type(C).alloc(C, 1)
                self.handed.add(out[-1])
                continue
            w = touched[self.first(touched)]
            C.free_stack.remove(w)
            out.append(w)
        return out


def weighted_bounds(width, blocks):
    """The generator's current cut rule, kept as the fallback for other widths and block counts."""
    weights = [400] + [200 + 9 * (2 * j - blocks) for j in range(1, blocks)]
    total, acc, bounds = sum(weights), 0, [0]
    for weight in weights[:-1]:
        acc += weight
        bounds.append((2 * width * acc + total) // (2 * total))
    bounds.append(width)
    for j in range(1, blocks):
        bounds[j] = max(bounds[j], bounds[j - 1] + 1)
    for j in reversed(range(1, blocks)):
        bounds[j] = min(bounds[j], bounds[j + 1] - 1)
    return bounds


def assign_wires(blocks, cin, clean, dirty, last):
    """Carry-ins and gs wires for blocks 1..blocks-1, handed out top block first.

    Every block j >= 1 needs one clean gs wire and one carry-in of any value. Clean wires beyond the
    gs wires serve as carry-ins first, then dirty wires make up the rest, so a dirty wire is borrowed
    only when the clean ones run out. Within that, each block takes whichever of the next clean and
    the next dirty wire came free first."""
    m = blocks - 1
    order = lambda w: last.get(w, -1)
    clean = sorted(clean, key=order)
    n_dirty = max(0, m - (len(clean) - m))
    dirty = sorted(dirty, key=order)[:n_dirty]
    assert len(clean) >= m and len(dirty) == n_dirty
    pool = clean[:2 * m - n_dirty]
    quota = m - n_dirty
    carry_ins, gs = [None] * blocks, [None] * blocks
    for j in reversed(range(1, blocks)):
        options = []
        if dirty:
            options.append((order(dirty[0]), "dirty"))
        if quota and pool:
            options.append((order(pool[0]), "clean"))
        _, kind = min(options)
        if kind == "dirty":
            carry_ins[j] = dirty.pop(0)
        else:
            carry_ins[j] = pool.pop(0)
            quota -= 1
        gs[j] = pool.pop(0)
    return carry_ins, gs


def compare_blocks(C, A, B, z, t, cin, clean, dirty, last, max_blocks=MAX_BLOCKS):
    width = len(A)
    assert len(B) == width and width > 0
    m = min(max_blocks - 1, width - 1, len(clean), (len(clean) + len(dirty)) // 2)
    blocks = m + 1
    if width == 255 and blocks in CUTS_255:
        bounds = list(CUTS_255[blocks])
    else:
        bounds = weighted_bounds(width, blocks)
    carry_ins, gs = assign_wires(blocks, cin, clean, dirty, last)
    carry_ins[0] = cin
    known_zero = {cin, *clean}

    # Block j's chain: its CX gates, then its CCX ladder, in the generator's interleaving.
    chains = []
    for j in range(blocks):
        lo, hi = bounds[j:j + 2]
        carry = carry_ins[j]
        chain = [(C.cx, (A[lo], B[lo])), (C.cx, (A[lo], carry))]
        for i in range(lo, hi):
            if i + 1 < hi:
                chain += [(C.cx, (A[i + 1], B[i + 1])), (C.cx, (A[i + 1], A[i]))]
            chain.append((C.ccx, (carry if i == lo else A[i - 1], B[i], A[i])))
        chains.append(chain)

    ops = []

    def chain(j, reverse=False, ladder=False):
        seq = [g for g in chains[j] if not ladder or g[0] == C.ccx]
        ops.extend(reversed(seq) if reverse else seq)

    tops = [A[bounds[j + 1] - 1] for j in range(blocks)]
    ops.extend((C.x, (w,)) for w in B)
    for j in range(blocks):
        chain(j)
    # Block j's top now holds v1 = g ^ e*d (d its carry-in's value); copy it, then flip d and redo the
    # ladder so the top holds v2 = g ^ e*(1 ^ d). Only the ladder reads the carry-in, so only the
    # ladder is undone and redone.
    for j in range(1, blocks):
        ops.append((C.cx, (tops[j], gs[j])))
    for j in range(1, blocks):
        chain(j, reverse=True, ladder=True)
    for j in range(1, blocks):
        ops.append((C.x, (carry_ins[j],)))
        chain(j, ladder=True)
    # Carry out of block j for carry in c is MAJ(v1, v2, c); with d = 0 it is v1 ^ (v1 ^ v2)*c,
    # which leaves c's wire alone.
    carry = tops[0]
    for j in range(1, blocks):
        g, top = gs[j], tops[j]
        ops.append((C.cx, (g, top)))
        if carry_ins[j] not in known_zero:
            ops.append((C.cx, (g, carry)))
        ops.append((C.ccx, (carry, top, g)))
        carry = g
    for gate, args in ops:
        gate(*args)
    C.ccx(z, carry, t)
    for gate, args in reversed(ops):
        gate(*args)


# Blocked-lookahead add-or-shift (design swarm, work/swarm/add/design2.py), same map as the ripple.
def _la_fan(C, src, dsts):
    sources, tree, filled = [src], [], 0
    while filled < len(dsts):
        layer = []
        for s in list(sources):
            if filled == len(dsts):
                break
            t = dsts[filled]
            filled += 1
            C.cx(s, t)
            tree.append((s, t))
            layer.append(t)
        sources.extend(layer)
    return tree


def _la_spread(m, k, ends=1):
    """Chain block sizes, bottom up: k blocks of 2 spread evenly among ripple steps, keeping
    `ends` ripple steps at each end."""
    k = max(0, min(k, (m - 2 * ends) // 2))
    if k == 0:
        return [1] * m
    inner = m - 2 * k - 2 * ends
    gaps = [inner * (i + 1) // (k + 1) - inner * i // (k + 1) for i in range(k + 1)]
    out = [1] * ends
    for i in range(k):
        out += [1] * gaps[i] + [2]
    out += [1] * gaps[k] + [1] * ends
    assert sum(out) == m
    return out


def _la_wide(m, big, k, ends=1):
    """Chain block sizes, bottom up: k blocks of `big` spread evenly among blocks of 2, keeping `ends`
    ripple steps at each end; a bit left over is one more ripple step in the middle."""
    inner = m - 2 * ends
    k = max(0, min(k, inner // big))
    n2, odd = divmod(inner - big * k, 2)
    mid = [big if (i + 1) * k // (k + n2) > i * k // (k + n2) else 2 for i in range(k + n2)]
    return [1] * ends + mid[:len(mid) // 2] + [1] * odd + mid[len(mid) // 2:] + [1] * ends


def _ramp_sizes(m, top, rep):
    """Chain block sizes for m bits, bottom up: 1, then each size from 2 to top - 1 `rep` times,
    blocks of `top` through the middle, and the same ramp down to 1 at the upper end. A block's
    products take time in proportion to its size before the chain can cross it, and its fixes take
    as long after, so a block near either end holds the sweep up unless it is short. A width too
    small for both ramps keeps as much of each as fits; what is left over is one middle block."""
    lo_rep, hi_rep = rep if isinstance(rep, tuple) else (rep, rep)
    ramps = [[1] + [k for k in range(2, top) for _ in range(r)] for r in (lo_rep, hi_rep)]
    ends, left = ([], []), m
    for i in range(max(map(len, ramps))):
        for end, ramp in zip(ends, ramps):
            if i < len(ramp) and len(end) == i and left >= ramp[i] + top:
                end.append(ramp[i])
                left -= ramp[i]
    low, high = ends
    mid = [top] * (left // top)
    left -= sum(mid)
    out = low + mid[:len(mid) // 2] + ([left] if left else []) + mid[len(mid) // 2:] + high[::-1]
    assert sum(out) == m and all(r >= 1 for r in out)
    return out


def _la_blocks(sizes):
    out, s = [], 0
    for r in sizes:
        out.append((s, r))
        s += r
    return out


class _LaTemps:
    """Clean temps in turn; a returned temp goes to the back. With idle, the temp that comes free
    first under the tick model, the earliest returned on a tie."""

    def __init__(self, wires, idle=None):
        self.free, self.idle = list(wires), idle

    def take(self):
        assert self.free, "out of temps"
        if self.idle:
            r = self.idle.ready()
            i = min(range(len(self.free)), key=lambda i: (r.get(self.free[i], 0), i))
            return self.free.pop(i)
        return self.free.pop(0)

    def give(self, w):
        self.free.append(w)


def add_or_shift_la(C, J, a, word, n, k2f=None, k2b=None, ntemps=5, ncopies=6, fsizes=None,
                    bsizes=None):
    """k2f / k2b: number of 2-blocks in F / B (None = as many as fit); fsizes / bsizes, when given,
    F's / B's block sizes bottom up instead. Wires: cin, ntemps temps, carry and t for the top,
    ncopies copies of J: 1 + ntemps + 2 + ncopies (14 by default)."""
    m = n - 1
    if k2f is None:
        k2f = m
    if k2b is None:
        k2b = m
    (cin,) = C.alloc(1)
    tw = C.alloc(ntemps)
    carry, t = C.alloc(2)
    K = C.alloc(ncopies)
    U = [cin] + a[1:]
    V = word[:m]

    idle = _Idle(C, len(C.prims)) if LA_IDLE_TEMPS else None
    # J holds ~J from here to the end; the copies hold ~J during F and J after it.
    C.x(J)
    tree = _la_fan(C, J, K)
    ki = [0]

    def copy(*others):
        """A copy of J for a gate on `others`: the copy that comes free last among those free by
        the time the others are, else the one free first; without the tick model, in turn. In turn,
        a bit whose gate is ready waits for the copy's last bit, which may sit deep in a large
        block."""
        if idle is not None and COPY_IDLE:
            r = idle.ready()
            t = max((r.get(w, 0) for w in others), default=0)
            fit = [g for g in K if r.get(g, 0) <= t]
            return max(fit, key=lambda g: r.get(g, 0)) if fit else min(K, key=lambda g: r.get(g, 0))
        return turn()

    def turn():
        w = K[ki[0] % len(K)]
        ki[0] += 1
        return w

    if WIDE_WAVE:
        # The same map as the loop below, whose CXs chain through U one tick a bit: every V[j] takes
        # the old U[j + 1], then U takes its adjacent XORs as a log-depth network.
        for j in range(m):
            C.cx(U[j + 1], V[j])
        wave = _Gates()
        # In segments bottom up, each a network on its own bits that reads the next one's lowest
        # bit before that segment changes it, so the low bits need not wait for the high ones.
        seg = WAVE_SEG or len(U)
        for lo in range(0, len(U) - 1, seg):
            adjacent_xor(wave, U[lo:min(lo + seg, len(U) - 1) + 1])
        _emit(C, wave.ops)
    else:
        for j in range(m):
            C.cx(U[j + 1], V[j])
            C.cx(U[j + 1], U[j])

    # carry and t are clean until the top, and every temp is returned clean, so they can serve as temps.
    temps = _LaTemps([carry, t] + tw if LA_TOP_TEMPS else tw, idle)
    for s, r in _la_blocks(fsizes or _la_spread(m, k2f)):
        if r == 1:
            C.ccx(U[s], V[s], U[s + 1])
        else:
            for i in range(1, r):                       # local: U[s+i+1] = its carry when f_s = 0
                C.ccx(U[s + i], V[s + i], U[s + i + 1])
            T = [V[s]]                                  # T[i] = p_s ... p_(s+i)
            for i in range(1, r):
                T.append(temps.take())
                C.ccx(T[i - 1], V[s + i], T[i])
            C.ccx(U[s], T[r - 1], U[s + r])             # the chain step
            # Every fix before any product is cleared: a fix and the next clear both hold T[i - 1],
            # and taking turns would make the fixes wait for the clears.
            for i in reversed(range(1, r)):
                C.ccx(U[s], T[i - 1], U[s + i])         # fix the inner bits, top down
            for i in reversed(range(1, r)):
                C.ccx(T[i - 1], V[s + i], T[i])
                temps.give(T[i])
        for j in range(s, s + r):
            g = copy(U[j], V[j])
            C.cx(g, V[j])
            C.ccx(U[j], g, V[j])                        # W_j: V[j] ^= ~J & ~f_j
    for g in K:
        C.x(g)

    # word[m..m+2] += ~J c_m, with c_m on U[m]: the generator's gates, J holding ~J.
    cm = U[m]
    C.ccx(cm, J, carry)
    C.ccx(carry, word[n - 1], t)
    C.ccx(t, word[n], word[n + 1])
    C.cx(t, word[n])
    C.cx(carry, word[n - 1])
    C.x(word[n - 1])
    C.ccx(carry, word[n - 1], t)
    C.x(word[n - 1])
    C.ccx(cm, J, carry)

    # B runs F's blocks backwards: f V = f p on every lane, so the same steps undo the carries.
    bb = _la_blocks(bsizes or _la_spread(m, k2b))
    # Each block's inner unfixes read U[s] = f_s, which the next block down overwrites. The first
    # needs no temp, so all of those go first; then the chain runs down one CCX a block.
    for s, r in bb:
        if r >= 2:
            C.ccx(U[s], V[s], U[s + 1])
    for s, r in reversed(bb):
        if r == 1:
            C.ccx(U[s], V[s], U[s + 1])
        else:
            # The products as one chain, then the unfixes, then the chain step: an unfix between
            # two products would hold the lower one. The carries unwind before the products clear,
            # since the next block down waits for the carries only.
            T = [V[s]]
            for i in range(1, r):
                T.append(temps.take())
                C.ccx(T[i - 1], V[s + i], T[i])
            for i in range(1, r - 1):
                C.ccx(U[s], T[i], U[s + i + 1])
            C.ccx(U[s], T[r - 1], U[s + r])
            for i in reversed(range(1, r)):
                C.ccx(U[s + i], V[s + i], U[s + i + 1])
            for i in reversed(range(1, r)):
                C.ccx(T[i - 1], V[s + i], T[i])
                temps.give(T[i])

    if WIDE_WAVE and WIDE_WAVE_BACK:
        _emit(C, reversed(wave.ops))
        for j in reversed(range(m)):
            C.cx(U[j + 1], V[j])
            C.cx(copy(V[j]), V[j])                      # V ^= ~J through J's copies, which hold J
            C.x(V[j])
    else:
        for j in reversed(range(m)):
            C.cx(U[j + 1], U[j])
            C.cx(U[j + 1], V[j])
            C.cx(J, V[j])                               # V ^= ~J

    def cswap(c, x, y):
        C.cx(y, x)
        C.ccx(c, x, y)
        C.cx(y, x)

    for i in reversed(range(n)):
        if i < (n + 1) // 2:
            cswap(turn(), word[i], word[n - i])
        if i < n // 2:
            cswap(turn(), word[i + 1], word[n - i])

    C.x(J)
    for s_, t_ in reversed(tree):
        C.cx(s_, t_)
    C.free(K)
    C.free(tw)
    C.free([cin])
    C.free([carry, t])      # top of the stack: idle first, so the tick's next alloc gets them


def add_or_shift_cdkm(C, J, a, word, n, add_copies=4, shift_copies=4):
    """2elliot depth/v2 (9d34aef) TwoBitMod._add_or_shift: on J, word = 2*word; off J,
    word += a >> 1 through word's top two bits. A CDKM ripple whose sum writes are gated on not J
    and moved off its carry chain, then the shift as two reversals whose swaps follow the down pass."""
    m = n - 1
    ops, tgt = a[1:], word[:m]
    (cin, carry, t) = C.alloc(3)
    carries = [cin] + ops[:-1]
    js = C.alloc(shift_copies)
    for w in js:
        C.cx(J, w)
    C.x(J)
    adds = [J] + C.alloc(add_copies - 1)
    fan = [(adds[j - (1 << (j.bit_length() - 1))], adds[j]) for j in range(1, add_copies)]
    for src, dst in fan:
        C.cx(src, dst)
    for i in range(m):
        C.cx(ops[i], tgt[i])
    C.cx(ops[0], cin)
    for i in range(m):
        if i + 1 < m:
            C.cx(ops[i + 1], ops[i])
        C.ccx(carries[i], tgt[i], ops[i])
        C.x(carries[i])
        C.ccx(adds[i % add_copies], carries[i], tgt[i])
        C.x(carries[i])
    top = adds[m % add_copies]
    C.ccx(top, ops[-1], carry)
    C.ccx(carry, word[n - 1], t)
    C.ccx(t, word[n], word[n + 1])
    C.cx(t, word[n])
    C.cx(carry, word[n - 1])
    C.x(word[n - 1])
    C.ccx(carry, word[n - 1], t)
    C.x(word[n - 1])
    C.ccx(top, ops[-1], carry)
    for src, dst in reversed(fan):
        C.cx(src, dst)
    swaps = 0
    for i in reversed(range(m)):
        C.ccx(carries[i], tgt[i], ops[i])
        if i + 1 < m:
            C.cx(ops[i + 1], ops[i])
        C.cx(J, tgt[i])
        C.cx(ops[i], tgt[i])
        if i < n - i:
            cswap(C, js[swaps % shift_copies], [word[i]], [word[n - i]])
            swaps += 1
        if i + 1 < n - i:
            cswap(C, js[swaps % shift_copies], [word[i + 1]], [word[n - i]])
            swaps += 1
    C.cx(ops[0], cin)
    for w in js[1:]:
        C.cx(js[0], w)
    C.x(J)
    C.cx(J, js[0])
    C.free([carry, t] + adds[1:] + js + [cin])


# Block-parallel prefix write for the tick's P (design swarm, work/swarm/prefix/design_r.py).
def tree_gates(x, u, seg=2):
    """Gates putting AND(x) on carriers u (len(x) - 1 of them, used in order): sub-chains of at most
    seg leaves, then a balanced tree over their tops. Returns (gates, root)."""
    m = len(x)
    if m == 1:
        return [], x[0]
    gates, tops, pos = [], [], 0
    for a in range(0, m, seg):
        part = x[a:a + seg]
        v = part[0]
        for i in range(1, len(part)):
            gates.append((v, part[i], u[pos]))
            v = u[pos]
            pos += 1
        tops.append(v)
    level = tops
    while len(level) > 1:
        nxt = []
        for k in range(0, len(level) - 1, 2):
            gates.append((level[k], level[k + 1], u[pos]))
            nxt.append(u[pos])
            pos += 1
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    assert pos == m - 1
    return gates, level[0]


class PrefixR:
    """sizes: block sizes after the first two bits (each at most 1 + the carriers below it);
    cap: today's rule (grow as fast as carriers allow, at most cap) when sizes is None;
    seg, host, reorder: W2 tree shape, hosting and nested carrier order (defaults are the design)."""

    def __init__(self, C, bits, sizes=None, cap=None, seg=2, host=True, reorder=True):
        self.C, self.bits, self.seg, self.host, self.reorder = C, bits, seg, host, reorder
        self.starts = []
        self.blocks = []
        h = len(bits)
        self.top = bits[0] if h == 1 else None
        if h < 2:
            return
        (self.first,) = C.alloc(1)
        C.ccx(bits[0], bits[1], self.first)
        tops, start = [self.first], 2
        sizes = list(sizes) if sizes is not None else None
        while start < h:
            carriers = bits[:start] + tops[:-1]
            if sizes:
                m = sizes.pop(0)
                assert 1 <= m <= len(carriers) + 1 and start + m <= h, (m, start)
            else:
                m = min(len(carriers) + 1, h - start, cap or h)
            x, u = bits[start:start + m], carriers[:m - 1]
            if HIGH_CARRIERS and start >= m - 1:
                # Any wire below start reads 1 when tops[-1] is, so the chain can sit on the bits just
                # below the block, which frees the lowest bits once the first blocks are built.
                u = bits[start - (m - 1):start]
            for w in u:
                C.x(w)
            v = x[0]
            for i in range(1, m):
                C.ccx(v, x[i], u[i - 1])
                v = u[i - 1]
            (p,) = C.alloc(1)
            C.ccx(tops[-1], v, p)
            self.blocks.append((tops[-1], p, x, u, v))
            self.starts.append(start)
            tops.append(p)
            start += m
        assert sizes is None or not sizes
        self.top = tops[-1]

    def write(self, t):
        C, bits, blocks = self.C, self.bits, self.blocks
        if not bits:
            return
        if len(bits) == 1:
            C.cx(t, bits[0])
            return
        # W1
        for a, p, x, u, v in reversed(blocks):
            C.cx(t, v)
            C.ccx(a, v, p)
            C.cx(t, v)
            m = len(x)
            for i in reversed(range(1, m)):
                below = x[0] if i == 1 else u[i - 2]
                C.x(below)
                C.ccx(p, below, x[i])
                C.x(below)
                C.ccx(below, x[i], u[i - 1])
                C.x(u[i - 1])
            for i in reversed(range(m)):
                C.cx(p, x[i])
        C.x(bits[0])
        C.ccx(t, bits[0], bits[1])
        C.x(bits[0])
        C.ccx(bits[0], bits[1], self.first)
        C.cx(t, bits[1])
        C.cx(t, bits[0])
        if not blocks:
            C.free([self.first])
            return
        # W2: F (t-free, recorded) with t's clears interleaved.
        F = []

        def f(op, *w):
            F.append((op, w))
            (C.x if op == "x" else C.ccx)(*w)

        leaves = bits[:len(bits) - len(blocks[-1][2])]      # the top block's bits are no leaf
        for w in leaves:
            f("x", w)
        f("ccx", bits[0], bits[1], self.first)
        # Plan every block's carriers first, so a deferred vv restore can be checked against
        # every later use of that wire.
        plan, order = [], [bits[0], bits[1]]
        gate = self.first
        for j, (a, p, x, u, v) in enumerate(blocks[:-1]):
            m = len(x)
            hosted = self.host and j >= 1 and len(blocks[j - 1][2]) >= m - 1
            if hosted:
                carriers = blocks[j - 1][2][:m - 1]
            else:
                carriers = (order if self.reorder else u)[:m - 1]
                assert set(carriers) <= set(u) | set(bits[:self.starts[j]]), "carrier above block"
            order = list(x) + order + [gate]
            gates, vv = tree_gates(x, carriers, self.seg)
            plan.append((gate, p, carriers, hosted, gates, vv))
            gate = p
        deferred = []
        for j, (gate, p, carriers, hosted, gates, vv) in enumerate(plan):
            # Nested carriers read 1 when N_j = 1 (H' bits, lower trees, lower N wires); a hosted
            # carrier is put back from H' to the raw new bit, which is 0 then. One X each.
            for w in carriers:
                f("x", w)
            for g in gates:
                f("ccx", *g)
            C.cx(t, vv)
            f("ccx", gate, vv, p)                # p: c_j -> N_(j+1); F's own part: 0 -> N_(j+1)
            later = set()
            for gate2, p2, carriers2, hosted2, gates2, vv2 in plan[j + 1:]:
                later |= set(carriers2) | {w for g in gates2 for w in g} | {gate2}
            if vv in later:
                C.cx(t, vv)
            else:
                deferred.append(vv)
        ntop = plan[-1][1] if plan else self.first
        C.ccx(t, ntop, blocks[-1][1])            # c_K ^= t*N_K
        for vv in deferred:
            C.cx(t, vv)
        # W3: F reversed.
        for op, w in reversed(F):
            (C.x if op == "x" else C.ccx)(*w)
        # free in the order the wires become idle: p_K first, then the tops W3 clears first
        C.free([blocks[-1][1]] + [b[1] for b in reversed(blocks[:-1])] + [self.first])


def sizes_for(h, dcap=31, top=30):
    """A block-size list for width h: today's doubling (each block as large as its carriers allow)
    up to a block of dcap bits, then sizes falling by one from `top` (at most one more than the
    block below, so every upper tree can sit on it) and never below top // 2, the last block taking
    what is left. At h = 221: 3, 7, 15, 31, 30, 29, 28, 27, 26, 23 (11 tops, one more than today)."""
    sizes, start, ntops = [], 2, 1
    while start < h:
        m = min(start + ntops, dcap, h - start)
        sizes.append(m); start += m; ntops += 1
        if m == dcap:
            break
    s, floor = top, max(2, top // 2)
    while start < h:
        m = max(1, min(s, h - start, start + ntops, (sizes[-1] + 1) if sizes else s))
        sizes.append(m); start += m; ntops += 1
        s = max(floor, s - 1)
    return sizes


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
        """Off J add a >> 1; on J rotate word's low n+1 bits left."""
        if CDKM_ADD and self.n > 3:
            return add_or_shift_cdkm(C, J, a, word, self.n)
        if WIDE_ADD and self.n > 3:
            sizes = _ramp_sizes(self.n - 1, WIDE_ADD, WIDE_ADD_RAMP)
            back = _ramp_sizes(self.n - 1, *WIDE_ADD_BACK) if WIDE_ADD_BACK else sizes
            if WIDE_ADD_FSIZES and sum(WIDE_ADD_FSIZES) == self.n - 1:
                sizes = list(WIDE_ADD_FSIZES)
            if WIDE_ADD_BSIZES and sum(WIDE_ADD_BSIZES) == self.n - 1:
                back = list(WIDE_ADD_BSIZES)
            temps = max(sum(r - 1 for r in sizes), sum(r - 1 for r in back))
            return add_or_shift_la(C, J, a, word, self.n, ntemps=max(temps, 1), ncopies=WIDE_ADD_COPIES,
                                   fsizes=sizes, bsizes=back)
        if LOOKAHEAD_ADD and self.n > 3:
            fsizes = _la_wide(self.n - 1, *LA_F_WIDE) if LA_F_WIDE else None
            bsizes = _la_wide(self.n - 1, *LA_B_WIDE) if LA_B_WIDE else None
            return add_or_shift_la(C, J, a, word, self.n, ntemps=LA_NTEMPS, ncopies=LA_NCOPIES,
                                   fsizes=fsizes, bsizes=bsizes)
        n = self.n
        (cin, D, carry, t) = C.alloc(4)
        extra = C.alloc(3)
        copies = [D] + extra
        sources, tree = [J], []
        filled = 0
        while filled < len(copies):
            layer = []
            for source in sources:
                if filled == len(copies):
                    break
                target = copies[filled]
                filled += 1
                C.cx(source, target)
                tree.append((source, target))
                layer.append(target)
            sources.extend(layer)
        C.x(J)
        ops, tgt = a[1:], word[:n - 1]
        carries = [cin] + ops[:-1]
        for i in range(n - 1):
            c, u = carries[i], tgt[i]
            _maj(C, c, u, ops[i])
            # x*(u ^ not-J ^ not-J*x) = x*u, so the ordinary undo still works.
            C.cx(J, u)
            C.cx(c, u)
            C.ccx(copies[i % len(copies)], c, u)
        C.ccx(J, ops[-1], carry)
        C.ccx(carry, word[n - 1], t)
        C.ccx(t, word[n], word[n + 1])
        C.cx(t, word[n])
        C.cx(carry, word[n - 1])
        C.x(word[n - 1])
        C.ccx(carry, word[n - 1], t)
        C.x(word[n - 1])
        C.ccx(J, ops[-1], carry)
        for i in reversed(range(n - 1)):
            c, u, o = carries[i], word[i], ops[i]
            C.ccx(c, u, o)
            C.cx(o, c)
            C.cx(J, u)
            C.cx(o, u)
            # Emit each reversal pair only after both bits are final.
            if i < (n + 1) // 2:
                cswap(C, copies[i % len(copies)], [word[i]], [word[n - i]])
            if i < n // 2:
                cswap(C, copies[(i + 1) % len(copies)], [word[i + 1]], [word[n - i]])
        C.x(J)
        for source, target in reversed(tree):
            C.cx(source, target)
        C.free(extra)
        C.free([cin, D, carry, t])

    def tick(self, C, prev, g, nxt, a, s, pebble=False, pin=None, pout=None):
        """One tick of a swept tape. prev is the cell before, or None on the first tick.
        nxt is the next cell's g0 wire, or None on the last tick, whose B1 reduces normally."""
        n, f, w = self.n, self.f, self.window
        tick_start = len(C.prims)
        io = [w for r in (prev or (), g, a, s) for w in r] + [w for w in (nxt, pin, pout) if w is not None]
        early_idle = _Idle(C, tick_start, io) if IDLE_FROM_START and IDLE_WIRES else None
        if early_idle:
            C.alloc = early_idle.alloc
        g0, b = g
        (spill, high) = C.alloc(2)
        word = s + [spill, high]
        # pin: the previous tick already left this tick's B1 mark on a wire.
        p = pin if pin is not None else self._after_b1(C, prev)
        if p is not None:
            cswap(C, p, [s[0]], [spill])
            cswap(C, p, [g0], [high])
            C.cx(p, g0)
        sw = self._begin(C, g0, b)
        if p is not None and WIDE_WAVE and P_FAN > 1:
            # p through copies on clean wires, so these CXs do not queue on one control; the add
            # no longer hides the queue behind a wave of its own.
            pc = [p] + C.alloc(4 * P_FAN - 1)
            ptree = [(pc[j - (1 << (j.bit_length() - 1))], pc[j]) for j in range(1, len(pc))]
            for src, dst in ptree:
                C.cx(src, dst)
            for i, x in enumerate(word):
                C.cx(pc[i * len(pc) // len(word)], x)
            for src, dst in reversed(ptree):
                C.cx(src, dst)
            C.free(pc[1:])
        elif p is not None:
            for x in word:
                C.cx(p, x)
        self._add_or_shift(C, b, a, word)
        add_end = len(C.prims)
        if p is not None:
            # The rotation finishes word's bits in mirror pairs from the middle out, so un-copy p in
            # that order. These CX all share the control p and commute, so the order is free.
            done, order = set(), [word[n + 1]]
            for i in reversed(range(n // 2 + 1)):
                order += [word[i], word[n - i]]
            order = [x for x in order + list(word) if not (x in done or done.add(x))]
            if P_FAN > 1:
                # p on clean wires the add just freed (the free stack's top, so no new wires), each
                # copy un-copying one contiguous slice, so the CX don't queue on one control.
                pc = [p] + C.alloc(P_FAN - 1)
                ptree = [(pc[j - (1 << (j.bit_length() - 1))], pc[j]) for j in range(1, len(pc))]
                for src, dst in ptree:
                    C.cx(src, dst)
                for i, x in enumerate(order):
                    C.cx(pc[i * len(pc) // len(order)], x)
                for src, dst in reversed(ptree):
                    C.cx(src, dst)
                C.free(pc[1:])
            else:
                for x in order:
                    C.cx(p, x)
        early_one = C.alloc(1) if ONE_BEFORE_DROP else None
        self._drop_after_b1(C, prev, p)
        (l0,) = [pout] if pout is not None else C.alloc(1)
        C.ccx(g0, a[0], l0)
        low = [l0] + word[:n - 1]
        L, H = low[:w], low[w:]
        k0, k1, k2 = word[n - 1:n + 2]

        # one = 0 only on a B1, which is the tick's inert lane. b holds J here, so this needs no CCX.
        # With ONE_BEFORE_DROP it is not the wire p's drop frees, which is busy until that ends.
        (one,) = early_one or C.alloc(1)
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
        idle = early_idle or (_Idle(C, tick_start) if IDLE_WIRES else None)
        if idle:
            C.alloc = idle.alloc
            C.alloc_after = idle.alloc_after
        if WIDE_P and len(H) > 1:
            P = wide.PrefixTree(C, H, WIDE_P, P_PEBBLE, P_FRESH - (P_FRESH_EDGE if pebble else 0))
        elif PARALLEL_PREFIX and len(H) > 1:
            P = PrefixR(C, H, sizes=list(P_SIZES_221) if len(H) == 221 else sizes_for(len(H)))
        else:
            P = Prefix(C, H, cap=37)
        M = wide.PrefixTree(C, L[ra:rb], WIDE_M, M_PEBBLE, M_FRESH) if WIDE_M and rb - ra > 1 else Prefix(C, L[ra:rb])
        if idle and not WIDE_IDLE:
            del C.alloc, C.alloc_after
        # Only M's skip reads M before its write, so the top below its last is cleared until then.
        M.release()

        # With at most one H bit, P reads a data wire or nothing, and there is nothing to borrow.
        spare = H[:] if len(H) > 1 and not WIDE_P else []

        def alloc():
            (c,) = C.alloc(1)
            return c

        def take():
            if not spare:
                return alloc()
            # Any H wire reads 1 when P's top is 1, so which one is borrowed is free.
            c = spare.pop(idle.first(spare) if idle else -1)
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
        aside = WIDE_LOW and WIDE_ASIDE
        if aside:
            # The trial chain runs on copies of its selectors and of L's bits outside the quiet run,
            # so the second chain need not wait for its undo: the selectors change and the second
            # chain starts as soon as e is known, and the trial unwinds beside it.
            act = [i for i in range(w) if not ra <= i < rb]
            lc, sc = C.alloc(len(act)), C.alloc(len(selectors))
            trial = list(L)
            for i, c in zip(act, lc):
                C.cx(L[i], c)
                trial[i] = c
            for s_, c in zip(selectors, sc):
                C.cx(s_, c)
            first = wide.LowAdd(C, trial, sc, masks, self.run, M.top, WIDE_LOW)
            carry_out = first.top
        elif WIDE_LOW:
            first = wide.LowAdd(C, L, selectors, masks, self.run, M.top, WIDE_LOW)
            carry_out = first.top
        else:
            chain = self._chain_up(C, L, selectors, masks, M, {}, take)
            carry_out = chain[w]
        (e,) = C.alloc(1)
        if H:
            C.ccx(P.top, carry_out, e)
        else:
            C.cx(carry_out, e)
        if aside:
            pass
        elif WIDE_LOW:
            first.undo()
        else:
            self._chain_down(C, L, selectors, masks, M, {}, chain, drop, False)
        if nxt is None:
            C.x(one)
        else:
            C.cx(b, one)
            C.cx(g0, one)
        if not aside:
            C.cx(sw, b)
            self._end(C, g0, b, sw)
        if KUPD_FAST and not spare:
            # k' = k + e, and the second chain's corner selector is k0' k1' = k0 k1 ^ e k1 in the
            # old bits: three holds k0 k1, so one CCX turns it into the corner selector, and three's
            # wire becomes one. one's own wire holds 0 here and is freed instead of three's.
            C.ccx(three, e, k2)
            C.ccx(k1, e, three)
            C.ccx(k0, e, k1)
            C.cx(e, k0)
            C.free([one])
            one = three
        else:
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
        if WIDE_LOW:
            def clear_sc():
                one_c, k0_c, k1_c, three_c = sc
                # The copies hold the old one, k0, k1 and k0 k1; k1 took k0 e and k0 took e since.
                C.ccx(k0_c, k1_c, three_c)
                C.cx(k1, k1_c)
                C.ccx(k0_c, e, k1_c)
                C.cx(k0, k0_c)
                C.cx(e, k0_c)
                if nxt is None:
                    C.x(one_c)
                else:
                    C.cx(b, one_c)
                    C.cx(g0, one_c)
                C.free(sc)

            if aside and ASIDE_UNDO_FIRST:
                # The trial on copies unwinds beside the second chain; emitted first, its first
                # reverse join reads M's top before the second chain's skip does, not after.
                first.undo()
            if aside and SECOND_FRESH and idle:
                # The trial's unwinding has just freed its wires, busy until it ends: the second
                # chain takes wires this tick has not touched, and P's write later takes those.
                C.alloc = lambda k=1: [idle._untouched() for _ in range(k)]
            second = wide.LowAdd(C, L, selectors, masks, self.run, M.top, WIDE_LOW)
            if aside and SECOND_FRESH and idle:
                C.alloc = idle.alloc
            carry_out = second.top
            if aside:
                if not ASIDE_UNDO_FIRST:
                    first.undo()
                clear_sc()
                # Until its write the second chain leaves L ^ K on each of these bits of L.
                for i, c in zip(act, lc):
                    C.cx(L[i], c)
                    for s_ in second._on(i):
                        C.cx(s_, c)
                C.free(lc)
                C.cx(sw, b)
                self._end(C, g0, b, sw)
        else:
            perch = self._perches(selectors, masks)
            kept = self._kept_carries(perch)
            schedule = self._pebble_schedule(perch) if pebble else None
            moves = None
            if schedule is not None:
                moves, low_kept = schedule
                kept = set(low_kept) | {self.run[0]}
            chain = self._chain_up(C, L, selectors, masks, M, perch, alloc, kept, moves)
            carry_out = chain[w]
        if aside and isinstance(P, wide.PrefixTree):
            # e = carry_out * AND(old H) is also carry_out * AND(~new H), which P's write holds
            # before its last use of carry_out. Clearing e there keeps this CCX, and the trial's
            # unwinding before it, out of the write's way.
            P.write(carry_out, lambda z: (C.ccx(carry_out, z, e), drop(e)))
        else:
            if H:
                C.ccx(P.top, carry_out, e)
            else:
                C.cx(carry_out, e)
            drop(e)
            P.write(carry_out)
        if idle and RESTORE_IDLE:
            # M's restore would take the wire chain[b]'s drop has just freed; take the one free first.
            M.alloc_restore = lambda: idle.alloc(1)[0]
        if WIDE_LOW and isinstance(M, wide.PrefixTree):
            second.write(lambda c, mid: M.write(c, mid))
        elif WIDE_LOW:
            second.write(lambda c, mid: (mid(M.top), M.restore(), M.write(c)))
        else:
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
        # Carry-select blocks borrow only clean wires already touched by this tick.
        cmp_a, cmp_b = a[1:], low[1:]
        width = len(cmp_a)
        if width and WIDE_CMP:
            wide.gt_tree(C, cmp_a, cmp_b, k1, z, WIDE_CMP, split=CMP_SPLIT)
        elif width:
            (cin,) = C.alloc(1)
            # Each wire's last op in this tick so far: the compare hands the top block the wires that
            # came free first. Touched free wires read 0.
            last = {}
            for pos in range(tick_start, len(C.prims), 4):
                last[C.prims[pos + 1]] = last[C.prims[pos + 2]] = last[C.prims[pos + 3]] = pos
            saved_free = list(C.free_stack)
            scratch = [wire for wire in saved_free if wire in last]
            # Wires the compare may borrow at any value; it restores them and nothing else runs meanwhile.
            dirty = [w for w in (a[0], nxt, b, l0, *(prev or ()), p, k0, g0, k2) if w is not None]
            C.free_stack[:] = [wire for wire in saved_free if wire not in last]
            compare_blocks(C, cmp_a, cmp_b, z, k1, cin, scratch, dirty, last)
            C.free_stack[:] = saved_free
            C.free([cin])
        mark_even(z)
        C.free([z])
        if idle and WIDE_IDLE:
            del C.alloc, C.alloc_after
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
        if LATE_FINAL_WIRES:
            # The next tick maps its scratch as this one does, so wires this tick first touched
            # after its add are ones the next tick also reaches only late: put sw there.
            first_touch = {}
            for pos in range(tick_start, len(C.prims), 4):
                for w in C.prims[pos + 1:pos + 4]:
                    first_touch.setdefault(w, pos)
            late = [w for w in C.free_stack if first_touch.get(w, -1) >= add_end]
            if late:
                C.free_stack.remove(late[-1])
                C.free_stack.append(late[-1])
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
        if pout is not None:
            # l0 is clean again and g0, b are final: leave the next tick's mark, NOR(g0, b), on it.
            C.x(l0)
            C.cx(g0, l0)
            C.cx(b, l0)
            C.ccx(g0, b, l0)
            C.free([spill, high])
        else:
            C.free([l0, spill, high])
        s[:] = s[-1:] + s[:-1]
        # Reuse only this tick's clean scratch, keeping its recorded arity fixed.
        touched = set(C.prims[tick_start + 1::4])
        touched.update(C.prims[tick_start + 2::4])
        touched.update(C.prims[tick_start + 3::4])
        saved_free = list(C.free_stack)
        copies = sorted(touched.intersection(saved_free))
        if LATE_FINAL_WIRES:
            # Only wires the next tick first reaches after its add, so its start doesn't wait on them.
            first_touch = {}
            for pos in range(tick_start, len(C.prims), 4):
                for w in C.prims[pos + 1:pos + 4]:
                    first_touch.setdefault(w, pos)
            late = [w for w in copies if first_touch[w] >= add_end]
            copies = late if len(late) >= LATE_MIN_COPIES else copies
        C.free_stack[:] = [wire for wire in saved_free if wire not in copies]
        sources, tree = [sw], []
        filled = 0
        pairs = list(zip(a, s))
        # The next tick starts on bit 0's swap: run it on sw before sw's fan, which waits on
        # copy wires the compare frees late. Both only read sw, so they commute.
        for u, v in pairs[:FINAL_EARLY]:
            C.cx(v, u)
            C.ccx(sw, u, v)
            C.cx(v, u)
        while filled < len(copies):
            layer = []
            for source in sources:
                if filled == len(copies):
                    break
                target = copies[filled]
                filled += 1
                C.cx(source, target)
                tree.append((source, target))
                layer.append(target)
            sources.extend(layer)
        for j, (u, v) in enumerate(pairs):
            if j < FINAL_EARLY:
                continue
            C.cx(v, u)
            C.ccx(sources[j % len(sources) if FINAL_INTERLEAVE else j * len(sources) // n], u, v)
            C.cx(v, u)
        for source, target in reversed(tree):
            C.cx(source, target)
        C.free_stack[:] = saved_free
        self._end(C, g0, b, sw)
