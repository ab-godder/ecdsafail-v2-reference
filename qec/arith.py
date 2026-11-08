"""Exact reversible integer and modular arithmetic. Registers are wire lists,
least significant bit first. Every routine returns its ancillas clean on every
input inside its stated domain."""


def _maj(C, c, b, a):
    C.cx(a, b)
    C.cx(a, c)
    C.ccx(c, b, a)


def _uma(C, c, b, a):
    C.ccx(c, b, a)
    C.cx(a, c)
    C.cx(c, b)


def _cuma(C, ctrl, c, b, a):
    C.ccx(c, b, a)
    C.ccx(ctrl, c, b)
    C.cx(a, c)
    # Left ungated because it restores b when ctrl is 0, so _cuma must not become a gated _uma.
    C.cx(a, b)


def _cx(C, ctrl, src, dst):
    if ctrl is None:
        C.cx(src, dst)
    else:
        C.ccx(ctrl, src, dst)


def add(C, a, b, cout=None, ctrl=None, cin=None):
    """b += ctrl*a mod 2^w (CDKM). cout ^= ctrl*carry out. A wire holding 1 as cin adds 1."""
    w = len(a)
    assert len(b) == w
    assert ctrl is None or ctrl not in (*a, *b, cout)
    if w == 0:
        return
    own = cin is None
    if own:
        (cin,) = C.alloc(1)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    if cout is not None:
        _cx(C, ctrl, a[-1], cout)
    for i in reversed(range(w)):
        if ctrl is None:
            _uma(C, carries[i], b[i], a[i])
        else:
            _cuma(C, ctrl, carries[i], b[i], a[i])
    if own:
        C.free([cin])


def sub(C, a, b, borrow=None, ctrl=None):
    """b -= ctrl*a mod 2^w. borrow ^= ctrl*[a > b]."""
    for w in b:
        C.x(w)
    add(C, a, b, borrow, ctrl)
    for w in b:
        C.x(w)


def gt_into(C, a, b, t, ctrl=None):
    """t ^= ctrl*[a > b], a and b restored. Carry out of (~b + a) is [a > b]."""
    assert ctrl is None or ctrl not in (*a, *b, t)
    w = len(a)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    for x in b:
        C.x(x)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    _cx(C, ctrl, a[-1], t)
    for i in reversed(range(w)):
        C.ccx(carries[i], b[i], a[i])
        C.cx(a[i], carries[i])
        C.cx(a[i], b[i])
    for x in b:
        C.x(x)
    C.free([cin])


def _gt_up(C, cin, a, b):
    """a[-1] = carry out of (b + a + cin), for b already complemented; _gt_down undoes it."""
    carries = [cin] + a[:-1]
    for i in range(len(a)):
        _maj(C, carries[i], b[i], a[i])


def _gt_down(C, cin, a, b):
    carries = [cin] + a[:-1]
    for i in reversed(range(len(a))):
        C.ccx(carries[i], b[i], a[i])
        C.cx(a[i], carries[i])
        C.cx(a[i], b[i])


def gt_into_blocked(C, a, b, t, ctrl, sizes, zeros=(), dirty=()):
    """t ^= ctrl*[a > b], a and b restored: gt_into's map, with the words cut into blocks of the
    given sizes, low first, whose ripples run side by side.

    Block 0 ripples once with carry in 0, so its top a wire holds the carry into block 1. Each
    other block j ripples with a carry in c, copies its carry out v1 onto a wire G, ripples back,
    and holds the carry out v2 for carry in not c. With g = [a_j > b_j] and e = [a_j == b_j],
    v1 = g ^ e*c and v2 = g ^ e*(1 - c), so for any c the carry into block j + 1 is
    MAJ(v1, v2, carry into j), one more ripple over the blocks onto G that ends in t. All of it is
    then undone, and a third ripple clears G.

    Block 0's carry in and each G read 0: a fresh wire, or one of zeros, idle wires that read 0
    whenever ctrl is 1. A higher block's carry in may be one of dirty, idle wires of any value,
    taken bottom block first. Every wire but t is restored on every input.
    """
    m = len(sizes)
    assert ctrl not in (*a, *b, t, *zeros, *dirty)
    assert sum(sizes) == len(a) == len(b) and min(sizes) > 0
    dirty = list(dirty)[:m - 1]
    zeros = list(zeros)[:2 * m - 1 - len(dirty)]
    fresh = C.alloc(2 * m - 1 - len(dirty) - len(zeros))
    # Roles in the order they are first needed take the wires idle longest: the allocator hands
    # out the last freed first, and zeros come last.
    clean = fresh[::-1] + zeros
    cins = [None] + [dirty.pop(0) if dirty else clean.pop(0) for _ in range(1, m)]
    cins[0] = clean.pop(0)
    gs = [None] + clean
    blocks, lo = [], 0
    for s in sizes:
        blocks.append((a[lo:lo + s], b[lo:lo + s]))
        lo += s
    tops = [a_j[-1] for a_j, _ in blocks]
    for x in b:
        C.x(x)
    for j in range(1, m):
        _gt_up(C, cins[j], *blocks[j])
        C.cx(tops[j], gs[j])
        _gt_down(C, cins[j], *blocks[j])
        C.x(cins[j])
        _gt_up(C, cins[j], *blocks[j])
    _gt_up(C, cins[0], *blocks[0])
    carries = tops[:1] + gs[1:]
    for j in range(1, m):
        _maj(C, carries[j - 1], tops[j], gs[j])
    _cx(C, ctrl, carries[-1], t)
    for j in reversed(range(1, m)):
        C.ccx(carries[j - 1], tops[j], gs[j])
        C.cx(gs[j], carries[j - 1])
        C.cx(gs[j], tops[j])
    _gt_down(C, cins[0], *blocks[0])
    for j in range(1, m):
        _gt_down(C, cins[j], *blocks[j])
        C.x(cins[j])
        _gt_up(C, cins[j], *blocks[j])
        C.cx(tops[j], gs[j])
        _gt_down(C, cins[j], *blocks[j])
    for x in b:
        C.x(x)
    C.free(fresh)


def carries_blocked(C, a, b, outs, ctrl, sizes):
    """outs[j] ^= ctrl*(carry out of the low blocks 0..j of b + a), a and b restored.

    gt_into_blocked's construction, read at every block boundary instead of only the top: block 0
    ripples with carry in 0, each higher block j ripples with carry in 0 and 1 and keeps both
    carry outs v1, v2, and a ripple of MAJ(v1, v2, carry into j) over the blocks gives each
    boundary's carry, which is read onto outs[j] right after its MAJ. The undo skips the reads.
    outs[-1] gets the carry out of the whole word. Every wire but outs is restored."""
    m = len(sizes)
    assert len(outs) == m and sum(sizes) == len(a) == len(b) and min(sizes) > 0
    assert ctrl not in (*a, *b, *outs)
    fresh = C.alloc(2 * m - 1)
    cins = [fresh[0]] + fresh[1:m]
    gs = [None] + fresh[m:]
    blocks, lo = [], 0
    for s in sizes:
        blocks.append((a[lo:lo + s], b[lo:lo + s]))
        lo += s
    tops = [a_j[-1] for a_j, _ in blocks]
    for j in range(1, m):
        _gt_up(C, cins[j], *blocks[j])
        C.cx(tops[j], gs[j])
        _gt_down(C, cins[j], *blocks[j])
        C.x(cins[j])
        _gt_up(C, cins[j], *blocks[j])
    _gt_up(C, cins[0], *blocks[0])
    _cx(C, ctrl, tops[0], outs[0])
    carries = tops[:1] + gs[1:]
    for j in range(1, m):
        _maj(C, carries[j - 1], tops[j], gs[j])
        _cx(C, ctrl, gs[j], outs[j])
    for j in reversed(range(1, m)):
        C.ccx(carries[j - 1], tops[j], gs[j])
        C.cx(gs[j], carries[j - 1])
        C.cx(gs[j], tops[j])
    _gt_down(C, cins[0], *blocks[0])
    for j in range(1, m):
        _gt_down(C, cins[j], *blocks[j])
        C.x(cins[j])
        _gt_up(C, cins[j], *blocks[j])
        C.cx(tops[j], gs[j])
        _gt_down(C, cins[j], *blocks[j])
    C.free(fresh)


def load_const(C, reg, value, ctrl=None):
    for i, w in enumerate(reg):
        if value >> i & 1:
            if ctrl is None:
                C.x(w)
            else:
                C.cx(ctrl, w)


def _carry_into(C, b, bit, c, t, ctrl=None):
    """t ^= MAJ(b, ctrl*bit, c), where c is a wire or None for 0. With ctrl, c
    must imply ctrl, which holds for every carry of reg + ctrl*value."""
    if not bit:
        if c is not None:
            C.ccx(b, c, t)
    elif c is None:
        _cx(C, ctrl, b, t)
    elif ctrl is None:
        C.ccx(b, c, t)
        C.cx(b, t)
        C.cx(c, t)
    else:
        # MAJ(b, ctrl, c) = c ^ b*(ctrl ^ c) when c implies ctrl.
        C.cx(c, ctrl)
        C.ccx(b, ctrl, t)
        C.cx(c, ctrl)
        C.cx(c, t)


def _carries(C, reg, value, ctrl=None):
    """cs[i] is the carry into reg[i] of reg + ctrl*value. Each entry is None
    when the carry is always 0, reg[i - 1] itself when it equals that bit, or a
    fresh wire. Release the fresh wires with _uncarry."""
    cs = [None]
    for i in range(len(reg) - 1):
        bit, c = value >> i & 1, cs[i]
        if c is None and not bit:
            cs.append(None)
        elif c is None and ctrl is None:
            cs.append(reg[i])
        else:
            (t,) = C.alloc(1)
            _carry_into(C, reg[i], bit, c, t, ctrl)
            cs.append(t)
    return cs


def _uncarry(C, reg, value, cs, i, ctrl=None):
    """Clear and free cs[i]. reg[i - 1] and cs[i - 1] must still hold what they held in _carries."""
    if i and cs[i] is not None and cs[i] != reg[i - 1]:
        _carry_into(C, reg[i - 1], value >> (i - 1) & 1, cs[i - 1], cs[i], ctrl)
        C.free([cs[i]])


def add_const(C, reg, value, ctrl=None, cout=None):
    """reg += ctrl*value mod 2^w, cout ^= ctrl*carry out. The carry chain reads
    the constant's bits, so 0 bits below the first 1 cost nothing."""
    w = len(reg)
    value %= 1 << w
    assert ctrl is None or ctrl not in (*reg, cout)
    cs = _carries(C, reg, value, ctrl)
    if cout is not None:
        _carry_into(C, reg[-1], value >> (w - 1) & 1, cs[-1], cout, ctrl)
    for i in reversed(range(w)):
        if value >> i & 1:
            if ctrl is None:
                C.x(reg[i])
            else:
                C.cx(ctrl, reg[i])
        if cs[i] is not None:
            C.cx(cs[i], reg[i])
        _uncarry(C, reg, value, cs, i, ctrl)


def gt_const_into(C, value, b, t, ctrl=None):
    """t ^= ctrl*[value > b], b restored. Carry out of (~b + value) is [value > b]."""
    assert ctrl is None or ctrl not in (*b, t)
    assert t not in b
    assert value < 1 << len(b)
    for x in b:
        C.x(x)
    cs = _carries(C, b, value, ctrl)
    _carry_into(C, b[-1], value >> (len(b) - 1) & 1, cs[-1], t, ctrl)
    for i in reversed(range(len(b))):
        _uncarry(C, b, value, cs, i, ctrl)
    for x in b:
        C.x(x)


def _reduce(C, r, value):
    """r holds s on m wires, value has bit m - 1 set, and s's top bit and the
    carry into it are never both 1. Leaves (s + flag*value) mod 2^(m-1) in
    r[:-1] and flag = [s + value >= 2^m] in r[-1]. One carry chain serves both
    the comparison and the conditional add."""
    assert value >> (len(r) - 1) & 1
    cs = _carries(C, r, value)
    flag = r[-1]
    if cs[-1] is not None:
        C.cx(cs[-1], flag)
    for i in reversed(range(len(r) - 1)):
        _uncarry(C, r, value, cs, i + 1)
        if cs[i] is not None:
            C.ccx(flag, cs[i], r[i])
        if value >> i & 1:
            C.cx(flag, r[i])


def cswap(C, ctrl, a, b):
    for ai, bi in zip(a, b):
        C.cx(bi, ai)
        C.ccx(ctrl, ai, bi)
        C.cx(bi, ai)


def crotate(C, ctrls, a):
    """Under ctrls[0]: a[i] <- a[i + 1] for i < len(a) - 1 and a[-1] <- a[0]. ctrls[1:] are
    clean wires that carry copies of the control while the swaps run, clean again after.

    The rotation is a reversal of a followed by a reversal of a[:-1]. Each reversal is a layer
    of disjoint swaps, so no swap waits on a chain through a or on one shared control: the
    layers are interleaved low pairs first, so low wires settle first, and dealt round robin
    to the copies."""
    m, k = len(a), len(ctrls)
    first = [(i, m - 1 - i) for i in range(m // 2)]
    second = [(i, m - 2 - i) for i in range((m - 1) // 2)]
    # second[i] reads the wires first[i] and first[i + 1] wrote, and no other first swap.
    order = first[:1]
    for i in range(1, len(first)):
        order += [first[i], second[i - 1]]
    order += second[max(len(first) - 1, 0):]
    # Copy j is made from copy j minus its top bit: a fan-out tree of depth log2(k).
    fan = [(ctrls[j - (1 << (j.bit_length() - 1))], ctrls[j]) for j in range(1, k)]
    for src, dst in fan:
        C.cx(src, dst)
    for t, (i, j) in enumerate(order):
        cswap(C, ctrls[t % k], [a[i]], [a[j]])
    for src, dst in reversed(fan):
        C.cx(src, dst)


def _selected_carry(C, bit, selectors, masks, i, incoming, target):
    """XOR a majority carry into target, restoring bit and the selectors.

    The selected addend bit is an XOR of constant masks. Callers use exclusive
    selectors, so that XOR is also the intended selected integer's bit.
    """
    active = [s for s, mask in zip(selectors, masks) if mask >> i & 1]
    if not active:
        if incoming is not None:
            C.ccx(bit, incoming, target)
        return
    b = active[0]
    for other in active[1:]:
        C.cx(other, b)
    if incoming is None:
        C.ccx(bit, b, target)
    else:
        # MAJ(x,b,c) = c XOR ((x XOR c) AND (b XOR c)).
        C.cx(incoming, bit)
        C.cx(incoming, b)
        C.ccx(bit, b, target)
        C.cx(incoming, target)
        C.cx(incoming, b)
        C.cx(incoming, bit)
    for other in reversed(active[1:]):
        C.cx(other, b)


def _selected_chain(C, low, selectors, masks):
    carries = [None]
    for i, bit in enumerate(low):
        (carry,) = C.alloc(1)
        _selected_carry(C, bit, selectors, masks, i, carries[-1], carry)
        carries.append(carry)
    return carries


def _add_selected_const(C, b, selectors, masks):
    """Add one of the constant masks under exclusive selectors, inline."""
    n = len(b)
    masks = [mask % (1 << n) for mask in masks]
    start = min(((mask & -mask).bit_length() - 1 for mask in masks if mask), default=n)
    low = b[start:]
    shifted = [mask >> start for mask in masks]
    chain = _selected_chain(C, low[:-1], selectors, shifted)
    for i in reversed(range(len(low))):
        if i + 1 < len(low):
            # Clear the carry against the old bit, before writing the sum.
            _selected_carry(C, low[i], selectors, shifted, i, chain[i], chain[i + 1])
            C.free([chain[i + 1]])
        for selector, mask in zip(selectors, shifted):
            if mask >> i & 1:
                C.cx(selector, low[i])
        if chain[i] is not None:
            C.cx(chain[i], low[i])


class Mod:
    """Arithmetic modulo an odd q < 2^n on canonical values in [0, q)."""

    def __init__(self, q):
        self.q = q
        self.n = q.bit_length()
        self.low_bits = (2 * ((1 << self.n) - q)).bit_length()

    def _wrap(self):
        """2^(n+1) - q, so s + _wrap() carries out of n + 1 bits exactly when s >= q."""
        return (1 << (self.n + 1)) - self.q

    def modadd(self, C, a, b, ctrl=None):
        """b = (b + ctrl*a) mod q, for b < q, and a < q when ctrl is 1."""
        (top,) = C.alloc(1)
        add(C, a, b, top, ctrl)
        # s = b + ctrl*a < 2q, so s's top bit and the carry into it are never both 1.
        _reduce(C, b + [top], self._wrap())
        # top = [q was subtracted] = ctrl*[a > b]; ctrl = 0 means s = b < q and top = 0.
        gt_into(C, a, b, top, ctrl)
        C.free([top])

    def moddbl(self, C, s):
        """s = 2s mod q. 2s is a relabel with a zero low bit, so the reduction
        runs on s against _wrap() >> 1, and the low bit of the result is the flag."""
        _reduce(C, s, self._wrap() >> 1)
        s[:] = s[-1:] + s[:-1]

    def dbl_add(self, C, ctrl, a, s, m=None):
        """s = (2s + ctrl*a) mod q, canonical inputs, with clean scratch.

        Like moddbl, rotate the output wire list right by one. The two quotient
        bits are cleared by parity and a controlled comparison, respectively.
        See DESIGN.md for the all-input reduction and cleanup argument.
        """
        n, q = self.n, self.q
        m = self.low_bits if m is None else m
        assert q >= 3 and q & 1 and self.low_bits <= m <= n
        assert len(a) == len(s) == n and len(set([ctrl] + a + s)) == 2 * n + 1
        f = (1 << n) - q
        l0, top = C.alloc(2)
        add(C, a[1:] + [l0], s, top, ctrl=ctrl)
        C.ccx(ctrl, a[0], l0)
        low = [l0] + s[:-1]
        L, H = low[:m], low[m:]
        # A[1] aliases H[0]; only A[2:] are owned scratch.
        prefix = [None] + (H[:1] if H else [])
        for bit in H[1:]:
            (p,) = C.alloc(1)
            C.ccx(prefix[-1], bit, p)
            prefix.append(p)
        k0, k1 = s[-1], top
        C.x(k1)
        C.cx(k0, k1)
        selectors, masks = [k1, k0], [f, 2 * f]
        chain = _selected_chain(C, L, selectors, masks)
        (e,) = C.alloc(1)
        if H:
            C.ccx(prefix[-1], chain[-1], e)
        else:
            C.cx(chain[-1], e)
        for i in reversed(range(m)):
            _selected_carry(C, L[i], selectors, masks, i, chain[i], chain[i + 1])
            C.free([chain[i + 1]])
        C.cx(k0, k1)
        C.x(k1)
        C.ccx(k0, e, k1)
        C.cx(e, k0)
        selectors = [k0, k1]
        chain = _selected_chain(C, L, selectors, masks)
        if H:
            C.ccx(prefix[-1], chain[-1], e)
        else:
            C.cx(chain[-1], e)
        C.free([e])
        for i in reversed(range(1, len(H))):
            # Clear the prefix before overwriting the bit that computed it.
            C.ccx(prefix[i], H[i], prefix[i + 1])
            C.ccx(chain[-1], prefix[i], H[i])
            C.free([prefix[i + 1]])
        if H:
            C.cx(chain[-1], H[0])
        for i in reversed(range(m)):
            _selected_carry(C, L[i], selectors, masks, i, chain[i], chain[i + 1])
            C.free([chain[i + 1]])
            if i:
                for selector, mask in zip(selectors, masks):
                    if mask >> i & 1:
                        C.cx(selector, L[i])
                C.cx(chain[i], L[i])
        (z,) = C.alloc(1)
        C.x(k0)
        C.ccx(ctrl, k0, z)
        C.x(k0)
        gt_into(C, a[1:], low[1:], k1, ctrl=z)
        C.x(k0)
        C.ccx(ctrl, k0, z)
        C.x(k0)
        C.cx(l0, k0)
        C.ccx(ctrl, a[0], l0)
        C.free([l0, top, z])
        s[:] = s[-1:] + s[:-1]

    def cmodadd_const(self, C, ctrl, value, b):
        """b = (b + ctrl*value) mod q on canonical b, with ctrl restored."""
        value %= self.q
        if value == 0:
            return
        (k,) = C.alloc(1)
        gt_const_into(C, self.q - value, b, k, ctrl)
        C.cx(ctrl, k)
        C.cx(k, ctrl)
        _add_selected_const(C, b, [ctrl, k], [value, value + (1 << self.n) - self.q])
        C.cx(k, ctrl)
        gt_const_into(C, value, b, k, ctrl)
        C.free([k])

    def creflect(self, C, ctrl, a, b):
        """b = (a - b) mod q when ctrl, valid on all canonical b."""
        a %= self.q
        (k,) = C.alloc(1)
        gt_const_into(C, a + 1, b, k, ctrl)
        C.cx(ctrl, k)
        for wire in b:
            C.cx(ctrl, wire)
        C.cx(k, ctrl)
        _add_selected_const(C, b, [ctrl, k], [a + 1, (a + 1 + self.q) % (1 << self.n)])
        C.cx(k, ctrl)
        C.cx(ctrl, k)
        gt_const_into(C, a + 1, b, k, ctrl)
        C.free([k])

    def cneg(self, C, ctrl, b):
        """b = q - b when ctrl, valid for b in [1, q). A bijection on all of
        {0,1}^n either way, since it is ~(b + 2^n - q - 1) mod 2^n."""
        add_const(C, b, (1 << self.n) - self.q - 1, ctrl=ctrl)
        for w in b:
            C.cx(ctrl, w)
