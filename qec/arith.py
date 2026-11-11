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
    # With no carry out and no control the top _maj and _uma cancel to b ^= a ^ c.
    top = w - 1 if cout is None and ctrl is None else w
    for i in range(top):
        _maj(C, carries[i], b[i], a[i])
    if cout is not None:
        _cx(C, ctrl, a[-1], cout)
    if top < w:
        C.cx(a[-1], b[-1])
        C.cx(carries[-1], b[-1])
    for i in reversed(range(top)):
        if ctrl is None:
            _uma(C, carries[i], b[i], a[i])
        else:
            _cuma(C, ctrl, carries[i], b[i], a[i])
    if own:
        C.free([cin])


def lookahead_sizes(w):
    """Block sizes for a w-bit lookahead chain, from the bottom: 1, 1, 2, 2, ... up to 8, and the
    same from the top down, so the blocks at both ends are short. A long block's precompute holds
    up the chain at the bottom, and its fixes and clears hold up the end of the sweep at the top."""
    bottom, top, left, b, d = [], [], w, 0, 0
    while left > 0:
        nb, nt = min(8, 1 + b // 2), min(8, 1 + d // 2)
        if nb <= nt:
            r = min(nb, left)
            bottom.append(r)
            b += 1
        else:
            r = min(nt, left)
            top.append(r)
            d += 1
        left -= r
    return bottom + top[::-1]


def lookahead_chain(C, U, V, temps, sizes):
    """f_j = x_j ^ p_j f_(j-1) in place on U for j = 1 .. len(U) - 1, from the chain input U[0],
    with U[j] holding x_j and V[j] holding p_j (a carry chain in CDKM's up-sweep form).
    For a block of r bits after U[i], first X_(i+j) = x_(i+j) ^ p_(i+j) X_(i+j-1) goes onto U[i+j]
    and the products T_j = p_(i+1) ... p_(i+j) onto temps, with no carry; then one CCX in series,
    U[i+r] ^= T_r U[i]; then U[i+j] ^= T_j U[i] for j < r, each temp cleared after its last use.
    One CCX a block is in series instead of one a bit. Temps are taken in turn and left clean; a
    block never needs more than len(temps) of them at once."""
    n = len(U)
    q = 0
    i = 0
    blk = 0
    while i < n - 1:
        r = min(sizes[min(blk, len(sizes) - 1)], len(temps) + 1, n - 1 - i)
        blk += 1
        if r == 1:
            C.ccx(U[i], V[i + 1], U[i + 1])
            i += 1
            continue
        for j in range(2, r + 1):
            C.ccx(U[i + j - 1], V[i + j], U[i + j])
        T = {1: V[i + 1]}
        for j in range(2, r + 1):
            T[j] = temps[q % len(temps)]
            q += 1
            C.ccx(T[j - 1], V[i + j], T[j])
        C.ccx(U[i], T[r], U[i + r])
        for j in range(r, 1, -1):
            C.ccx(T[j - 1], V[i + j], T[j])
            C.ccx(U[i], T[j - 1], U[i + j - 1])
        i += r


class _Gates:
    """Records cx and ccx calls, so a sweep can be emitted forward and then in reverse."""

    def __init__(self):
        self.ops = []

    def cx(self, a, b):
        self.ops.append((a, b))

    def ccx(self, a, b, c):
        self.ops.append((a, b, c))


def _emit(C, ops):
    for o in ops:
        (C.cx if len(o) == 2 else C.ccx)(*o)


def adjacent_xor(C, a):
    """a[i] ^= a[i + 1] for i < len(a) - 1, from the old values, in about 2 log2 len(a) CX layers:
    the inverse of an in-place suffix XOR done as a Brent-Kung prefix network (Brent and Kung,
    1982), so that network's CXs in reverse order. The ripple it replaces is len(a) - 1 CXs in
    series, since each shares a wire with the next."""
    n = len(a)
    ops = []
    d = 1
    while d < n:
        for p in range(2 * d - 1, n, 2 * d):
            ops.append((p - d, p))
        d *= 2
    while d > 1:
        d //= 2
        for p in range(3 * d - 1, n, 2 * d):
            ops.append((p - d, p))
    # prefix position p is suffix index n - 1 - p: a[n - 1 - p] ^= a[n - 1 - p + d]
    for s, t in reversed(ops):
        C.cx(a[n - 1 - s], a[n - 1 - t])


def prefix_xor(C, a):
    """a[i] ^= a[0] ^ ... ^ a[i - 1] for every i, from the old values, as a Brent-Kung prefix network
    (Brent and Kung, 1982) in about 2 log2 len(a) CX layers, instead of len(a) - 1 CXs in series."""
    n = len(a)
    d = 1
    while d < n:
        for p in range(2 * d - 1, n, 2 * d):
            C.cx(a[p - d], a[p])
        d *= 2
    while d > 1:
        d //= 2
        for p in range(3 * d - 1, n, 2 * d):
            C.cx(a[p - d], a[p])


class _TempQueue:
    """Clean temps handed out in turn; a freed temp goes to the back, so it is reused as late as
    possible. Counts the most taken at once."""

    def __init__(self, temps):
        self.free, self.live, self.peak = list(temps), 0, 0

    def take(self):
        self.live += 1
        self.peak = max(self.peak, self.live)
        return self.free.pop(0)

    def give(self, t):
        self.live -= 1
        self.free.append(t)


def lookahead_tree(C, U, V, temps, r):
    """f_j = x_j ^ p_j f_(j-1) in place on U for j = 1 .. len(U) - 1 from the chain input U[0],
    with U[j] holding x_j and V[j] holding p_j, as lookahead_chain, but with the blocks joined by a
    Brent-Kung prefix tree (Brent and Kung, 1982; in place as in Draper et al., quant-ph/0406142)
    instead of a chain: about 2 log2 of the block count CCX layers in series, for
    tree_temps(len(U), r) temps. Blocks of r bits: local prefix in place and the block's product on
    a temp (inner products cleared and made again for the fixes when r >= 3); up-sweep G_R ^=
    P_R G_L with each segment's product on a new temp unless the segment starts at the input;
    down-sweep; products cleared in reverse; then each block's fixes U[i+j] ^= T_j U[i]."""
    _lookahead_tree(C, U, V, _TempQueue(temps), r)


def tree_temps(n, r):
    """The most temps lookahead_tree holds at once for len(U) = n and blocks of r."""
    q = _TempQueue(range(4 * n))
    _lookahead_tree(_Gates(), list(range(n)), list(range(n, 2 * n)), q, r)
    return q.peak


def _lookahead_tree(C, U, V, q, r):
    n = len(U)
    blocks = []
    i = 0
    while i < n - 1:
        b = min(r, n - 1 - i)
        for j in range(2, b + 1):
            C.ccx(U[i + j - 1], V[i + j], U[i + j])
        T = {1: V[i + 1]}
        for j in range(2, b + 1):
            T[j] = q.take()
            C.ccx(T[j - 1], V[i + j], T[j])
        for j in range(b - 1, 1, -1):
            C.ccx(T[j - 1], V[i + j], T[j])
            q.give(T.pop(j))
        blocks.append((i, b, T))
        i += b
    G = [U[0]] + [U[i + b] for i, b, _ in blocks]
    P = [None] + [T[b] for _, b, T in blocks]
    made = []
    m = len(G)
    d = 1
    while d < m:
        for R in range(2 * d - 1, m, 2 * d):
            L = R - d
            # The product first: it is the next level's P, and both gates hold P[R].
            own = P[R]
            if P[L] is None:
                P[R] = None
            else:
                t = q.take()
                C.ccx(P[L], own, t)
                made.append((P[L], own, t))
                P[R] = t
            C.ccx(G[L], own, G[R])
        d *= 2
    while d > 1:
        d //= 2
        for R in range(3 * d - 1, m, 2 * d):
            C.ccx(G[R - d], P[R], G[R])
    for a, c, t in reversed(made):
        C.ccx(a, c, t)
        q.give(t)
    for i, b, T in blocks:
        for j in range(2, b):
            T[j] = q.take()
            C.ccx(T[j - 1], V[i + j], T[j])
        for j in range(1, b):
            C.ccx(U[i], T[j], U[i + j])
        for j in range(b, 1, -1):
            C.ccx(T[j - 1], V[i + j], T[j])
            q.give(T.pop(j))


# Block sizes lookahead_sweep tries, in order. Blocks of 1 make the tree a plain Brent-Kung prefix
# over the bits, with each p on its V wire: fewer CCX and fewer ticks than blocks of 3, for about
# len(U) temps instead of 2 len(U) / 3.
SWEEP_BLOCKS = (1, 2, 3, 4, 5, 6, 8)


def lookahead_sweep(C, U, V, temps):
    """lookahead_chain's map (f_j = x_j ^ p_j f_(j-1) in place on U from U[0]): lookahead_tree with
    the first block size in SWEEP_BLOCKS whose tree_temps fits len(temps), else lookahead_chain with
    lookahead_sizes. Larger blocks need fewer temps and are slower."""
    n = len(U)
    for r in SWEEP_BLOCKS:
        if tree_temps(n, r) <= len(temps):
            lookahead_tree(C, U, V, temps, r)
            return
    lookahead_chain(C, U, V, temps, lookahead_sizes(n - 1))


def add_lookahead(C, a, b, cout, temps, tree=False, cin=None):
    """b += a mod 2^w and cout ^= the carry out (none if cout is None), as add, with both sweeps in
    lookahead form. A wire holding 1 as cin adds 1, as in add; the sweeps restore it.
    The up-sweep leaves add's state after its last _maj (b_i = p_i, the carry wire c_0 ^ a_0,
    a_i = c_(i+1) ^ a_(i+1), a_(w-1) = c_w). The carries of a + ~s are those of a + b, so that
    state for (a, ~s) differs only on b, by b_i = ~(p_i ^ c_i ^ a_i): one CX from the wire below
    and an X a bit. The same sweep run backwards then leaves (a, ~s). Temps: clean wires, any
    number, left clean. With tree, the up-sweep's adjacent XORs are adjacent_xor's network and its
    carries are lookahead_sweep's, the same maps in log depth when the temps fit a tree."""
    w = len(a)
    assert len(b) == w and w >= 1 and temps
    own = cin is None
    if own:
        (cin,) = C.alloc(1)
    up = _Gates()
    for i in range(w):
        up.cx(a[i], b[i])
    up.cx(a[0], cin)
    if tree:
        adjacent_xor(up, a)
        lookahead_sweep(up, [cin] + list(a), [None] + list(b), temps)
    else:
        for i in range(w - 1):
            up.cx(a[i + 1], a[i])
        lookahead_chain(up, [cin] + list(a), [None] + list(b), temps, lookahead_sizes(w))
    _emit(C, up.ops)
    if cout is not None:
        C.cx(a[-1], cout)
    C.cx(cin, b[0])
    for i in range(1, w):
        C.cx(a[i - 1], b[i])
    for x in b:
        C.x(x)
    _emit(C, reversed(up.ops))
    for x in b:
        C.x(x)
    if own:
        C.free([cin])


def _up_sweep(a, b, cin, temps):
    """add_lookahead's tree up-sweep on (a, b) as recorded gates: b_i = p_i, cin = c_0 ^ a_0,
    a_i = c_(i+1) ^ a_(i+1), a_(w-1) = c_w."""
    up = _Gates()
    for i in range(len(a)):
        up.cx(a[i], b[i])
    up.cx(a[0], cin)
    adjacent_xor(up, a)
    lookahead_sweep(up, [cin] + list(a), [None] + list(b), temps)
    return up


def cadd_lookahead(C, a, b, cout, ctrl, temps, copies):
    """b += ctrl*a mod 2^w and cout ^= ctrl*carry out: add_lookahead's tree with its middle (the
    sum's CXs and both X layers) under ctrl, through ctrl's copies on clean wires. With ctrl = 0 the
    up-sweep meets its own reverse and leaves (a, b)."""
    w = len(a)
    assert len(b) == w and w >= 1 and temps and ctrl not in (*a, *b, cout)
    (cin,) = C.alloc(1)
    up = _up_sweep(a, b, cin, temps)
    _emit(C, up.ops)
    fan = _Gates()
    _doubling(fan, ctrl, copies)
    _emit(C, fan.ops)
    src = [ctrl] + list(copies)
    if cout is not None:
        C.ccx(src[-1], a[-1], cout)
    C.ccx(src[0], cin, b[0])
    for i in range(1, w):
        C.ccx(src[i % len(src)], a[i - 1], b[i])
    for i, x in enumerate(b):
        C.cx(src[i % len(src)], x)
    _emit(C, reversed(up.ops))
    for i, x in enumerate(b):
        C.cx(src[i % len(src)], x)
    _emit(C, reversed(fan.ops))
    C.free([cin])


def gt_lookahead(C, a, b, t, ctrl, temps):
    """gt_into's map, t ^= ctrl*[a > b] with a and b restored, reading the carry out of a + ~b from
    add_lookahead's tree up-sweep."""
    assert len(a) == len(b) >= 1 and temps and ctrl not in (*a, *b, t)
    for x in b:
        C.x(x)
    (cin,) = C.alloc(1)
    up = _up_sweep(a, b, cin, temps)
    _emit(C, up.ops)
    _cx(C, ctrl, a[-1], t)
    _emit(C, reversed(up.ops))
    C.free([cin])
    for x in b:
        C.x(x)


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


def _doubling(C, src, copies):
    """CX src into each clean copy, each copy a source once it holds src: about log2 depth."""
    have, todo = [src], list(copies)
    while todo:
        new = []
        for s in have:
            if not todo:
                break
            d = todo.pop(0)
            C.cx(s, d)
            new.append(d)
        have += new


def fan_into(C, src, targets, copies):
    """targets[i] ^= src for every i, through src's copies on clean wires, which end clean: about
    2 log2 len(copies) + len(targets) / (len(copies) + 1) CX layers instead of len(targets)."""
    g = _Gates()
    _doubling(g, src, copies)
    _emit(C, g.ops)
    srcs = [src] + list(copies)
    for i, t in enumerate(targets):
        C.cx(srcs[i % len(srcs)], t)
    _emit(C, reversed(g.ops))


def _carry_tree(g, G, P, q):
    """G[-1] ^= the carry out of c_(j+1) = G[j] ^ P[j] c_j from c_0 = 0, where G[j] and P[j] are never
    both 1, by a balanced tree of segment joins G_R ^= P_R G_L, P_R := P_L P_R, about log2 len(G) CCX
    layers. The G wires keep their segments' generates; products go on temps from the queue q. P[0] is
    not read."""
    segs = [(G[0], None)] + list(zip(G[1:], P[1:]))
    while len(segs) > 1:
        nxt = []
        for k in range(0, len(segs) - 1, 2):
            (gl, pl), (gr, pr) = segs[k], segs[k + 1]
            g.ccx(gl, pr, gr)
            if pl is None:
                nxt.append((gr, None))
            else:
                t = q.take()
                g.ccx(pl, pr, t)
                nxt.append((gr, t))
        if len(segs) % 2:
            nxt.append(segs[-1])
        segs = nxt


def gt_const_tree(C, value, b, t, ctrl=None):
    """gt_const_into's map, t ^= ctrl*[value > b] with b restored, with the carry out of ~b + value
    joined by _carry_tree instead of a chain. Leaf i: g = ~b_i and p = b_i where value's bit i is 1,
    g = 0 and p = ~b_i where it is 0. Bits below value's lowest 1 carry nothing."""
    assert ctrl is None or ctrl not in (*b, t)
    assert t not in b
    assert value < 1 << len(b)
    if value == 0:
        return
    s = (value & -value).bit_length() - 1
    idx = range(s, len(b))
    ones = [i for i in idx if value >> i & 1]
    flips = [b[i] for i in idx if not value >> i & 1]
    G = C.alloc(len(idx))
    where = dict(zip(idx, G))
    temps = C.alloc(len(idx))
    for x in flips:
        C.x(x)
    for i in ones:
        C.cx(b[i], where[i])
        C.x(where[i])
    g = _Gates()
    _carry_tree(g, G, [b[i] for i in idx], _TempQueue(temps))
    _emit(C, g.ops)
    _cx(C, ctrl, G[-1], t)
    _emit(C, reversed(g.ops))
    for i in ones:
        C.x(where[i])
        C.cx(b[i], where[i])
    for x in flips:
        C.x(x)
    C.free(temps)
    C.free(G)


def _add_selected_tree(C, b, selectors, masks, max_temps=None):
    """_add_selected_const's map, b += the mask its one set selector picks mod 2^len(b), with the
    addend copied onto clean wires by doubling fans and added by add_lookahead's tree, on at most
    max_temps temps (lookahead_sweep then takes the block size that fits)."""
    n = len(b)
    masks = [mask % (1 << n) for mask in masks]
    start = min(((mask & -mask).bit_length() - 1 for mask in masks if mask), default=n)
    low = b[start:]
    if not low:
        return
    A = C.alloc(len(low))
    classes = {}
    for i, w in enumerate(A):
        key = tuple(s for s, mask in zip(selectors, masks) if mask >> (start + i) & 1)
        if key:
            classes.setdefault(key, []).append(w)
    g = _Gates()
    for key, wires in classes.items():
        if len(key) == 1:
            _doubling(g, key[0], wires)
        else:
            for s in key:
                g.cx(s, wires[0])
            _doubling(g, wires[0], wires[1:])
    _emit(C, g.ops)
    want = max(1, tree_temps(len(low) + 1, 1))
    temps = C.alloc(want if max_temps is None else min(want, max_temps))
    add_lookahead(C, A, low, None, temps, tree=True)
    C.free(temps)
    _emit(C, reversed(g.ops))
    C.free(A)


# Run the top-level constant adds' carry chains as log-depth trees on clean wires.
CONST_TREE = True
# Copies of the control for creflect's complement of b.
REFLECT_COPIES = 31


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
        gt, add_sel = (gt_const_tree, _add_selected_tree) if CONST_TREE else (gt_const_into, _add_selected_const)
        (k,) = C.alloc(1)
        gt(C, self.q - value, b, k, ctrl)
        C.cx(ctrl, k)
        C.cx(k, ctrl)
        add_sel(C, b, [ctrl, k], [value, value + (1 << self.n) - self.q])
        C.cx(k, ctrl)
        gt(C, value, b, k, ctrl)
        C.free([k])

    def creflect(self, C, ctrl, a, b):
        """b = (a - b) mod q when ctrl, valid on all canonical b."""
        a %= self.q
        gt, add_sel = (gt_const_tree, _add_selected_tree) if CONST_TREE else (gt_const_into, _add_selected_const)
        (k,) = C.alloc(1)
        gt(C, a + 1, b, k, ctrl)
        C.cx(ctrl, k)
        if CONST_TREE:
            copies = C.alloc(REFLECT_COPIES)
            fan_into(C, ctrl, b, copies)
            C.free(copies)
        else:
            for wire in b:
                C.cx(ctrl, wire)
        C.cx(k, ctrl)
        add_sel(C, b, [ctrl, k], [a + 1, (a + 1 + self.q) % (1 << self.n)])
        C.cx(k, ctrl)
        C.cx(ctrl, k)
        gt(C, a + 1, b, k, ctrl)
        C.free([k])

    def cneg(self, C, ctrl, b):
        """b = q - b when ctrl, valid for b in [1, q). A bijection on all of
        {0,1}^n either way, since it is ~(b + 2^n - q - 1) mod 2^n."""
        add_const(C, b, (1 << self.n) - self.q - 1, ctrl=ctrl)
        for w in b:
            C.cx(ctrl, w)
