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


def _maj_into(C, c, b, a, t):
    """t ^= MAJ(c, b, a), with c, b and a restored."""
    C.cx(a, b)
    C.cx(a, c)
    C.ccx(c, b, t)
    C.cx(a, t)
    C.cx(a, c)
    C.cx(a, b)


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
    """b += ctrl*a mod 2^w (CDKM), plus cin when ctrl is None. cout ^= ctrl*carry
    out. cin is a wire restored on return."""
    w = len(a)
    assert len(b) == w
    assert ctrl is None or ctrl not in (*a, *b, cout)
    assert ctrl is None or cin is None
    if w == 0:
        return
    own_cin = cin is None
    if own_cin:
        (cin,) = C.alloc(1)
    carries = [cin] + a[:-1]
    if ctrl is None:
        # The top bit's carry goes straight to cout, so a[-1] never has to hold it.
        for i in range(w - 1):
            _maj(C, carries[i], b[i], a[i])
        if cout is not None:
            _maj_into(C, carries[-1], b[-1], a[-1], cout)
        C.cx(a[-1], b[-1])
        C.cx(carries[-1], b[-1])
        for i in reversed(range(w - 1)):
            _uma(C, carries[i], b[i], a[i])
    else:
        for i in range(w):
            _maj(C, carries[i], b[i], a[i])
        if cout is not None:
            _cx(C, ctrl, a[-1], cout)
        for i in reversed(range(w)):
            _cuma(C, ctrl, carries[i], b[i], a[i])
    if own_cin:
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
    top = w if ctrl is not None else w - 1
    for i in range(top):
        _maj(C, carries[i], b[i], a[i])
    if ctrl is None:
        _maj_into(C, carries[-1], b[-1], a[-1], t)
    else:
        C.ccx(ctrl, a[-1], t)
    for i in reversed(range(top)):
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


def _carries(C, reg, value, ctrl=None, cin=None):
    """cs[i] is the carry into reg[i] of reg + ctrl*value + cin. Each entry is
    None when the carry is always 0, reg[i - 1] itself when it equals that bit,
    or a fresh wire. Release the fresh wires with _uncarry."""
    cs = [cin]
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


def add_const(C, reg, value, ctrl=None, cout=None, cin=None):
    """reg += ctrl*value mod 2^w, plus cin when ctrl is None, cout ^= ctrl*carry
    out. The carry chain reads the constant's bits, so 0 bits below the first 1
    cost nothing. cin is a wire restored on return."""
    w = len(reg)
    value %= 1 << w
    assert ctrl is None or ctrl not in (*reg, cout)
    assert ctrl is None or cin is None
    cs = _carries(C, reg, value, ctrl, cin)
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


SQUARE_DEPTH = 1


def _add_short(C, a, b, sub_=False):
    """b += a or b -= a, for a that fits in fewer wires than b. CDKM keeps carries
    on a's wires, so the shortfall takes clean wires for just this add."""
    zeros = C.alloc(max(0, len(b) - len(a)))
    (sub if sub_ else add)(C, (a + zeros)[:len(b)], b)
    C.free(zeros)


def _square_wires(n, depth):
    """Wires _square takes from its pool for an n-bit square."""
    if depth == 0 or n < 4:
        return 2 * n
    m = n // 2
    t = n - m + 1
    return _square_wires(m, depth - 1) + _square_wires(n - m, depth - 1) + _square_wires(t, depth - 1)


def _square(C, y, take, depth):
    """Wires holding y^2, 2n of them, from clean wires handed out by take. y is
    restored. A row adds only the y[j] above its y[i], and the Karatsuba split
    uses y^2 = a^2 + 4^m*b^2 + 2^m*((a + b)^2 - a^2 - b^2)."""
    n = len(y)
    if depth == 0 or n < 4:
        p = take(2 * n)
        C.cx(y[0], p[0])
        for i in range(n - 1):
            add(C, y[i + 1:], p[2 * i + 2:n + i + 1], p[n + i + 1], ctrl=y[i])
        zeros = C.alloc(n - 1)
        add(C, [w for pair in zip(y[1:], zeros) for w in pair], p[2:])
        C.free(zeros)
        return p
    m = n // 2
    a, b = y[:m], y[m:]
    low, high = _square(C, a, take, depth - 1), _square(C, b, take, depth - 1)
    # a + b is cheap to rebuild, so its wires go back before the folds need wires.
    t = C.alloc(n - m + 1)
    for w, u in zip(b, t):
        C.cx(w, u)
    _add_short(C, a, t)
    mid = _square(C, t, take, depth - 1)
    _add_short(C, a, t, sub_=True)
    for w, u in zip(b, t):
        C.cx(w, u)
    C.free(t)
    for part in (low, high):
        _add_short(C, part, mid, sub_=True)
    r = low + high
    _add_short(C, mid, r[m:])
    return r


class Mod:
    """Arithmetic modulo an odd q < 2^n on canonical values in [0, q)."""

    def __init__(self, q):
        self.q = q
        self.n = q.bit_length()

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

    def square_plan(self):
        """(folds, flags) for reducing y^2 by folding 2^n to f = 2^n - q. Each fold
        is (bound on the high part, wires it adds above the low n, wires for its
        product with f when that is built apart). flags counts the conditional
        subtractions of q << k that finish the reduction."""
        n, f = self.n, (1 << self.n) - self.q
        bound, folds = ((1 << n) - 1) ** 2, []
        while True:
            folded = (1 << n) - 1 + (bound >> n) * f
            if not bound >> n or folded >= bound:
                break
            h = bound >> n
            # Only a high part narrower than n pays for summing its product first:
            # then just one add ripples the whole width of the low bits.
            folds.append((h, folded.bit_length() - n, (h * f).bit_length() if h.bit_length() < n else 0))
            bound = folded
        return folds, (bound // self.q).bit_length()

    def square_pool_size(self):
        folds, flags = self.square_plan()
        return _square_wires(self.n, SQUARE_DEPTH) + sum(w + p for _, w, p in folds) + flags

    def square_into(self, C, y, pool):
        """Leaves y^2 mod q in pool[:n] for any n-bit y, with pool clean on entry
        and holding the garbage for the reverse pass. pool has square_pool_size()
        wires. Each fold adds the high part into the low n bits in place, and f
        only sets which shifts run, so a sparse f is cheap and a dense f is merely
        slower."""
        n, f = self.n, (1 << self.n) - self.q
        folds, flags = self.square_plan()
        wires = iter(pool)

        def take(k):
            return [next(wires) for _ in range(k)]

        cur = _square(C, y, take, SQUARE_DEPTH)
        low, top = cur[:n], cur[n:]
        for h, extra, apart in folds:
            low = low[:n] + take(extra)
            acc = take(apart) if apart else low
            done = 0
            for j in range(f.bit_length()):
                if f >> j & 1:
                    done |= 1 << j
                    span = (h * done + (0 if apart else (1 << n) - 1)).bit_length() - j
                    _add_short(C, top[:span], acc[j:j + span])
            if apart:
                _add_short(C, acc, low)
            top = low[n:]
        for k, t in zip(reversed(range(flags)), take(flags)):
            gt_const_into(C, self.q << k, low, t)
            C.x(t)
            add_const(C, low, -(self.q << k), ctrl=t)
        pool[:] = low[:n] + [w for w in pool if w not in set(low[:n])]

    def offset(self):
        """c = (2^n - q - 1)/2. Offset residues store x in [0, q) as x + c, which
        centres them in [0, 2^n): the flags of a modular double and a modular
        add become a top bit and a carry out, with no comparison against q."""
        return ((1 << self.n) - self.q - 1) >> 1

    def offset_add(self, C, a, b):
        """b = (b + a) mod q on offset residues. b + a + 1 carries out of n bits
        exactly when the residues sum to at least q, because 2c + 1 = 2^n - q."""
        (t, cin) = C.alloc(2)
        C.x(cin)
        add(C, a, b, t, cin=cin)
        C.x(cin)
        # t = 1 leaves the plain residue, which needs +c. t = 0 needs -(c + 1),
        # which is the complement of (~b + c + 1). Both run on one carry chain.
        C.x(t)
        for w in b:
            C.cx(t, w)
        add_const(C, b, self.offset(), cin=t)
        for w in b:
            C.cx(t, w)
        C.x(t)
        # The residue wrapped exactly when the result is below a; the offsets cancel.
        gt_into(C, a, b, t)
        C.free([t, cin])

    def offset_dbl(self, C, s):
        """s = 2s mod q on offset residues. The top bit is the flag [2x >= q], and
        rotating it in as the low bit leaves s + c under the flag and s - c
        otherwise, which is one chain on the upper bits and a complement frame."""
        c = self.offset()
        s[:] = s[-1:] + s[:-1]
        t, h = s[0], s[1:]
        C.x(t)
        for w in h:
            C.cx(t, w)
        add_const(C, h, (c + 1) >> 1)
        for w in h:
            C.cx(t, w)
        if not c & 1:
            C.x(t)

    def cmodadd_const(self, C, ctrl, value, b):
        """b = (b + ctrl*value) mod q."""
        value %= self.q
        if value == 0:
            return
        (top,) = C.alloc(1)
        add_const(C, b, value, ctrl, top)
        _reduce(C, b + [top], self._wrap())
        gt_const_into(C, value, b, top, ctrl)
        C.free([top])

    def cneg(self, C, ctrl, b):
        """b = q - b when ctrl, valid for b in [1, q). A bijection on all of
        {0,1}^n either way, since it is ~(b + 2^n - q - 1) mod 2^n."""
        add_const(C, b, (1 << self.n) - self.q - 1, ctrl=ctrl)
        for w in b:
            C.cx(ctrl, w)
