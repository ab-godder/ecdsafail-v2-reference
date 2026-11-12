"""Y^2 mod q on a fold of 2^n = f, for q = 2^n - f with f*f <= q.

Write ~Y = 2^n - 1 - Y. Each pair i < j gives y_i*~y_j + y_j*~y_i = y_i ^ y_j, and the
diagonal gives 0, so Y*~Y = sum over i < j of (y_i ^ y_j)*2^(i + j), a sum of rows with no
control on any of them. Q = Y + Y*~Y = Y*(2^n - Y), so Y^2 = 2^n*(Y - hi) - lo for
Q = lo + 2^n*hi, and Y - hi = D is at least 0. The rest folds 2^n to f:

    Z = (2^n - lo) + f*D = lo1 + 2^n*E, E <= f
    V = lo1 + f*E, V = Y^2 + f mod q, 0 <= V - f, V < 2^n + f*f

If V >= 2^n, dropping that bit subtracts q and f together and leaves V - 2^n <= f*f - 1 < q.
Otherwise V - f < q, so subtracting f is enough. So no comparison with q is needed.
"""
from .arith import _carries, _uncarry, add, add_const, gt_into, sub
from .builder import Circuit
from .two_bit import Prefix

# The longest block of the square's Prefixes. Shorter blocks take ticks off each write for a few
# more top wires, which fit under the multiplier's peak.
PREFIX_CAP = 37


def _lend(C, wires, fn):
    """Run fn with wires, which hold 0, handed out by C's allocator before any other. fn returns
    every ancilla it takes clean, so they hold 0 again after, and they leave the free list."""
    C.free(wires)
    fn()
    held = set(wires)
    assert held <= set(C.free_stack), "a lent wire was kept"
    C.free_stack = [w for w in C.free_stack if w not in held]


def _keep_add(C, a, b, cout, K):
    """b += a, cout ^= carry out, each inner carry left on a junk wire popped from K."""
    c = None
    for k, (x, y) in enumerate(zip(a, b)):
        t = cout if k == len(a) - 1 else K.pop(0)
        if c is None:
            C.ccx(x, y, t)
            C.cx(x, y)
        else:
            # MAJ(x, y, c) = c ^ (x ^ c)(y ^ c); the sum x ^ y ^ c lands on y.
            C.cx(c, x)
            C.cx(c, y)
            C.ccx(x, y, t)
            C.cx(c, t)
            C.cx(x, y)
            C.cx(c, y)
            C.cx(c, x)
        c = t


def _kept_inc(C, bits, t, K, cout=None):
    """bits += t mod 2^len(bits), cout ^= the carry out. The carry into bits[i], t AND bits[:i]
    all 1, is left on K[i - 1]; K's len(bits) - 1 wires hold 0 on entry. One CCX a bit."""
    carries = [t] + list(K)
    assert len(carries) == max(1, len(bits))
    for i in range(len(bits) - 1):
        C.ccx(carries[i], bits[i], carries[i + 1])
    if cout is not None:
        C.ccx(carries[-1], bits[-1], cout)
    for x, y in zip(carries, bits):
        C.cx(x, y)


def _rows(C, Y, Q, copy=True, K=None, cin=None):
    """Q = Y*(2^w - Y) on 2w clean wires, w = len(Y), or Y*(2^w - 1 - Y) without the copy.
    After row i, Q < 2^(w + i + 1), so each carry lands on a wire no earlier row reached; the
    uncopied sums are smaller. Y's wires walk through their Gray code. K, if given, is a list
    of junk wires, and a row whose inner carries fit there keeps them, at one CCX a bit. cin, if
    given, is a clean wire every add takes as its carry in, so two rows on their own cins share
    no wire."""
    n = len(Y)
    if copy:
        for y, q in zip(Y, Q):
            C.cx(y, q)
    for i in range(n - 1):
        for w in Y[i + 1:]:
            C.cx(Y[i], w)
        if K is not None and len(K) >= n - 2 - i:
            _keep_add(C, Y[i + 1:], Q[2 * i + 1:n + i], Q[n + i], K)
        else:
            add(C, Y[i + 1:], Q[2 * i + 1:n + i], cout=Q[n + i], cin=cin)
    for i in range(1, n):
        C.cx(Y[i - 1], Y[i])


def _change(C, a, b, sign=1, cin=None, signed=False):
    """b += sign*a (+ cin) mod 2^len(b); a zero-extended, or sign-extended when signed."""
    assert len(a) <= len(b)
    pad = C.alloc(len(b) - len(a))
    if signed:
        for p in pad:
            C.cx(a[-1], p)
    if sign < 0:
        for w in b:
            C.x(w)
    add(C, list(a) + pad, b, cin=cin)
    if sign < 0:
        for w in b:
            C.x(w)
    if signed:
        for p in pad:
            C.cx(a[-1], p)
    C.free(pad)


def karatsuba_junk(w, thr):
    """Junk wires _karatsuba takes for width w with leaves of width at most thr."""
    if w <= thr:
        return 0
    h = (w + 1) // 2
    ws = h + 1
    return karatsuba_junk(h, thr) + karatsuba_junk(w - h, thr) + 3 * ws + karatsuba_junk(ws, thr)


def karatsuba_carries(w, thr):
    """Inner row carries _karatsuba's leaves can keep for width w: a leaf of width v has rows
    of v - 2, ..., 1, 0 inner carries."""
    if w <= thr:
        return (w - 1) * (w - 2) // 2 if w > 1 else 0
    h = (w + 1) // 2
    return karatsuba_carries(h, thr) + karatsuba_carries(w - h, thr) + karatsuba_carries(h + 1, thr)


def _karatsuba(C, Y, Q, J, thr, kind="R", K=None, Z=None, cin=None):
    """Q = R(Y) = Y*(2^w - Y) (kind R) or T(Y) = Y*(2^w - 1 - Y) (kind T) on 2w clean wires.

    Y = A + 2^h*B with h = ceil(w/2), m = w - h, and S = A + B on ws = h + 1 bits. Then
        R(Y) = R(A) + 2^(2h)*T(B) + 2^h*V,   V = T(S) - R(A) - T(B) - [h != m]*2^m*S
        T(Y) = T(A) + 2^(2h)*T(B) + 2^h*V,   V = R(~S) - T(A) - T(B) - [h != m]*2^m*S
                                                  - (2^ws - 1)
    where ~S = 2^ws - 1 - S, so R(~S) = T(S) - S + 2^ws - 1. V is A*~B + B*~A, plus B for R,
    so 0 <= V < 2^(w + 1), and V is formed mod 2^(w + 1). R(A) and T(B) go straight onto Q's
    halves; S and the S child's output take wires from J, a list of clean junk wires, and
    keep values fixed by Y. Leaves of width at most thr are the XOR rows. Z, if given, gets the
    S child's output wires above V that its value leaves 0: R(x) <= 2^(2ws - 2) clears the top
    bit and T(x) < 2^(2ws - 2) the top two, and nothing above V is written after the child.
    cin, if given, is a clean carry in for the leaves' rows."""
    w = len(Y)
    if w <= thr:
        _rows(C, Y, Q, kind == "R", K, cin)
        return
    h = (w + 1) // 2
    m = w - h
    A, B = Y[:h], Y[h:]
    _karatsuba(C, A, Q[:2 * h], J, thr, kind, K, Z, cin)
    _karatsuba(C, B, Q[2 * h:], J, thr, "T", K, Z, cin)
    ws = h + 1
    S = [J.pop(0) for _ in range(ws)]
    for x, s in zip(A, S):
        C.cx(x, s)
    pad = C.alloc(h - m)
    add(C, B + pad, S[:h], cout=S[h])
    C.free(pad)
    W = [J.pop(0) for _ in range(2 * ws)]
    if kind == "T":
        for x in S:
            C.x(x)
    child = "R" if kind == "T" else "T"
    _karatsuba(C, S, W, J, thr, child, K, Z, cin)
    if Z is not None:
        Z.extend(W[max(w + 1, 2 * ws - (1 if child == "R" else 2)):])
    if kind == "T":
        for x in S:
            C.x(x)
    V = W[:w + 1]
    _change(C, Q[:2 * h], V, -1)
    _change(C, Q[2 * h:], V, -1)
    if h != m:
        _change(C, S, V[m:], -1)
    if kind == "R":
        _change(C, V, Q[h:])
        return
    # V holds the cross term plus 2^ws - 1. Take 2^ws off its top bits, which leaves at
    # least -1 on w + 1 bits, then add it sign-extended with a carry in of 1.
    add_const(C, V[ws:], (1 << (w + 1 - ws)) - 1)
    (one,) = C.alloc(1)
    C.x(one)
    _change(C, V, Q[h:], cin=one, signed=True)
    C.x(one)
    C.free([one])



def signed_digits(v):
    """Non-adjacent form of v, as (position, sign)."""
    out, s = [], 0
    while v:
        if v & 1:
            d = 2 - (v & 3)
            out.append((s, d))
            v -= d
        v >>= 1
        s += 1
    return out


class SquareFold:
    # Junk wires that keep leaf row carries, capped at what the leaves have. Each saves 2 CCX,
    # one in the square and one in its reverse, and adds a wire at the square's peak, which is
    # the circuit's: 650 puts secp256k1 at 4,096 wires.
    KEEP = 650
    def __init__(self, mod):
        self.mod, self.n = mod, mod.n
        self.f = f = (1 << mod.n) - mod.q
        self.b, self.m = f.bit_length(), (f * f).bit_length()
        # The low digit is +1 whatever f is, so its carry in can hold Z's + 1.
        self.digits = [(0, 1)] + signed_digits(f - 1)
        # Karatsuba leaves of at most n // 6 bits: 42 at n = 256, and 3 at the toy sizes so
        # the recursion runs there. KEEP junk wires hold row carries, one CCX saved per wire.
        self.thr = max(3, mod.n // 6)
        # The two increments above the fold's low bits keep their carries, n - m - 1 and
        # n - b - 1 wires, 1 CCX a bit in place of a Prefix's 3.
        self.incs = (max(0, mod.n - self.m - 1), max(0, mod.n - self.b - 1))
        self.fixed_junk = 2 * mod.n + self.b + self.m + 3 + sum(self.incs)
        self.keep = min(self.KEEP, karatsuba_carries(mod.n, self.thr))
        self.junk = self.fixed_junk + karatsuba_junk(mod.n, self.thr) + self.keep

    @staticmethod
    def applies(q):
        f = (1 << q.bit_length()) - q
        return f * f <= q

    def modadd(self, C, a, b, ctrl, Z=()):
        _lend(C, list(Z), lambda: self._modadd(C, a, b, ctrl))

    def _modadd(self, C, a, b, ctrl):
        """Mod.modadd with _reduce's carries above f's window held by a Prefix, so it runs
        beside the wires compute leaves live. The same flag, the same write, the same CCX.
        Z lists wires holding 0, handed to its ancillas first and returned holding 0."""
        f, w = self.f, self.b
        top, t = C.alloc(2)
        add(C, a, b, top, ctrl)
        cs = _carries(C, b[:w + 1], f)
        carry = Prefix(C, b[w:], cap=PREFIX_CAP)
        # s's top bit and the carry into it are never both 1, so their XOR is the flag.
        C.ccx(cs[w], carry.top, top)
        C.ccx(top, cs[w], t)
        carry.write(t)
        C.ccx(top, cs[w], t)
        for i in reversed(range(w)):
            _uncarry(C, b, f, cs, i + 1)
            if cs[i] is not None:
                C.ccx(top, cs[i], b[i])
            if f >> i & 1:
                C.cx(top, b[i])
        gt_into(C, a, b, top, ctrl)
        C.free([top, t])

    def _multiple(self, C, a, T, a_max, t_max=0, cin=None):
        """T += f*a mod 2^len(T), plus 1 if cin holds 1, for a <= a_max and T <= t_max on entry.

        Around each digit T holds t + d*a, plus 1 once cin is in, with d the digits so far.
        Those values at the corners of the box bound every input. A subtract adds onto ~T, so
        it bounds -v - 1 instead. Where every bound lies in [-2^i, 2^i) for a bit i at or above
        a's top, the digit's bits from i up all hold the sign. A sign that is always 0 never
        moves, so the digit adds below bit i only. Otherwise bit i takes the carry out and the
        bits above copy it. A digit onto a T that is 0, with no carry in, copies a."""
        inc, done, d = int(cin is not None), 0, 0
        for s, sign in self.digits:
            box = [(t, x) for t in (0, t_max) for x in (0, a_max)]
            before = [t + done + x * d for t, x in box]
            d += sign << s
            after = [t + inc + x * d for t, x in box]
            done = inc
            target = T[s:]
            if not target:
                continue
            if not any(before) and cin is None and sign > 0:
                for x, w in zip(a, target):
                    C.cx(x, w)
                continue
            ends = [v >> s for v in before + after]
            lo, hi = (min(ends), max(ends)) if sign > 0 else (-max(ends) - 1, -min(ends) - 1)
            top = next((i for i in range(len(a), len(target)) if -(1 << i) <= lo and hi < 1 << i),
                       len(target))
            carry = lo < 0 and top < len(target)
            above = target[top + 1:] if carry else []
            flips = [] if sign > 0 else target if carry else target[:top]
            for w in flips:
                C.x(w)
            for w in above:
                C.cx(target[top], w)
            pad = C.alloc(max(0, top - len(a)))
            add(C, a[:top] + pad, target[:top], cout=target[top] if carry else None, cin=cin)
            C.free(pad)
            for w in above:
                C.cx(target[top], w)
            for w in flips:
                C.x(w)
            cin = None

    def compute(self, C, Y, J):
        """J clean on entry. On exit J[:n] holds Y^2 mod q, canonical, for any n-bit Y.
        The other wires of J hold values fixed by Y, for the reverse of this block."""
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:2 * n + b + m + 3]
        inc = J[2 * n + b + m + 3:self.fixed_junk]
        inc_t, inc_k = inc[:self.incs[0]], inc[self.incs[0]:]
        junk = list(J[self.fixed_junk:self.junk - self.keep])
        keep = list(J[self.junk - self.keep:]) if self.keep else None
        # E, T, t, c, k and inc are first written after the product, so its ancillas can borrow
        # them: the root's add onto its upper half pads V with w - h - 1 zero wires.
        zeros = []
        _lend(C, E + T + [t, c, k] + inc,
              lambda: _karatsuba(C, Y, lo + hi, junk, self.thr, "R", keep, zeros))
        assert not junk and not keep
        # J's places that hold 0 from here to the exit: the fold and the modadd beside this
        # block borrow them.
        self.zero_slots = [J.index(z) for z in zeros]
        _lend(C, zeros, lambda: self._fold(C, Y, lo, hi, E, T, t, c, k, inc_t, inc_k))

    def _fold(self, C, Y, lo, hi, E, T, t, c, k, inc_t, inc_k):
        """lo + 2^n*hi = Y*(2^n - Y) on entry, the other registers 0. On exit lo = Y^2 mod q."""
        n, f, b, m = self.n, self.f, self.b, self.m
        # hi = Y - hi = ~(hi + ~Y), and lo = ~lo.
        for w in Y:
            C.x(w)
        add(C, Y, hi)
        for w in Y:
            C.x(w)
        for w in hi + lo:
            C.x(w)
        (one,) = C.alloc(1)
        C.x(one)
        # hi is any n-bit D, and lo + E holds ~lo, below 2^n.
        self._multiple(C, hi, lo + E, (1 << n) - 1, (1 << n) - 1, one)
        C.x(one)
        C.free([one])
        self._multiple(C, E, T, f)
        high = lo[m:]
        add(C, T, lo[:m], cout=t)
        if high:
            _kept_inc(C, high, t, inc_t, cout=c)
        else:
            C.cx(t, c)
        # With c = 0, lo >= f. Take f off the low bits, and the borrow off the rest.
        C.x(c)
        add_const(C, lo[:b], (1 << b) - f, c, k)
        C.cx(c, k)
        high = lo[b:]
        for w in high:
            C.x(w)
        _kept_inc(C, high, k, inc_k)
        for w in high:
            C.x(w)
        C.x(c)


def _replayed(C, build, use):
    """Run build(), then use(result), then build's gates again in reverse order, which undoes
    them: every gate is self-inverse. use must leave every wire build touched as it found it.
    build must emit primitives only. Returns build's result."""
    start, calls = len(C.prims), len(C.calls) if C.calls is not None else 0
    r = build()
    end = len(C.prims)
    assert C.calls is None or len(C.calls) == calls, "build emitted a call"
    use(r)
    ops = C.prims[start:end]
    for j in reversed(range(0, len(ops), 4)):
        C.prims.extend(ops[j:j + 4])
    return r


def _and_top(C, bits):
    """Prefix's top, AND(bits), built on about log2(len(bits)) fresh wires and bits themselves.
    Returns (top, wires to free once the build is replayed backwards)."""
    p = Prefix(C, bits, cap=PREFIX_CAP)
    own = [] if len(bits) < 2 else [p.first] + [blk[1] for blk in p.blocks]
    return p.top, own


class HalfSquares:
    """X -= gen*(Y^2 + K) mod q by one Karatsuba level taken mod q, for even n and small f.

    Write h = n/2, beta = 2^h, Y = L + beta*H and S = L + H, so beta^2 = 2^n = f mod q and
        Y^2 = L^2 + beta*(S^2 - L^2) - (beta - f)*H^2  (mod q).
    Each half square Z^2 is the rows of SquareFold followed by its conversion, which leaves
    Z' = Z^2 + 2^k - 1 on 2k wires, k = len(Z). Two groups each build one canonical residue
    from Y, which a gen-controlled modadd takes into X, and are then unbuilt:
        AC = (L' + beta*(S' - L')) mod q,   B = (beta*B0 + f*(2^n - 1 - P)) mod q,
    with H' = B0 + beta*B1 and P = H' - B1. AC is subtracted and B added, which leaves
    X - gen*(Y^2 - K) for the constant K below; the caller adds gen*K. Every bound the folds
    need is checked on the integers in __init__ and asserted there."""

    # Leaf widths for a Karatsuba inside each group's half squares, or None for plain rows. A
    # level's sum and middle product stay on the group's junk until its unbuild.
    THR_AC = None
    THR_B = None

    def __init__(self, mod):
        self.mod, self.q, self.n = mod, mod.q, mod.n
        b = self.bounds(self.q)
        assert b is not None, "HalfSquares does not apply"
        self.__dict__.update(b)
        n, h, f = self.n, self.h, self.f
        self.helper = SquareFold(mod)
        beta = 1 << h
        # AC carries beta - 1 + 2^n beyond its share of Y^2 and B carries (beta - f)*(beta - 1)
        # + f*(2^n - 1), so X - gen*AC + gen*B = X - gen*Y^2 + gen*K.
        K = (beta - f) * (beta - 1) + f * ((1 << n) - 1) - (beta - 1 + (1 << n))
        self.K = K % self.q

    @staticmethod
    def bounds(q):
        """The widths the folds use, or None where a bound they need fails."""
        n = q.bit_length()
        f = (1 << n) - q
        if n % 2 or n < 8 or f <= 0:
            return None
        h, fb = n // 2, f.bit_length()
        beta = 1 << h
        half_max = (beta - 1) ** 2 + beta - 1             # Z' for a k = h half square
        # G = S' - L' = H*(2L + H) + 2^h, largest at L = H = beta - 1.
        g1_max = (3 * (beta - 1) ** 2 + beta) >> h
        tb = (f * g1_max).bit_length()
        # AC before its last fold: L' + f*G1 + beta*G0 < 3 * 2^n, so its two carries are never
        # both 1 and e = c1 + 2*c2 is at most 2.
        v_ac = half_max + f * g1_max + beta * (beta - 1)
        # B: T = beta*B0 + f*~P on n + hib bits, then lo + f*hi carries at most once.
        t_b = beta * (beta - 1) + f * ((1 << n) - 1)
        hi_max = t_b >> n
        hib, e2b = hi_max.bit_length(), (f * hi_max).bit_length()
        ok = (v_ac < 3 << n and tb <= n and ((1 << n) - 1) + f * hi_max < 2 << n
              and e2b <= n and 3 * f < 1 << (fb + 2) and fb + 2 <= n and 2 * f < q
              and t_b < 1 << (n + hib))
        # The finish: R + f*e < 2^n + 2f, and when it carries, R + f < 3f < 2^(fb + 2) <= q.
        ok = ok and 3 * f < q
        return dict(h=h, f=f, fb=fb, g1_max=g1_max, tb=tb, hi_max=hi_max, hib=hib,
                    e2b=e2b) if ok else None

    @staticmethod
    def applies(q):
        return HalfSquares.bounds(q) is not None

    def _half(self, C, Z, W, cin, J=None, thr=None, lent=()):
        """W, 2k clean wires, becomes Z^2 + 2^k - 1. Rows, or a Karatsuba with leaves of at most
        thr bits whose junk comes from J, then hi = Z - hi and lo = ~lo. lent are wires that
        hold 0 and are written only later, which the product's ancillas may take first."""
        k = len(Z)
        if thr is None or k <= thr:
            _rows(C, Z, W, cin=cin)
        else:
            junk = [J.pop(0) for _ in range(karatsuba_junk(k, thr))]
            _lend(C, list(lent), lambda: _karatsuba(C, Z, W, junk, thr, "R", cin=cin))
            assert not junk
        lo, hi = W[:k], W[k:]
        for w in Z:
            C.x(w)
        add(C, Z, hi, cin=cin)
        for w in Z:
            C.x(w)
        for w in hi + lo:
            C.x(w)

    @staticmethod
    def _add_into(C, a, b, J, cout):
        """b += a for len(a) <= len(b), cout ^= the carry out of b. The carry out of b's low
        len(a) bits stays on a junk wire from J; a Prefix takes it through b's high bits."""
        k = len(a)
        high = b[k:]
        if not high:
            add(C, a, b, cout=cout)
            return
        t = J.pop(0)
        add(C, a, b[:k], cout=t)
        carry = Prefix(C, high, cap=PREFIX_CAP)
        C.ccx(t, carry.top, cout)
        carry.write(t)

    def _canon(self, C, R, J):
        """R < 2^n becomes R mod q. k = [R >= q] = [R + f >= 2^n] stays on a junk wire."""
        f, fb = self.f, self.fb
        low, high = R[:fb], R[fb:]
        k = J.pop(0)

        cs = _carries(C, R[:fb + 1], f)
        _, own = _replayed(C, lambda: _and_top(C, high), lambda r: C.ccx(cs[fb], r[0], k))
        C.free(own)
        for i in reversed(range(fb + 1)):
            _uncarry(C, R[:fb + 1], f, cs, i)
        (t2,) = C.alloc(1)
        add_const(C, low, f, ctrl=k, cout=t2)
        # k = 1 means high was all ones and the low add carried out: R - q = R + f - 2^n.
        for w in high:
            C.cx(k, w)
        C.cx(k, t2)
        C.free([t2])

    def _finish(self, C, R, c1, c2, J):
        """R + 2^n*(c1 + 2*c2), c1 and c2 never both 1, becomes its canonical residue on R."""
        f, fb = self.f, self.fb
        F = C.alloc(fb + 2)

        def load():
            for i in range(fb):
                if f >> i & 1:
                    C.cx(c1, F[i])
                    if c2 is not None:
                        C.cx(c2, F[i + 1])

        load()
        c3 = J.pop(0)
        self._add_into(C, F, R, J, c3)
        load()
        C.free(F)
        # R + f*c3 < 3f fits fb + 2 bits when c3 = 1.
        add_const(C, R[:fb + 2], f, ctrl=c3)
        self._canon(C, R, J)

    def build_ac(self, C, Y, J):
        """Returns the wires of AC, canonical; every other wire it leaves set is from J."""
        n, h, f = self.n, self.h, self.f
        L, H = Y[:h], Y[h:]
        s = J.pop(0)
        add(C, L, H, cout=s)
        S = H + [s]
        RA = [J.pop(0) for _ in range(n)]
        RC = [J.pop(0) for _ in range(n + 2)]
        cins = C.alloc(2)
        # T is written after both halves; each half's product ancillas take half of it first.
        ja, jc = (karatsuba_junk(k, self.THR_AC) if self.THR_AC else 0 for k in (h, h + 1))
        T = J[ja + jc:ja + jc + self.tb]
        self._half(C, L, RA, cins[0], J, self.THR_AC, T[:self.tb // 2])
        self._half(C, S, RC, cins[1], J, self.THR_AC, T[self.tb // 2:])
        C.free(cins)
        # G = S' - L' on RC; it is at least 2^h.
        pad = C.alloc(2)
        sub(C, RA + pad, RC)
        C.free(pad)
        G0, G1 = RC[:h], RC[h:]
        T = [J.pop(0) for _ in range(self.tb)]
        self.helper._multiple(C, G1, T, self.g1_max)
        c1, c2 = J.pop(0), J.pop(0)
        self._add_into(C, T, RA, J, c1)
        (pg,) = C.alloc(1)
        add(C, G0 + [pg], RA[h:] + [c1], cout=c2)
        C.free([pg])
        self._finish(C, RA, c1, c2, J)
        return RA

    def build_b(self, C, Y, J):
        n, h, f = self.n, self.h, self.f
        H = Y[h:]
        RB = [J.pop(0) for _ in range(n)]
        (cin,) = C.alloc(1)
        jb = karatsuba_junk(n - h, self.THR_B) if self.THR_B else 0
        self._half(C, H, RB, cin, J, self.THR_B, J[jb:jb + n + self.hib])
        C.free([cin])
        T = [J.pop(0) for _ in range(n + self.hib)]
        for x, w in zip(RB[:h], T[h:]):
            C.cx(x, w)
        # RB becomes P = B0 - B1 + beta*(B1 - b), b = [B1 > B0], and then ~P.
        b = J.pop(0)
        sub(C, RB[h:], RB[:h], borrow=b)
        for w in RB[h:]:
            C.x(w)
        Prefix(C, RB[h:], cap=PREFIX_CAP).write(b)
        for w in RB[h:]:
            C.x(w)
        for w in RB:
            C.x(w)
        self.helper._multiple(C, RB, T, (1 << n) - 1, (1 << h) * ((1 << h) - 1))
        lo, hi = T[:n], T[n:]
        E = [J.pop(0) for _ in range(self.e2b)]
        self.helper._multiple(C, hi, E, self.hi_max)
        c1 = J.pop(0)
        self._add_into(C, E, lo, J, c1)
        self._finish(C, lo, c1, None, J)
        return lo

    def _recorded(self, C, build):
        """The build as a block on (Y, J), with the junk count and where its residue sits in J."""
        probe = Circuit()
        Y = probe.alloc(self.n)
        kj = sum(karatsuba_junk(k, t) for k, t in ((self.h, self.THR_AC), (self.h + 1, self.THR_AC),
                                                    (self.n - self.h, self.THR_B)) if t)
        J = probe.alloc(4 * self.n + 64 + kj)
        size = len(J)
        R = build(probe, Y, J)
        count = size - len(J)
        n = self.n
        holder = {}

        def fn(y, j):
            rest = list(j)
            out = build(C, y, rest)
            assert not rest
            holder["at"] = [list(j).index(w) for w in out]

        block = C.record_shape(fn, [n, count])
        return block, count, holder["at"]

    def mulsub(self, C, gen, Y, X):
        """X = X - gen*(Y^2 - K) mod q: the caller takes gen*K off."""
        n = self.n
        modadd = C.record_shape(lambda k, a, b: self.helper._modadd(C, a, b, k[0]), [1, n, n])
        for build, inverse in ((self.build_ac, True), (self.build_b, False)):
            block, count, at = self._recorded(C, build)
            J = C.alloc(count)
            C.apply(block, [Y, J])
            C.apply(modadd, [[gen], [J[i] for i in at], X], inverse=inverse)
            C.apply(block, [Y, J], inverse=True)
            C.free(J)
