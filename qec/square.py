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
from .arith import _carries, _uncarry, add, add_const, gt_into
from .two_bit import Prefix


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


def _rows(C, Y, Q, copy=True, K=None):
    """Q = Y*(2^w - Y) on 2w clean wires, w = len(Y), or Y*(2^w - 1 - Y) without the copy.
    After row i, Q < 2^(w + i + 1), so each carry lands on a wire no earlier row reached; the
    uncopied sums are smaller. Y's wires walk through their Gray code. K, if given, is a list
    of junk wires, and a row whose inner carries fit there keeps them, at one CCX a bit."""
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
            add(C, Y[i + 1:], Q[2 * i + 1:n + i], cout=Q[n + i])
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


def _karatsuba(C, Y, Q, J, thr, kind="R", K=None):
    """Q = R(Y) = Y*(2^w - Y) (kind R) or T(Y) = Y*(2^w - 1 - Y) (kind T) on 2w clean wires.

    Y = A + 2^h*B with h = ceil(w/2), m = w - h, and S = A + B on ws = h + 1 bits. Then
        R(Y) = R(A) + 2^(2h)*T(B) + 2^h*V,   V = T(S) - R(A) - T(B) - [h != m]*2^m*S
        T(Y) = T(A) + 2^(2h)*T(B) + 2^h*V,   V = R(~S) - T(A) - T(B) - [h != m]*2^m*S
                                                  - (2^ws - 1)
    where ~S = 2^ws - 1 - S, so R(~S) = T(S) - S + 2^ws - 1. V is A*~B + B*~A, plus B for R,
    so 0 <= V < 2^(w + 1), and V is formed mod 2^(w + 1). R(A) and T(B) go straight onto Q's
    halves; S and the S child's output take wires from J, a list of clean junk wires, and
    keep values fixed by Y. Leaves of width at most thr are the XOR rows."""
    w = len(Y)
    if w <= thr:
        _rows(C, Y, Q, kind == "R", K)
        return
    h = (w + 1) // 2
    m = w - h
    A, B = Y[:h], Y[h:]
    _karatsuba(C, A, Q[:2 * h], J, thr, kind, K)
    _karatsuba(C, B, Q[2 * h:], J, thr, "T", K)
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
    _karatsuba(C, S, W, J, thr, "R" if kind == "T" else "T", K)
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
    # the circuit's: 636 puts secp256k1 at 4,096 wires.
    KEEP = 636
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

    def modadd(self, C, a, b, ctrl):
        """Mod.modadd with _reduce's carries above f's window held by a Prefix, so it runs
        beside the wires compute leaves live. The same flag, the same write, the same CCX."""
        f, w = self.f, self.b
        top, t = C.alloc(2)
        add(C, a, b, top, ctrl)
        cs = _carries(C, b[:w + 1], f)
        carry = Prefix(C, b[w:])
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
        _lend(C, E + T + [t, c, k] + inc,
              lambda: _karatsuba(C, Y, lo + hi, junk, self.thr, "R", keep))
        assert not junk and not keep
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
