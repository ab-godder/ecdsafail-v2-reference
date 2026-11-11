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
from .arith import (_carries, _uncarry, add, add_const, add_lookahead, cadd_lookahead, fan_into, gt_into,
                    gt_lookahead, prefix_xor)
from .two_bit import Prefix

# Add the square's rows as lookahead trees, on ROW_TEMPS clean wires beside the fold's own.
ROW_TREE = True
ROW_TEMPS = 30
# Fan each row's Gray-code step out through at most this many copies on clean temps, and undo the walk
# by a prefix network; 0 runs both as CX chains.
ROW_FAN_COPIES = 63
# Run modadd's controlled add and comparison as lookahead trees, on MODADD_TEMPS clean wires and
# MODADD_COPIES copies of the control.
MODADD_TREE = True
MODADD_TEMPS = 95
MODADD_COPIES = 15
# Run the fold's adds as add_lookahead trees, on FOLD_TEMPS clean wires and the fold's clean junk.
FOLD_TREE = True
FOLD_TEMPS = 126
# Split the rows once (Karatsuba): the pairs inside each half of Y add side by side on the two
# halves of Q, and the pairs across the halves come from the rows of the halves' sum, a quarter
# fewer row bits in all. The sum's rows and the sum stay on junk wires.
SPLIT_ROWS = True


def _add(C, a, b, cin=None, temps=()):
    """b += a mod 2^len(b), plus 1 if cin holds 1: add_lookahead's tree on temps when FOLD_TREE."""
    if FOLD_TREE:
        add_lookahead(C, a, b, None, list(temps), tree=True, cin=cin)
    else:
        add(C, a, b, cin=cin)


def _fan_copies(targets, cap):
    """The number of copies, 2^j - 1 <= cap, that fans src into targets in the fewest CX layers:
    j to make them and ceil(targets / 2^j) to fan out (unmaking them overlaps the add)."""
    best, j = (targets, 0), 1
    while (1 << j) - 1 <= cap:
        best = min(best, (j + -(-targets // (1 << j)), (1 << j) - 1))
        j += 1
    return best[1]


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
    def __init__(self, mod):
        self.mod, self.n = mod, mod.n
        self.f = f = (1 << mod.n) - mod.q
        self.b, self.m = f.bit_length(), (f * f).bit_length()
        # The low digit is +1 whatever f is, so its carry in can hold Z's + 1.
        self.digits = [(0, 1)] + signed_digits(f - 1)
        self.base_junk = 2 * mod.n + self.b + self.m + 3
        self.half = (mod.n + 1) // 2
        # The cross rows' register, 2 * half + 1 wires.
        self.split = SPLIT_ROWS and mod.n >= 4
        self.junk = self.base_junk + (2 * self.half + 1 if self.split else 0)

    @staticmethod
    def applies(q):
        f = (1 << q.bit_length()) - q
        return f * f <= q

    def modadd(self, C, a, b, ctrl):
        """Mod.modadd with _reduce's carries above f's window held by a Prefix, so it runs
        beside the wires compute leaves live. The same flag and the same write; with MODADD_TREE the
        add and the comparison run as lookahead trees."""
        f, w = self.f, self.b
        top, t = C.alloc(2)
        if MODADD_TREE:
            temps, copies = C.alloc(MODADD_TEMPS), C.alloc(MODADD_COPIES)
            cadd_lookahead(C, a, b, top, ctrl, temps, copies)
            C.free(copies)
            C.free(temps)
        else:
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
        if MODADD_TREE:
            temps = C.alloc(MODADD_TEMPS)
            gt_lookahead(C, a, b, top, ctrl, temps)
            C.free(temps)
        else:
            gt_into(C, a, b, top, ctrl)
        C.free([top, t])

    def _pair_rows(self, C, Y, Q, temps, cin=None):
        """Q += the sum over i < j of (y_i ^ y_j)*2^(i + j), with Q below 2^len(Y) on entry. After
        row i, Q < 2^(len(Y) + i + 1), so each carry lands on a wire no earlier row reached. Y's wires
        walk through their Gray code. The rows add in series, each with add_lookahead on clean temps.
        cin: a clean wire for the adds' carry in, so that two runs side by side share no wire."""
        n = len(Y)
        for i in range(n - 1):
            if ROW_FAN_COPIES:
                # The tree takes temps from the front, so the copies sit at the back.
                k = _fan_copies(n - 1 - i, min(ROW_FAN_COPIES, len(temps)))
                fan_into(C, Y[i], Y[i + 1:], temps[len(temps) - k:])
            else:
                for w in Y[i + 1:]:
                    C.cx(Y[i], w)
            add_lookahead(C, Y[i + 1:], Q[2 * i + 1:n + i], Q[n + i], temps, tree=ROW_TREE, cin=cin)
        if ROW_FAN_COPIES:
            prefix_xor(C, Y)
        else:
            for i in range(1, n):
                C.cx(Y[i - 1], Y[i])

    def _rows(self, C, Y, Q, temps):
        """Q = Y*(2^n - Y) on clean wires, as one run of pair rows on Y with ROW_TEMPS more clean
        wires when ROW_TREE."""
        extra = C.alloc(ROW_TEMPS) if ROW_TREE else []
        for y, q in zip(Y, Q):
            C.cx(y, q)
        self._pair_rows(C, Y, Q, list(temps) + extra)
        C.free(extra)

    def _rows_split(self, C, Y, Q, temps, K):
        """Q = Y*(2^n - Y) on clean wires, with K clean on entry and fixed by Y on exit.

        Write U(Z) for Z's pair rows, the sum over i < j of (z_i ^ z_j)*2^(i + j) = Z*(2^k - 1 - Z)
        for a k-bit Z, and Y = Y0 + 2^h*Y1 with Y0 of h bits and Y1 of h1 = n - h <= h. The pairs
        inside Y0 and inside Y1 are U(Y0) and 2^(2h)*U(Y1), on disjoint halves of Q, so the two runs
        of rows share no wire and run side by side. The pairs across the halves are 2^h*X with
        X = Y0*(2^h1 - 1 - Y1) + Y1*(2^h - 1 - Y0), and for S = Y0 + Y1 as an (h + 1)-bit number
        U(S) = S*(2^(h+1) - 1) - Y0^2 - Y1^2 - 2*Y0*Y1 gives
        X = U(S) - U(Y0) - U(Y1) - S*(2^h - 2^h1 + 1).
        Q starts as Y0, so after the halves' rows Q = Y0 + U(Y0) + 2^(2h)*U(Y1), short of
        Y + U(Y) by 2^h*(X + Y1). K takes U(S), less Q's low 2h bits, less Q's high bits, less
        S*2^h1 when h1 < h: X + Y1, which is at least 0 and below 2^(2h + 1). It then adds into Q
        from bit h. S sits on h + 1 scratch wires and is cleared again; K keeps X + Y1."""
        n, h = self.n, self.half
        S = C.alloc(h + 1)
        Y0, Y1, lo, hi = Y[:h], Y[h:], Q[:2 * h], Q[2 * h:]
        extra = C.alloc(ROW_TEMPS) if ROW_TREE else []
        cin, *temps = list(temps) + extra
        spare = temps + list(K) + list(S)
        cut = len(spare) // 2
        for y, q in zip(Y0, lo):
            C.cx(y, q)
        self._pair_rows(C, Y0, lo, spare[1:cut], spare[0])
        self._pair_rows(C, Y1, hi, spare[cut + 1:], spare[cut])
        start = len(C.prims)
        for y, s in zip(Y0, S):
            C.cx(y, s)
        pad = C.alloc(h - len(Y1))
        add_lookahead(C, Y1 + pad, S[:h], S[h], temps, tree=True, cin=cin)
        made = C.prims[start:]
        self._pair_rows(C, S, K, temps, cin)

        def take(addend, target):
            fill = C.alloc(len(target) - len(addend))
            for w in target:
                C.x(w)
            add_lookahead(C, list(addend) + fill, target, None, temps, tree=True, cin=cin)
            for w in target:
                C.x(w)
            C.free(fill)

        take(lo, K)
        take(hi, K)
        if len(Y1) < h:
            take(S, K[len(Y1):])
        for i in reversed(range(0, len(made), 4)):
            op, a, b, c = made[i:i + 4]
            (C.x(a) if op == 0 else C.cx(a, b) if op == 1 else C.ccx(a, b, c))
        C.free(pad)
        # S is clean again: it pads K to the width of Q above bit h, and the rest joins the temps.
        target = Q[h:]
        fill = len(target) - len(K)
        add_lookahead(C, list(K) + S[:fill], target, None, temps + S[fill:], tree=True, cin=cin)
        C.free(extra)
        C.free(S)

    def _multiple(self, C, a, T, cin=None, temps=()):
        """T += f*a mod 2^len(T), plus 1 if cin holds 1. temps: clean wires for the adds' trees."""
        for s, sign in self.digits:
            target = T[s:]
            if not target:
                continue
            pad = C.alloc(max(0, len(target) - len(a)))
            b = a[:len(target)] + pad
            if sign > 0:
                _add(C, b, target, cin, temps)
            else:
                for w in target:
                    C.x(w)
                _add(C, b, target, None, temps)
                for w in target:
                    C.x(w)
            C.free(pad)
            cin = None

    def compute(self, C, Y, J):
        """J clean on entry. On exit J[:n] holds Y^2 mod q, canonical, for any n-bit Y.
        The other wires of J hold values fixed by Y, for the reverse of this block."""
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:self.base_junk]
        own = J[2 * n:self.base_junk]
        # E, T, t, c and k are clean until the rows are done
        if self.split:
            self._rows_split(C, Y, lo + hi, own, J[self.base_junk:])
        else:
            self._rows(C, Y, lo + hi, own)
        fold = C.alloc(FOLD_TEMPS) if FOLD_TREE else []
        # hi = Y - hi = ~(hi + ~Y), and lo = ~lo.
        for w in Y:
            C.x(w)
        _add(C, Y, hi, None, own + fold)
        for w in Y:
            C.x(w)
        for w in hi + lo:
            C.x(w)
        (one,) = C.alloc(1)
        C.x(one)
        self._multiple(C, hi, lo + E, one, T + [t, c, k] + fold)
        C.x(one)
        C.free([one])
        self._multiple(C, E, T, None, [t, c, k] + fold)
        high = lo[m:]
        if FOLD_TREE:
            add_lookahead(C, T, lo[:m], t, [c, k] + fold, tree=True)
        else:
            add(C, T, lo[:m], cout=t)
        C.free(fold)
        if high:
            carry = Prefix(C, high)
            C.ccx(t, carry.top, c)
            carry.write(t)
        else:
            C.cx(t, c)
        # With c = 0, lo >= f. Take f off the low bits, and the borrow off the rest.
        C.x(c)
        add_const(C, lo[:b], (1 << b) - f, c, k)
        C.cx(c, k)
        high = lo[b:]
        for w in high:
            C.x(w)
        Prefix(C, high).write(k)
        for w in high:
            C.x(w)
        C.x(c)
