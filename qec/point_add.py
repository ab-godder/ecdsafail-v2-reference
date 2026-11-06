"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, _maj, add_const, cswap, load_const
from .builder import Circuit
from .curve import Curve
from .model import iteration_bound, width_schedule


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for x != 0, M(0, y) = (0, y); see model.py."""

    def __init__(self, C, mod):
        self.C, self.mod = C, mod
        q, n = mod.q, mod.n
        self.iters = iteration_bound(q)
        self.widths = width_schedule(q, self.iters)
        self._euclid_steps = {}
        self._bezout_step = None
        self.block = None

    def _euclid_step(self, w):
        """One Euclid step on (U, Z, e) with v = Z + e and u odd. b0 = v mod 2 is
        Z[0] ^ e, and that wire keeps it. When v is odd, Z's upper bits are v's
        whatever e is, so one carry chain over U[1:w] and ~Z[1:w] gives the swap
        bit, the difference and the smaller value. |u - v|/2 is the complement of
        that sum when u <= v, and the sum plus one otherwise, so the + 1 stays
        pending in e. An even v just shifts, carrying e along."""
        if w not in self._euclid_steps:
            C = self.C

            def step(U, Z, E, D):
                top, g1, e_next = D
                e, g0 = E[0], Z[0]
                a, b = U[1:w], Z[1:w]
                C.cx(e, g0)
                for x in b:
                    C.x(x)
                (cin,) = C.alloc(1)
                carries = [cin] + a[:-1]
                for c, x, y in zip(carries, b, a):
                    _maj(C, c, x, y)
                if a:
                    C.ccx(g0, a[-1], g1)
                for c, x, y in reversed(list(zip(carries, b, a))):
                    C.ccx(c, x, y)
                    C.ccx(g0, c, x)
                    C.cx(y, c)
                    C.cx(y, x)
                    # g1 implies g0, so x holds the sum and x ^ c is y ^ ~Z[i].
                    C.cx(c, x)
                    C.x(x)
                    C.ccx(g1, x, y)
                    C.x(x)
                    C.cx(c, x)
                C.free([cin])
                for x in b:
                    C.x(x)
                    C.cx(g1, x)
                C.x(g0)
                C.ccx(g0, e, e_next)
                C.x(g0)
                C.cx(g1, e_next)
                Z[:] = Z[1:] + [top]
                E[:] = [e_next]
                D[:] = [g0, g1, e]

            n = self.mod.n
            self._euclid_steps[w] = C.record_shape(step, [n, n, 1, 3])
        return self._euclid_steps[w]

    def _bezout(self):
        """One inverse Euclid step on (R, W) = y*(u, 2v + u) mod q, in offset residues.
        These run from (y, y), since Euclid ends on (u, v) = (1, 0). They end on
        (0, 2xy), or stay at (y, y) when x = 0."""
        if self._bezout_step is None:
            C, mod, n = self.C, self.mod, self.mod.n

            def step(g, R, W):
                g0, g1 = g
                cswap(C, g1, R, W)
                mod.offset_dbl(C, W)
                # g0 = 0 inverts v' = v/2, so W = 2W - R. Complementing an
                # offset residue negates it, and the complement of 0's residue
                # still works as an addend.
                C.x(g0)
                for r in R:
                    C.cx(g0, r)
                mod.offset_add(C, R, W)
                for r in R:
                    C.cx(g0, r)
                C.x(g0)

            self._bezout_step = C.record_shape(step, [2, n, n])
        return self._bezout_step

    def _is_zero_into(self, bits, t):
        C = self.C
        (one,) = C.alloc(1)
        C.x(one)
        eq_const_into(C, one, bits, 0, t)
        C.x(one)
        C.free([one])

    def _init_u(self, U, z):
        """U = q, or 1 when x = 0 so that every Euclid step idles on registers
        that fit the shrinking widths."""
        C, q = self.C, self.mod.q
        load_const(C, U, q)
        load_const(C, U, q ^ 1, ctrl=z[0])

    def _body(self, Xr, Yr):
        C, n, T = self.C, self.mod.n, self.iters
        (z,) = C.alloc(1)
        self._is_zero_into(Xr, z)
        U = C.alloc(n)
        Z = list(Xr)
        E = C.alloc(1)
        I = C.run(self._init_u, [U, [z]])
        G = C.alloc(3 * T)

        def euclid(U, Z, E, G):
            for i in range(T):
                D = G[3 * i:3 * i + 3]
                C.apply(self._euclid_step(self.widths[i]), [U, Z, E, D])
                G[3 * i:3 * i + 3] = D

        Euc = C.run(euclid, [U, Z, E, G])
        W = C.alloc(n)
        mod, offset = self.mod, self.mod.offset()
        add_const(C, Yr, offset)
        for y, w in zip(Yr, W):
            C.cx(y, w)
        bz = self._bezout()
        for i in reversed(range(T)):
            C.apply(bz, [G[3 * i:3 * i + 2], Yr, W])
        # (R + W)/2 is xy for x != 0 and y for x = 0.
        mod.offset_add(C, Yr, W)
        C.apply(C.record_shape(lambda s: mod.offset_dbl(C, s), [n]), [W], inverse=True)
        for y, w in zip(Yr, W):
            C.ccx(z, w, y)
        C.x(z)
        load_const(C, Yr, offset, ctrl=z)
        C.x(z)
        add_const(C, W, -offset)
        C.apply(Euc, [U, Z, E, G], inverse=True)
        C.apply(I, [U, [z]], inverse=True)
        self._is_zero_into(Xr, z)
        C.free(U + E + G + Yr + [z])
        Yr[:] = W

    def apply(self, Xr, Yr, inverse=False):
        if self.block is None:
            self.block = self.C.record(self._body, [Xr, Yr])
        self.C.apply(self.block, [Xr, Yr], inverse)


def eq_const_into(C, ctrl, bits, value, target):
    """target ^= ctrl AND [bits == value]."""
    flips = [w for i, w in enumerate(bits) if not value >> i & 1]
    for w in flips:
        C.x(w)
    chain = []
    acc = ctrl
    for w in bits[:-1]:
        (a,) = C.alloc(1)
        C.ccx(acc, w, a)
        chain.append((acc, w, a))
        acc = a
    C.ccx(acc, bits[-1], target)
    for acc, w, a in reversed(chain):
        C.ccx(acc, w, a)
    C.free([a for _, _, a in chain])
    for w in flips:
        C.x(w)


def encode(P, n):
    x, y = (0, 0) if P is None else P
    return x | y << n


class ControlledPointAdd:
    def __init__(self, curve: Curve, A):
        assert A is not None and curve.on_curve(A)
        self.curve, self.A = curve, A
        self.mod = Mod(curve.q)
        self.n = self.mod.n
        neg = curve.neg
        twoA = curve.add(A, A)
        # (exceptional input, its sum with A); the generic kernel divides by zero on these.
        self.exceptions = [(None, A), (A, twoA), (neg(A), None), (neg(twoA), neg(A))]
        enc = [encode(e, self.n) for e, _ in self.exceptions]
        assert len(set(enc)) == 4, "addend order must exceed 3"

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * Y^2 mod q, Y untouched. Y^2 is computed once as an integer
        with no modular steps, reduced, subtracted once, and reversed."""
        C, mod, n = self.C, self.mod, self.n
        pool_size = mod.square_pool_size()
        Z, pool = C.alloc(n), C.alloc(pool_size)
        if self._modadd_blk is None:
            self._modadd_blk = C.record_shape(lambda k, a, b: mod.modadd(C, a, b, ctrl=k[0]), [1, n, n])
            self._square = C.record_shape(lambda y, z, p: mod.square_into(C, y, z, p), [n, n, pool_size])
        C.apply(self._square, [Yr, Z, pool])
        C.apply(self._modadd_blk, [[gen], Z, Xr], inverse=True)
        C.apply(self._square, [Yr, Z, pool], inverse=True)
        C.free(Z + pool)

    def _kernel(self, gen, Xr, Yr):
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A
        mod.cmodadd_const(C, gen, q - a, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        self.M.apply(Xr, Yr, inverse=True)
        mod.cmodadd_const(C, gen, 3 * a, Xr)
        self._mulsub(gen, Yr, Xr)
        self.M.apply(Xr, Yr)
        mod.cneg(C, gen, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        mod.cmodadd_const(C, gen, a, Xr)

    def build(self, hierarchical=False):
        self.C = C = Circuit(hierarchical)
        n = self.n
        (c,) = C.alloc(1)
        Xr, Yr = C.alloc(n), C.alloc(n)
        declared = Xr + Yr
        self.M = DialogMultiplier(C, self.mod)
        self._modadd_blk = None

        flags = C.alloc(4)
        for f, (e, _) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e, n), f)
        e, gen = C.alloc(2)
        for f in flags:
            C.cx(f, e)
        C.x(e)
        C.ccx(c, e, gen)
        for f, (e_in, e_out) in zip(flags, self.exceptions):
            load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f)

        self._kernel(gen, Xr, Yr)

        C.ccx(c, e, gen)
        C.x(e)
        for f in flags:
            C.cx(f, e)
        C.free([e, gen])
        for f, (_, e_out) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e_out, n), f)
        C.free(flags)

        _restore_layout(C, Xr + Yr, declared)
        self.layout = {"control": [c], "x": declared[:n], "y": declared[n:]}
        return C


def _restore_layout(C, cur, declared):
    """SWAP logical bit i from wire cur[i] onto declared[i]. A displaced wire
    holds either another live bit, which is tracked, or a zero ancilla."""
    cur = list(cur)
    where = {w: i for i, w in enumerate(cur)}
    for i, d in enumerate(declared):
        s = cur[i]
        if s == d:
            continue
        C.swap(s, d)
        j = where.pop(d, None)
        del where[s]
        if j is not None:
            cur[j] = s
            where[s] = j
        cur[i] = d
        where[d] = i
