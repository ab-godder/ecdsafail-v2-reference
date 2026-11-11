"""Prove the K1 cut: in the kernel's cmodadd_const(gen, 3a, X) at secp256k1's G, the CCX pair the
table in qec/point_add.py names may become CX from gen, on every canonical X (b < q) and either gen.

    python3.12 -I -B scripts/k1_check.py

A small reduced ordered BDD package of its own runs the call symbolically over its 257 inputs, gen
then X's bits low to high, once as built and once with the cut. For every wire, the two final
functions must agree wherever b < q holds. Controls, each of which must find a disagreement:
either half of the pair alone; the pair without the b < q bound; and one other CCX pair on the
same two controls, which is not a valid cut.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec.arith import Mod
from qec.builder import CCX, CX, X, Circuit
from qec.curve import SECP256K1, SECP256K1_G
from qec.point_add import K1_CUTS

FALSE, TRUE = 0, 1


class BDD:
    def __init__(self):
        self.nodes = [(1 << 30, 0, 0), (1 << 30, 1, 1)]
        self.unique = {}
        self.memo = {}

    def mk(self, v, lo, hi):
        if lo == hi:
            return lo
        k = (v, lo, hi)
        n = self.unique.get(k)
        if n is None:
            n = self.unique[k] = len(self.nodes)
            self.nodes.append(k)
        return n

    def var(self, v):
        return self.mk(v, FALSE, TRUE)

    def apply(self, op, a, b):
        if op == "and":
            if a == FALSE or b == FALSE:
                return FALSE
            if a == TRUE:
                return b
            if b == TRUE or a == b:
                return a
        else:
            if a == FALSE:
                return b
            if b == FALSE:
                return a
            if a == b:
                return FALSE
        if a > b:
            a, b = b, a
        k = (op, a, b)
        r = self.memo.get(k)
        if r is not None:
            return r
        va, la, ha = self.nodes[a]
        vb, lb, hb = self.nodes[b]
        v = min(va, vb)
        a0, a1 = (la, ha) if va == v else (a, a)
        b0, b1 = (lb, hb) if vb == v else (b, b)
        r = self.mk(v, self.apply(op, a0, b0), self.apply(op, a1, b1))
        self.memo[k] = r
        return r

    def xor(self, a, b):
        return self.apply("xor", a, b)

    def land(self, a, b):
        return self.apply("and", a, b)

    def lnot(self, a):
        return self.xor(a, TRUE)

    def sat(self, f):
        """One satisfying assignment {var: bit} of f, or None."""
        if f == FALSE:
            return None
        out = {}
        while f > TRUE:
            v, lo, hi = self.nodes[f]
            if hi != FALSE:
                out[v], f = 1, hi
            else:
                out[v], f = 0, lo
        return out


def call_ops(q, value):
    """The ops of mod.cmodadd_const(ctrl, value, b) alone, its ctrl wire and b's wires."""
    C = Circuit()
    (ctrl,), b = C.alloc(1), C.alloc(q.bit_length())
    Mod(q).cmodadd_const(C, ctrl, value, b)
    return list(C.prims), ctrl, b, C.num_wires


def run(bdd, ops, num, ctrl, b):
    s = [FALSE] * num
    s[ctrl] = bdd.var(0)
    for i, w in enumerate(b):
        s[w] = bdd.var(i + 1)
    for p in range(0, len(ops), 4):
        op, x, y, t = ops[p:p + 4]
        if op == X:
            s[x] = bdd.lnot(s[x])
        elif op == CX:
            s[y] = bdd.xor(s[y], s[x])
        elif op == CCX:
            s[t] = bdd.xor(s[t], bdd.land(s[x], s[y]))
        else:
            s[x], s[y] = s[y], s[x]
    return s


def below(bdd, b, q):
    """b < q as a BDD, built from the top bit down."""
    lt, eq = FALSE, TRUE
    for i in reversed(range(len(b))):
        x = bdd.var(i + 1)
        if q >> i & 1:
            lt = bdd.xor(lt, bdd.land(eq, bdd.lnot(x)))
            eq = bdd.land(eq, x)
        else:
            eq = bdd.land(eq, bdd.lnot(x))
    return lt


def eval_at(bdd, f, g, x):
    while f > TRUE:
        v, lo, hi = bdd.nodes[f]
        bit = g if v == 0 else x >> (v - 1) & 1
        f = hi if bit else lo
    return f


def concrete(ops, num, ctrl, b):
    """(ctrl, b) -> every wire's final value, by the shared-nothing lane simulator."""
    from array import array
    from qec import sim
    arr = array("i", ops)

    def at(g, x):
        init = {ctrl: g, **{w: x >> i & 1 for i, w in enumerate(b)}}
        return tuple(sim.run(arr, num, init, 1))
    return at


def cut(ops, positions, ctrl):
    out = list(ops)
    for pos in positions:
        op, x, y, t = out[4 * pos:4 * pos + 4]
        assert op == CCX and ctrl in (x, y), (pos, op, x, y)
        out[4 * pos:4 * pos + 4] = [CX, ctrl, t, t]
    return out


def differs(bdd, base, ops, num, ctrl, b, dom):
    """A lane (ctrl, b) in dom where some wire differs, or None."""
    s0 = base
    s1 = run(bdd, ops, num, ctrl, b)
    for w in range(num):
        d = bdd.land(bdd.xor(s0[w], s1[w]), dom)
        if d != FALSE:
            a = bdd.sat(d)
            return a.get(0, 0), sum(a.get(i + 1, 0) << i for i in range(len(b)))
    return None


def main():
    q, a = SECP256K1.q, SECP256K1_G[0]
    value = 3 * a % q
    positions = K1_CUTS[(q, value)]
    ops, ctrl, b, num = call_ops(q, value)
    bdd = BDD()
    base = run(bdd, ops, num, ctrl, b)
    dom = below(bdd, b, q)
    rows = {}
    rows["pair on b < q"] = differs(bdd, base, cut(ops, positions, ctrl), num, ctrl, b, dom)
    for k, pos in enumerate(positions):
        rows[f"half {pos} alone"] = differs(bdd, base, cut(ops, [pos], ctrl), num, ctrl, b, dom)
    rows["pair on every b"] = differs(bdd, base, cut(ops, positions, ctrl), num, ctrl, b, TRUE)
    pairs = [p for p in range(len(ops) // 4)
             if ops[4 * p] == CCX and {ops[4 * p + 1], ops[4 * p + 2]} == {ctrl, b[254]} and p not in positions]
    other = (pairs[0], pairs[-1])
    rows[f"other pair {other}"] = differs(bdd, base, cut(ops, other, ctrl), num, ctrl, b, dom)
    ok = rows["pair on b < q"] is None and all(v is not None for k, v in rows.items() if k != "pair on b < q")
    for k, v in rows.items():
        print(f"{k}: {'no difference' if v is None else f'differs at ctrl={v[0]}, b=0x{v[1]:x}'}")
    # The BDD run is the call: on random lanes it is b + ctrl*value mod q. Each counterexample is
    # one: the plain and cut ops, simulated concretely, differ on it.
    import random
    rng = random.Random(1)
    lanes = [(rng.randrange(2), rng.randrange(q)) for _ in range(64)] + [(1, q - 1), (1, q - value), (1, q - value - 1)]
    for g, x in lanes:
        out = sum(int(eval_at(bdd, base[w], g, x)) << i for i, w in enumerate(b))
        assert out == (x + g * value) % q, (g, x)
    plain = concrete(ops, num, ctrl, b)
    for k, v in rows.items():
        if v is not None:
            g, x = v
            assert k == "pair on every b" or x < q
            mut = cut(ops, positions if k.startswith("pair") else other if k.startswith("other")
                      else [int(k.split()[1])], ctrl)
            assert concrete(mut, num, ctrl, b)(g, x) != plain(g, x), k
    print(f"positions {positions} of {len(ops) // 4} ops; {len(bdd.nodes)} BDD nodes; {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
