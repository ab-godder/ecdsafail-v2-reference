"""Offset-residue checks that the gate's pinned component list does not reach:
every value at gate level for every odd prime below 130, including the
complement of the offset residue of 0, and boundary inputs at secp256k1."""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.arith import Mod
from scripts.component_check import run_case


def cases(mod, values, operands):
    """Each case runs on registers holding residue + c; an operand of q is the
    complement of the residue of 0."""
    q, n, c = mod.q, mod.n, mod.offset()
    out = {}

    def b_add(C):
        a, b = C.alloc(n), C.alloc(n)
        mod.offset_add(C, a, b)
        return [a, b], [a, b]
    out["add"] = run_case(b_add, [(a + c, s + c) for a in operands for s in values],
                          lambda a, s: (a, (s + a - 2 * c) % q + c))

    def b_cadd(C):
        g, a, b = C.alloc(1), C.alloc(n), C.alloc(n)
        mod.offset_add(C, a, b, ctrl=g[0])
        return [g, a, b], [g, a, b]
    out["add_ctrl"] = run_case(b_cadd, [(g, a + c, s + c) for g in (0, 1) for a in operands for s in values],
                               lambda g, a, s: (g, a, (s - c + g * (a - c)) % q + c))

    def b_dbl(C):
        s = C.alloc(n)
        s_in = list(s)
        mod.offset_dbl(C, s)
        return [s_in], [s]
    out["dbl"] = run_case(b_dbl, [(s + c,) for s in values], lambda s: (2 * (s - c) % q + c,))
    return out


def main():
    out, primes = {}, [p for p in range(3, 130, 2) if all(p % d for d in range(3, int(p**0.5) + 1, 2))]
    for q in primes:
        for name, r in cases(Mod(q), range(q), range(q + 1)).items():
            out[f"q={q}:{name}"] = r
    q = 2**256 - 2**32 - 977
    f, c = 2**256 - q, Mod(q).offset()
    rng = random.Random(7)
    edge = [0, 1, 2, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2, 2**255, 2**255 - 1, f, q - f, c, q - c, c + 1, q - c - 1, 2**255 - c, 2**255 - c - 1]
    edge = [e % q for e in edge] + [rng.randrange(q) for _ in range(30)]
    for name, r in cases(Mod(q), edge, edge + [q]).items():
        out[f"secp256k1:{name}"] = r
    fails = {k: v for k, v in out.items() if v["failures"]}
    print(json.dumps({"cases": len(out), "lanes": sum(v["lanes"] for v in out.values()), "failing": fails}))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
