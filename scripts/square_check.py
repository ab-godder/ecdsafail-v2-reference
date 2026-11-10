"""Gate-level check of the square path X -= gen * Y^2 mod q at several KEEP values.

Builds the square block, the gen-controlled modadd (borrowing the junk wires the square
leaves 0) and the square's reverse exactly as ControlledPointAdd._square_mulsub does, and checks on every lane that X becomes
(X - gen*Y^2) mod q, that gen and Y are unchanged and that every other wire is zero. Small
widths run every (gen, Y, X) with Y any n-bit value and X canonical; larger widths run every Y
(or random Y) with random canonical X. KEEP sweeps no kept carries, a few, part of the leaves'
carries, all of them, and more than all (capped)."""
import argparse
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.curve import SECP256K1
from qec.square import SquareFold, karatsuba_carries


def is_prime(v):
    if v < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % p == 0:
            return v == p
    d, s = v - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, v)
        if x in (1, v - 1):
            continue
        for _ in range(s - 1):
            x = x * x % v
            if x == v - 1:
                break
        else:
            return False
    return True


def fold_primes(n, count):
    """Primes q = 2^n - f with f*f <= q: the smallest f and the largest f that apply."""
    fs = [f for f in range(1, 1 << (n // 2 + 1), 2) if is_prime((1 << n) - f)
          and SquareFold.applies((1 << n) - f)]
    pick = fs[:count] + fs[-count:]
    return sorted({(1 << n) - f for f in pick})


def run(q, keep, lanes):
    """lanes: list of (gen, Y, X). Returns (failures, wires, ccx, kept)."""
    mod = Mod(q)
    n = mod.n
    Fold = type("Fold", (SquareFold,), {"KEEP": keep})
    fold = Fold(mod)
    C = Circuit()
    (gen,) = C.alloc(1)
    Y, X = C.alloc(n), C.alloc(n)
    sq = C.record_shape(lambda y, j: fold.compute(C, y, j), [n, fold.junk])
    zeros = fold.zero_slots
    madd = C.record_shape(lambda k, a, b, z: fold.modadd(C, a, b, k[0], z), [1, n, n, len(zeros)])
    J = C.alloc(fold.junk)
    C.apply(sq, [Y, J])
    C.apply(madd, [[gen], J[:n], X, [J[i] for i in zeros]], inverse=True)
    C.apply(sq, [Y, J], inverse=True)
    C.free(J)
    k = len(lanes)
    init = {}
    init.update(sim.pack([g for g, _, _ in lanes], [gen]))
    init.update(sim.pack([y for _, y, _ in lanes], Y))
    init.update(sim.pack([x for _, _, x in lanes], X))
    state = sim.run(C.ops, C.num_wires, init, k)
    got_g = sim.unpack(state, [gen], k)
    got_y = sim.unpack(state, Y, k)
    got_x = sim.unpack(state, X, k)
    live = {gen, *Y, *X}
    dirty = 0
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    fails = [i for i, (g, y, x) in enumerate(lanes)
             if (got_g[i], got_y[i], got_x[i]) != (g, y, (x - g * y * y) % q) or dirty >> i & 1]
    return fails, C.num_wires, C.ops[0::4].count(2), fold.keep


def keeps_for(n):
    full = karatsuba_carries(n, max(3, n // 6))
    return sorted({0, 1, 2, 3, full // 2, max(0, full - 1), full, full + 5})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int, nargs="*", default=[4, 5, 6, 7, 8])
    ap.add_argument("--all-y", type=int, nargs="*", default=[10, 12, 14])
    ap.add_argument("--sampled", type=int, nargs="*", default=[24, 30, 36, 48])
    ap.add_argument("--secp", type=int, default=64, help="random lanes at secp256k1 per KEEP")
    ap.add_argument("--primes", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    bad = 0
    t0 = time.time()
    rows = []
    for n in args.exhaustive + args.all_y + args.sampled:
        for q in fold_primes(n, args.primes):
            if n in args.exhaustive:
                lanes = [(g, y, x) for g in (0, 1) for y in range(1 << n) for x in range(q)]
            elif n in args.all_y:
                lanes = [(g, y, rng.randrange(q)) for g in (0, 1) for y in range(1 << n)]
            else:
                lanes = [(g, rng.randrange(1 << n), rng.randrange(q)) for g in (0, 1) for _ in range(512)]
                lanes += [(1, v, rng.randrange(q)) for v in (0, 1, q - 1, q, (1 << n) - 1)]
            for keep in keeps_for(n):
                fails, wires, ccx, kept = run(q, keep, lanes)
                bad += len(fails)
                rows.append((n, q, keep, kept, len(lanes), len(fails), wires, ccx))
                print(f"n={n} q={q} KEEP={keep} kept={kept} lanes={len(lanes)} failures={len(fails)} "
                      f"wires={wires} ccx={ccx}", flush=True)
    if args.secp:
        q = SECP256K1.q
        lanes = [(g, rng.randrange(1 << 256), rng.randrange(q)) for g in (0, 1) for _ in range(args.secp // 2)]
        lanes += [(1, v, rng.randrange(q)) for v in (0, 1, q - 1, q, (1 << 256) - 1)]
        for keep in (0, 950, SquareFold.KEEP, karatsuba_carries(256, 42)):
            fails, wires, ccx, kept = run(q, keep, lanes)
            bad += len(fails)
            print(f"secp256k1 KEEP={keep} kept={kept} lanes={len(lanes)} failures={len(fails)} "
                  f"wires={wires} ccx={ccx}", flush=True)
    print(f"total failures {bad}, {time.time() - t0:.1f}s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
