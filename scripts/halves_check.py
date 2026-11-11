"""Gate-level check of HalfSquares, the square path X -= gen*(Y^2 - K) mod q by one Karatsuba
level taken mod q.

Builds the two groups, their gen-controlled modadds and the unbuilds exactly as
HalfSquares.mulsub does, and checks on every lane that X becomes (X - gen*Y^2 + gen*K) mod q,
that gen and Y are unchanged and that every other wire is zero. Y is any n-bit value and X
canonical. Small widths run every (gen, Y, X), middle widths every Y with random X, larger
widths and secp256k1 random lanes plus edge values of Y. Every prime of each small width where
the class applies is run; larger widths take the smallest and largest f. Each width also runs
with Karatsuba leaves inside the half squares (HalfSquares.THR_AC, THR_B) small enough to recurse
there, and secp256k1 with the leaf widths the trade curve uses."""
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
from qec.square import HalfSquares


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


def primes(n, count=None):
    """Every prime 2^n - f the class applies to, or the count smallest and largest f. Its bounds
    need f below 2^(n/2 + 1), so the search for the largest starts there."""
    ok = lambda f: is_prime((1 << n) - f) and HalfSquares.applies((1 << n) - f)
    # The bounds fail once f*2^(n/2) nears 2^n; find the last f where they hold by bisection,
    # checking each candidate below it again in full.
    lo_f, hi_f = 1, 1 << (n // 2 + 1)
    while hi_f - lo_f > 2:
        mid = (lo_f + hi_f) // 2 | 1
        lo_f, hi_f = (mid, hi_f) if HalfSquares.bounds((1 << n) - mid) else (lo_f, mid)
    top = hi_f + 64
    if count is None:
        fs = [f for f in range(1, top, 2) if ok(f)]
    else:
        lo = [f for f in range(1, top, 2) if ok(f)][:count] if top < 1 << 16 else []
        if not lo:
            f = 1
            while len(lo) < count and f < top:
                if ok(f):
                    lo.append(f)
                f += 2
        hi, f = [], top | 1
        while len(hi) < count and f > 0:
            if ok(f):
                hi.append(f)
            f -= 2
        fs = lo + hi
    return sorted({(1 << n) - f for f in fs})


def edge_ys(n, q):
    h = n // 2
    vals = {0, 1, 2, q - 1, q, q + 1, (1 << n) - 1, (1 << n) - 2, (1 << h) - 1, 1 << h,
            ((1 << h) - 1) << h, (((1 << h) - 1) << h) | ((1 << h) - 1), (1 << (n - 1)),
            (1 << (n - 1)) - 1}
    return sorted(v % (1 << n) for v in vals)


def pack(values, wires):
    """sim.pack in linear time: each wire's lane mask from one binary string."""
    return {w: int("".join("1" if v >> i & 1 else "0" for v in reversed(values)) or "0", 2)
            for i, w in enumerate(wires)}


def unpack(state, wires, lanes):
    """sim.unpack in linear time."""
    vals = [0] * lanes
    for i, w in enumerate(wires):
        bits = bin(state[w])[2:].zfill(lanes)[::-1]
        for lane in range(lanes):
            if bits[lane] == "1":
                vals[lane] |= 1 << i
    return vals


def run(q, lanes, drop=None, thr=(None, None)):
    """drop, a fraction, removes that CCX of the circuit: a control that must fail. thr is
    (THR_AC, THR_B)."""
    mod = Mod(q)
    n = mod.n
    hs = HalfSquares(mod)
    hs.THR_AC, hs.THR_B = thr
    C = Circuit()
    (gen,) = C.alloc(1)
    Y, X = C.alloc(n), C.alloc(n)
    hs.mulsub(C, gen, Y, X)
    if drop is not None:
        ccx = [k for k in range(len(C.prims) // 4) if C.prims[4 * k] == 2]
        k = ccx[int(len(ccx) * drop)]
        del C.prims[4 * k:4 * k + 4]
    k = len(lanes)
    init = {}
    init.update(pack([g for g, _, _ in lanes], [gen]))
    init.update(pack([y for _, y, _ in lanes], Y))
    init.update(pack([x for _, _, x in lanes], X))
    state = sim.run(C.ops, C.num_wires, init, k)
    got_g = unpack(state, [gen], k)
    got_y = unpack(state, Y, k)
    got_x = unpack(state, X, k)
    live = {gen, *Y, *X}
    dirty = 0
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    fails = [i for i, (g, y, x) in enumerate(lanes)
             if (got_g[i], got_y[i], got_x[i]) != (g, y, (x - g * (y * y - hs.K)) % q)
             or dirty >> i & 1]
    return fails, C.num_wires, C.ops[0::4].count(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int, nargs="*", default=[8, 10])
    ap.add_argument("--all-y", type=int, nargs="*", default=[12, 14, 16])
    ap.add_argument("--sampled", type=int, nargs="*", default=[18, 20, 24, 32, 48, 64])
    ap.add_argument("--secp", type=int, default=256, help="random lanes at secp256k1")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    bad, t0 = 0, time.time()
    for n in args.exhaustive + args.all_y + args.sampled:
        for q in primes(n, None if n in args.exhaustive + args.all_y else 2):
            if n in args.exhaustive:
                lanes = [(g, y, x) for g in (0, 1) for y in range(1 << n) for x in range(q)]
            elif n in args.all_y:
                lanes = [(g, y, rng.randrange(q)) for g in (0, 1) for y in range(1 << n)]
                lanes += [(1, y, x) for y in edge_ys(n, q) for x in (0, 1, q - 1)]
            else:
                lanes = [(g, rng.randrange(1 << n), rng.randrange(q)) for g in (0, 1) for _ in range(512)]
                lanes += [(1, y, x) for y in edge_ys(n, q) for x in (0, 1, q - 1, rng.randrange(q))]
            h = n // 2
            thrs = [(None, None), (3, 3), (None, max(3, h // 2)), (max(3, h // 2 + 1), None)]
            if n in args.exhaustive and n > 8:
                thrs = thrs[:2]
            for thr in thrs:
                fails, wires, ccx = run(q, lanes, thr=thr)
                bad += len(fails)
                print(f"n={n} q={q} f={(1 << n) - q} thr={thr} lanes={len(lanes)} "
                      f"failures={len(fails)} wires={wires} ccx={ccx}", flush=True)
    if args.secp:
        q = SECP256K1.q
        lanes = [(g, rng.randrange(1 << 256), rng.randrange(q)) for g in (0, 1) for _ in range(args.secp // 2)]
        lanes += [(1, y, x) for y in edge_ys(256, q) for x in (0, 1, q - 1, rng.randrange(q))]
        for thr in ((None, None), (None, 65), (66, 65), (43, 43)):
            fails, wires, ccx = run(q, lanes, thr=thr)
            bad += len(fails)
            print(f"secp256k1 thr={thr} lanes={len(lanes)} failures={len(fails)} wires={wires} "
                  f"ccx={ccx}", flush=True)
    # Controls: one CCX dropped at each of several places, run on every (Y, X) with gen = 1 at a
    # 10-bit prime, so every gate of both builds and the modadds runs on its whole domain. A
    # control that passes there is an equivalent circuit, which is reported, not counted as caught.
    caught = controls = 0
    q = primes(10)[1]
    lanes = [(1, y, x) for y in range(1 << 10) for x in range(q)]
    for k in range(1, 20):
        fails, _, _ = run(q, lanes, drop=k / 20)
        controls += 1
        caught += bool(fails)
        if not fails:
            print(f"control {k}/20 at q={q}: equivalent on the whole domain", flush=True)
    print(f"controls caught {caught} of {controls}")
    print(f"total failures {bad}, {time.time() - t0:.1f}s")
    sys.exit(1 if bad or caught < controls - 2 else 0)


if __name__ == "__main__":
    main()
