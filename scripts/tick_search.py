"""Local search over the add's explicit block sizes on the tick bench (steady ticks + LAM * CCX).

    python3.12 -I -B scripts/tick_search.py fwd|back ITER [LAM] [SEED]
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qec import two_bit
import tick_bench


def measure():
    C, block, bounds = tick_bench.build(3)
    C, bounds, _ = tick_bench.post_pass(C, block, bounds)
    end, _ = tick_bench.schedule(C.prims, C.num_wires)
    fins = [max(end[bounds[j]:bounds[j + 1]]) for j in range(3)]
    ccx = sum(1 for k in C.prims[4 * bounds[0]:4 * bounds[1]:4] if k == 2)
    return (fins[-1] - fins[0]) / 2, ccx, C.num_wires


def mutate(sizes, rng):
    s = list(sizes)
    op = rng.random()
    i = rng.randrange(len(s))
    if op < 0.4 and len(s) > 1:
        j = rng.randrange(len(s))
        if i != j and s[i] > 1:
            s[i] -= 1
            s[j] += 1
    elif op < 0.7 and s[i] > 1:
        k = rng.randrange(1, s[i])
        s[i:i + 1] = [k, s[i] - k]
    elif len(s) > 1:
        j = i + 1 if i + 1 < len(s) else i - 1
        a, b = min(i, j), max(i, j)
        s[a:b + 1] = [s[a] + s[b]]
    return s


def main():
    which, iters = sys.argv[1], int(sys.argv[2])
    lam = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    n = 255
    names = {"fwd": ["WIDE_ADD_FSIZES"], "back": ["WIDE_ADD_BSIZES"],
             "both": ["WIDE_ADD_FSIZES", "WIDE_ADD_BSIZES"]}[which]
    cur = {}
    for name in names:
        v = getattr(two_bit, name)
        if not v:
            v = (two_bit._ramp_sizes(n, two_bit.WIDE_ADD, two_bit.WIDE_ADD_RAMP) if name == "WIDE_ADD_FSIZES"
                 else two_bit._ramp_sizes(n, *two_bit.WIDE_ADD_BACK))
        cur[name] = list(v)
        setattr(two_bit, name, tuple(v))
    t, c, w = measure()
    best = t + lam * c
    print("start", t, c, w, best, flush=True)
    for it in range(iters):
        name = rng.choice(names)
        cand = mutate(cur[name], rng)
        setattr(two_bit, name, tuple(cand))
        t2, c2, w2 = measure()
        score = t2 + lam * c2
        if w2 <= w and score <= best:
            if score < best:
                print(it, name, t2, c2, round(score, 1), tuple(cand), flush=True)
            cur[name], best = cand, score
        else:
            setattr(two_bit, name, tuple(cur[name]))
    for name in names:
        print("final", name, tuple(cur[name]))
    print("final", measure())


if __name__ == "__main__":
    main()
