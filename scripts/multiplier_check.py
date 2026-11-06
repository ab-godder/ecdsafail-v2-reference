"""The multiplier on its whole domain, which the toy point-add check cannot reach:
that check sees about 2q of the q^2 inputs to M. Runs M and its inverse at gate
level on every canonical (x, y) for every odd prime up to 131, and rejects
single-gate mutations at q = 59. There the offset is not 0, so the offset
residues' constant chains have gates in them, unlike at q = 127."""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod
from qec.builder import CCX, CX, Circuit
from qec.point_add import DialogMultiplier

PRIME_LIMIT, MUTATION_Q, TRIALS = 131, 59, 200


def build(q, inverse):
    mod = Mod(q)
    C = Circuit()
    x, y = C.alloc(mod.n), C.alloc(mod.n)
    x_in, y_in = list(x), list(y)
    DialogMultiplier(C, mod).apply(x, y, inverse=inverse)
    return C, x_in, y_in, x, y


def failures(q, inverse, C, ops, x_in, y_in, x, y):
    inputs = [(a, b) for a in range(q) for b in range(q)]
    lanes = len(inputs)
    init = {**sim.pack([a for a, _ in inputs], x_in), **sim.pack([b for _, b in inputs], y_in)}
    state = sim.run(ops, C.num_wires, init, lanes)
    got_x, got_y = sim.unpack(state, x, lanes), sim.unpack(state, y, lanes)
    live = {*x, *y}
    dirty = 0
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    bad = []
    for i, (a, b) in enumerate(inputs):
        want = b if a == 0 else (pow(a, -1, q) if inverse else a) * b % q
        if got_x[i] != a or got_y[i] != want or dirty >> i & 1:
            bad.append((a, b))
    return lanes, bad


def main():
    primes = [p for p in range(3, PRIME_LIMIT + 1, 2) if all(p % d for d in range(3, int(p**0.5) + 1, 2))]
    lanes, failing = 0, {}
    for q in primes:
        for inverse in (False, True):
            C, *regs = build(q, inverse)
            k, bad = failures(q, inverse, C, C.ops, *regs)
            lanes += k
            if bad:
                failing[f"q={q},inverse={inverse}"] = {"failures": len(bad), "first": bad[:3]}
    rng = random.Random(7)
    C, *regs = build(MUTATION_Q, False)
    ccx = [i for i in range(len(C)) if C.ops[4 * i] == CCX]
    caught = 0
    for _ in range(TRIALS):
        ops = C.ops[:]
        i = rng.choice(ccx)
        ops[4 * i] = CX
        ops[4 * i + 2] = ops[4 * i + 3]
        caught += bool(failures(MUTATION_Q, False, C, ops, *regs)[1])
    print(json.dumps({"primes": len(primes), "lanes": lanes, "failing": failing,
                      "mutation_q": MUTATION_Q, "trials": TRIALS, "caught": caught}))
    sys.exit(1 if failing else 0)


if __name__ == "__main__":
    main()
