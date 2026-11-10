"""Gate-level check of the multiplier at every KEEP budget, so that each boundary tick keeps all,
part or none of its chain.

For each prime q given (default: every odd prime below 600), KEEP runs from 0 to the full kept
plan; each build runs M forward and inverse on every x in 1..q-1 against six y values and counts
lanes where Y is wrong, X changed or any other wire is not zero."""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.emit import write_circuit
from qec.point_add import DialogMultiplier


def bad_lanes(q: int, ys, inverse: bool) -> int:
    pairs = [(x, y) for x in range(1, q) for y in ys]
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    DialogMultiplier(C, mod).apply(X, Y, inverse=inverse)
    assert len(set(X + Y)) == 2 * mod.n, "aliased output register"
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y})
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    xs, values = zip(*pairs)
    init = sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1])
    state = sim.run(ops, C.num_wires, init, len(pairs))
    expected = [y * (pow(x, -1, q) if inverse else x) % q for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def full_counts(m: DialogMultiplier):
    w = m.widths
    return [DialogMultiplier._kept_count(w[i], i >= 1 and w[i] >= 3 and w[i - 1] == w[i],
                                         i + 1 < m.iters and w[i + 1] == w[i] - 1, i == m.iters - 1)
            for i in range(m.iters)]


primes = [int(a) for a in sys.argv[1:]] or [
    p for p in range(3, 600, 2) if all(p % d for d in range(3, int(p ** 0.5) + 1, 2))]
saved = DialogMultiplier.KEEP
total_bad = cases = partial_cases = 0
for q in primes:
    DialogMultiplier.KEEP = 10 ** 9
    probe = DialogMultiplier(Circuit(True), Mod(q))
    full, counts = sum(probe.kept_plan()), full_counts(probe)
    ys = sorted({0, 1, 2, q - 1, q // 2, 3 % q})
    for keep in range(full + 1):
        DialogMultiplier.KEEP = keep
        plan = DialogMultiplier(Circuit(True), Mod(q)).kept_plan()
        bad = bad_lanes(q, ys, False) + bad_lanes(q, ys, True)
        total_bad += bad
        cases += 1
        partial_cases += any(0 < k < c for k, c in zip(plan, counts))
        if bad:
            print("FAIL", q, keep, bad, flush=True)
DialogMultiplier.KEEP = saved
print(f"{len(primes)} primes, {cases} builds, {partial_cases} with a partial tick, {total_bad} bad lanes")
sys.exit(1 if total_bad else 0)
