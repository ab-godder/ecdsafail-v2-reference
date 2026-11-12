"""The Bezout tick's low-window runs, quiet and ones, at small primes that have them.

    python3.12 -I -B scripts/runs_check.py

TwoBitMod._runs finds no run at model_check's primes, so neither the quiet run's Prefix nor a ones
run's complemented Prefix is reached there. Each prime here has a ones run, a quiet run or both.
The dialog multiplier runs forward and inverse on every x, against every y below 300 and 6
sampled y above, and every output and scratch wire must be exact.

Controls, which must fail at some prime: the bridge without its CCX, the bridge without its
s term, and every ones run treated as a quiet run.
"""
import io
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.emit import write_circuit
from qec.point_add import DialogMultiplier
from qec.two_bit import TwoBitMod

PRIMES = [97, 113, 193, 223, 241, 449, 479, 769, 1567, 1777, 1789]


def bad_lanes(q, inverse):
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    X0, Y0 = list(X), list(Y)
    # apply relabels X and Y to the wires that hold the outputs.
    DialogMultiplier(C, mod).apply(X, Y, inverse=inverse)
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y}, depth_pass=False)
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    ys = range(q) if q < 300 else random.Random(q).sample(range(q), 6)
    xs, values = zip(*[(x, y) for x in range(1, q) for y in ys])
    state = sim.run(ops, C.num_wires, sim.pack(xs, X0) | sim.pack(values, Y0), len(xs))
    expected = [y * (pow(x, -1, q) if inverse else x) % q for x, y in zip(xs, values)]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def sweep(label):
    total = 0
    for q in PRIMES:
        mod = TwoBitMod(q)
        assert mod.merged and mod.runs, q
        for inverse in (False, True):
            bad = bad_lanes(q, inverse)
            total += bad
            if label == "build":
                print(f"q={q} runs={mod.runs} inverse={inverse}: {bad} bad lanes")
    return total


def main():
    failed = sweep("build") > 0
    real_bridge, real_runs = TwoBitMod._bridge, TwoBitMod._runs

    def no_ccx(self, C, selectors, masks, a, carry, top, out):
        active = [sel for sel, mask in zip(selectors, masks) if mask >> a & 1]
        for sel in active:
            C.cx(sel, out)

    def no_s(self, C, selectors, masks, a, carry, top, out):
        C.ccx(carry, top, out)

    def all_quiet(self):
        return [(a, b, False) for a, b, _ in real_runs(self)]

    caught = 0
    controls = (("bridge without its CCX", "_bridge", no_ccx),
                ("bridge without s", "_bridge", no_s),
                ("ones runs as quiet", "_runs", all_quiet))
    for label, name, fake in controls:
        setattr(TwoBitMod, name, fake)
        hit = sweep(label) > 0
        setattr(TwoBitMod, name, real_bridge if name == "_bridge" else real_runs)
        caught += hit
        print(f"control {label}: {'caught' if hit else 'NOT caught'}")
    print(f"runs_check: {'FAIL' if failed else 'PASS'}; controls caught {caught} of {len(controls)}")
    sys.exit(1 if failed or caught != len(controls) else 0)


if __name__ == "__main__":
    main()
