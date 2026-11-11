"""Check the circuit-level post-passes on toy curves, exhaustively.

    python3 -I scripts/check_passes.py [--toy-bits 4 6 7 8 9 10] [--toy-per-size 6]

For each toy circuit: emit the hierarchical JSON, run qec/split_dirs.py (every
definition that runs in both directions), qec/gate_pass.py and qec/reorder_pass.py
on it in memory, expand the result with qec/flatten.py's rules and check every lane
(c, R) as scripts/check.py does. Also reports the board ticks before and after.
"""
import argparse
import io
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "qec"))

import check as basis
import gate_pass
import reorder_pass
import split_dirs
from qec import flatten
from qec.curve import toy_secp_curves
from qec.emit import write_circuit
from qec.point_add import ControlledPointAdd


def passes(circuit):
    split_dirs.split(circuit, split_dirs.both_directions(circuit))
    gate_pass.apply(circuit)
    return reorder_pass.apply(circuit, jobs=1, log=False, rounds=3, rounds_min_runs=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy-bits", type=int, nargs="*", default=[4, 6, 7, 8, 9, 10])
    ap.add_argument("--toy-per-size", type=int, default=6)
    args = ap.parse_args()
    bad = 0
    for bits in args.toy_bits:
        for curve in toy_secp_curves(bits)[: args.toy_per_size]:
            pts = curve.points()
            for A in [P for P in pts[1:] if P[1] != 0][:2]:
                t0 = time.time()
                gen = ControlledPointAdd(curve, A)
                C = gen.build(hierarchical=True)
                buf = io.StringIO()
                write_circuit(buf, C, gen.layout)
                circuit = json.loads(buf.getvalue())
                before = reorder_pass.schedule(circuit)
                stats = passes(circuit)
                doc = flatten.load(io.StringIO(json.dumps(circuit)))
                ops = basis.expand(doc)
                cases = [(c, R) for c in (0, 1) for R in pts]
                fails = basis.check(gen, ops, cases, circuit["num_wires"])
                bad += len(fails)
                print(json.dumps({"q": curve.q, "addend": list(A), "lanes": len(cases), "failures": len(fails),
                                  "ticks_before": before, "ticks_after": stats["ticks_out"],
                                  "seconds": round(time.time() - t0, 1)}), flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
