"""Circuit-level post-passes, in order: split_dirs, gate_pass, reorder_pass.

    python3 -I -B scripts/post_pass.py <emitted.json> <out.json> [--jobs N] [--rounds N]

Each pass keeps the circuit's function; see each module's docstring for why.
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "qec"))

import gate_pass
import reorder_pass
import split_dirs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--rounds", type=int, default=0, help="reorder_pass iterated local search rounds")
    a = ap.parse_args()
    with open(a.src) as fh:
        doc = json.load(fh)
    circuit = doc.get("circuit", doc)
    t0 = time.time()
    twins = split_dirs.split(circuit, split_dirs.default_targets(circuit))
    print(json.dumps({"split": len(twins), "definitions": len(circuit["definitions"])}), flush=True)
    print(json.dumps({"gate_pass": gate_pass.apply(circuit)}), flush=True)
    print(json.dumps({"reorder_pass": reorder_pass.apply(circuit, jobs=a.jobs, rounds=a.rounds)}), flush=True)
    with open(a.dst, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(json.dumps({"seconds": round(time.time() - t0)}))


if __name__ == "__main__":
    main()
