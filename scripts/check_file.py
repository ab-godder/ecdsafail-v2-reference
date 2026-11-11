"""Simulate a secp256k1 circuit file at A = G, as scripts/check.py --secp does for an emission.

    python3.12 -I -B scripts/check_file.py <circuit-or-method.json> [--secp-random 120] [--seed 1]

For a post-passed circuit: the file is expanded by qec/flatten.py's rules and every lane (c, R),
c in {0, 1}, R in inf, A, -A, -2A, 2A, 3A and random multiples of G, must give R + c*A with the
control kept and every other wire zero.
"""
import argparse
import io
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from check import check, expand, special_points, structural
from qec import flatten
from qec.curve import SECP256K1, SECP256K1_G
from qec.point_add import ControlledPointAdd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", type=Path)
    ap.add_argument("--secp-random", type=int, default=120)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    doc = json.loads(a.file.read_text())
    circuit = doc.get("circuit", doc)
    parsed = flatten.load(io.StringIO(json.dumps(circuit)))
    ops = expand(parsed)
    counts = structural(ops, circuit["num_wires"])
    gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
    gen.build(hierarchical=True)
    rng = random.Random(a.seed)
    curve, A = SECP256K1, SECP256K1_G
    rand = [curve.mul(rng.randrange(1, 2 ** 256), A) for _ in range(a.secp_random)]
    cases = [(c, R) for c in (0, 1) for R in special_points(curve, A) + rand]
    fails = check(gen, ops, cases, circuit["num_wires"])
    print(json.dumps({"file": a.file.name, "wires": circuit["num_wires"], "ops": len(ops) // 4, "counts": counts,
                      "lanes": len(cases), "failures": len(fails), "first_failures": fails[:3]}))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
