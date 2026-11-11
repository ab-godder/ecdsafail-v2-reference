"""Expand each hierarchical circuit with qec/flatten.py and compare the SHA-256 of
its flat `operations` array with receipts/baseline/flat_ops.sha256 (toy-127,
secp256k1) or with a freshly emitted flat stream (--all-toys)."""
import argparse
import hashlib
import io
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec import flatten
from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves
from qec.emit import _op, write_circuit
from qec.point_add import ControlledPointAdd

TOY_BITS, TOY_PER_SIZE = (4, 6, 7, 8, 9, 10), 6
PINNED = ("toy127", "secp256k1")


def baseline():
    lines = (ROOT / "receipts/baseline/flat_ops.sha256").read_text().split("\n")
    return {name: digest for digest, name, _ in (line.split() for line in lines if line)}


def flat_digest(gen):
    ops = gen.build().ops
    it = iter(ops)
    return hashlib.sha256(b"[" + ",".join(_op(*t) for t in zip(it, it, it, it)).encode() + b"]").hexdigest()


def hier_text(gen):
    buf = io.StringIO()
    # The pin is the builder's flat order; the depth pass reorders it, so it is checked without.
    write_circuit(buf, gen.build(hierarchical=True), gen.layout, depth_pass=False)
    return buf.getvalue()


def hier_digest(text):
    doc = flatten.load(io.StringIO(text))
    h = hashlib.sha256()
    for chunk in flatten.operations_json(doc):
        h.update(chunk.encode())
    instrs = [i for d in doc["definitions"] for i in d["operations"]] + doc["operations"]
    calls = [i for i in instrs if type(i) is flatten.Call]
    shape = {"definitions": len(doc["definitions"]), "call_sites": len(calls), "primitive_ops": len(instrs) - len(calls),
             "call_wire_args": sum(len(c.wires) for c in calls), "bytes": len(text.encode())}
    return h.hexdigest(), shape


def instances(all_toys):
    toy = next(c for c in toy_secp_curves(7) if c.q == 127)
    yield "toy127", toy, toy.points()[5]
    yield "secp256k1", SECP256K1, SECP256K1_G
    if all_toys:
        for bits in TOY_BITS:
            for curve in toy_secp_curves(bits)[:TOY_PER_SIZE]:
                for A in [P for P in curve.points()[1:] if P[1] != 0][:2]:
                    yield f"toy q={curve.q} A={A}", curve, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-toys", action="store_true", help="also every toy curve and addend that scripts/check.py runs")
    ap.add_argument("--only", choices=PINNED)
    ap.add_argument("--file", type=Path, help="with --only, check this hierarchical JSON instead of emitting one")
    args = ap.parse_args()
    if args.file and not args.only:
        ap.error("--file needs --only")
    pinned = baseline()
    bad = 0
    for name, curve, A in instances(args.all_toys):
        if args.only and name != args.only:
            continue
        want = None
        try:
            gen = ControlledPointAdd(curve, A)
            text = args.file.read_text() if args.file else hier_text(gen)
            want = pinned[name] if name in PINNED else flat_digest(gen)
            got, shape = hier_digest(text)
            ok = got == want
        except Exception as e:
            traceback.print_exc()
            ok, got, shape = False, None, {"error": repr(e)[:200]}
        bad += not ok
        print(f"{name}: hier {'==' if ok else '!='} flat {json.dumps({'digest': got, 'want': want, **shape})}", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
