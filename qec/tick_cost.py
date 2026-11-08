"""Costs of a circuit under the board's model, written independently of any entry's code.

    python3 -I qec/tick_cost.py <circuit-or-method.json>

T = 7 per CCX/CCZ. qubit_time = num_wires * ticks. Ticks: ASAP, one tick per lowered gate,
SWAP = 3 CX, CCX/CCZ lowered to the standard 15-gate Clifford+T sequence
(H t; CX b,t; Tdg t; CX a,t; T t; CX b,t; Tdg t; CX a,t; T b; T t; H t; CX a,b; T a; Tdg b; CX a,b).
"""
import json
import sys
from collections import defaultdict


def load(path):
    with open(path) as fh:
        doc = json.load(fh)
    return doc.get("circuit", doc)


def ccx_ticks(ra, rb, rt):
    # Same gate sequence as the docstring, one tick each, ASAP per wire.
    rt += 1                       # H t
    x = max(rb, rt) + 1; rb = rt = x   # CX b,t
    rt += 1                       # Tdg t
    x = max(ra, rt) + 1; ra = rt = x   # CX a,t
    rt += 1                       # T t
    x = max(rb, rt) + 1; rb = rt = x   # CX b,t
    rt += 1                       # Tdg t
    x = max(ra, rt) + 1; ra = rt = x   # CX a,t
    rb += 1                       # T b
    rt += 1                       # T t
    rt += 1                       # H t
    x = max(ra, rb) + 1; ra = rb = x   # CX a,b
    ra += 1; rb += 1              # T a, Tdg b
    x = max(ra, rb) + 1; ra = rb = x   # CX a,b
    return ra, rb, rt


class Coster:
    def __init__(self, circuit):
        self.c = circuit
        self.bodies = [d["operations"] for d in circuit.get("definitions", [])]
        self.ready = [0] * circuit["num_wires"]
        self.counts = defaultdict(int)
        self.calls = defaultdict(int)

    def run(self, ops, wmap, inverse, body_id=None):
        ready = self.ready
        counts = self.counts
        seq = reversed(ops) if inverse else ops
        for r in seq:
            if "call" in r:
                sub = [wmap[w] for w in r["wires"]] if wmap is not None else r["wires"]
                self.calls[r["call"]] += 1
                self.run(self.bodies[r["call"]], sub, inverse ^ bool(r.get("inverse", False)), r["call"])
                continue
            o = r["op"]
            counts[o] += 1
            m = wmap.__getitem__ if wmap is not None else (lambda w: w)
            if o in ("X", "Z"):
                ready[m(r["target"])] += 1
            elif o in ("CX", "CZ"):
                a, t = m(r["control"]), m(r["target"])
                ready[a] = ready[t] = max(ready[a], ready[t]) + 1
            elif o == "SWAP":
                a, b = m(r["wire1"]), m(r["wire2"])
                ready[a] = ready[b] = max(ready[a], ready[b]) + 3
            elif o in ("CCX", "CCZ"):
                a, b, t = m(r["control1"]), m(r["control2"]), m(r["target"])
                ready[a], ready[b], ready[t] = ccx_ticks(ready[a], ready[b], ready[t])
            elif o == "CONTROLLED_NEG":
                ready[m(r["control"])] += 1
            elif o == "NEG":
                pass
            else:
                raise SystemExit(f"unknown op {o}")

    def total(self):
        self.run(self.c["operations"], None, False)
        ticks = max(self.ready)
        ccx = self.counts["CCX"] + self.counts["CCZ"]
        return {"wires": self.c["num_wires"], "ticks": ticks, "T": 7 * ccx, "ccx": ccx,
                "qubit_time": self.c["num_wires"] * ticks, "ops": sum(self.counts.values()),
                "counts": dict(self.counts)}


if __name__ == "__main__":
    c = load(sys.argv[1])
    res = Coster(c).total()
    print(json.dumps(res))
