"""Ticks of a chain of Bezout ticks on secp256k1, alone, under the board's tick model.

    python3.12 -I -B scripts/tick_bench.py [N] [--crit]

Records the steady-state tick (consume and produce, a prev and a nxt) once and applies it N times
in a row as the pass does, on a flat circuit. Prints total ticks, ticks per extra tick, CCX per
tick, and with --crit the critical path of the last tick, by op index within the tick body.
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from qec.builder import Circuit
from qec.two_bit import TwoBitMod

Q = 2**256 - 2**32 - 977


def ccx_ticks(ra, rb, rt):
    rt += 1
    x = max(rb, rt) + 1; rb = rt = x
    rt += 1
    x = max(ra, rt) + 1; ra = rt = x
    rt += 1
    x = max(rb, rt) + 1; rb = rt = x
    rt += 1
    x = max(ra, rt) + 1; ra = rt = x
    rb += 1
    rt += 1
    rt += 1
    x = max(ra, rb) + 1; ra = rb = x
    ra += 1; rb += 1
    x = max(ra, rb) + 1; ra = rb = x
    return ra, rb, rt


def schedule(prims, nw):
    ready = [0] * nw
    last = [-1] * nw
    n = len(prims) // 4
    end = [0] * n
    pred = [-1] * n
    for g in range(n):
        k, a, b, c = prims[4 * g:4 * g + 4]
        if k == 0:
            ws = (a,)
        elif k == 2:
            ws = (a, b, c)
        else:
            ws = (a, b)
        bw = max(ws, key=lambda w: ready[w])
        pred[g] = last[bw]
        if k == 0:
            ready[a] += 1
            e = ready[a]
        elif k == 1:
            e = max(ready[a], ready[b]) + 1
            ready[a] = ready[b] = e
        elif k == 3:
            e = max(ready[a], ready[b]) + 3
            ready[a] = ready[b] = e
        else:
            ra, rb, rt = ccx_ticks(ready[a], ready[b], ready[c])
            ready[a], ready[b], ready[c] = ra, rb, rt
            e = max(ra, rb, rt)
        end[g] = e
        for w in ws:
            last[w] = g
    return end, pred


PHASES = []
ORIG = []
CMP_IO = {}


def _instrument(rec):
    """Wrap the tick's routines to record (name, first op, end op) within the tick body."""
    from qec import two_bit, wide

    def wrap(owner, attr, name):
        f = getattr(owner, attr)

        def g(*a, **k):
            s = len(rec.prims) // 4
            out = f(*a, **k)
            PHASES.append((name, s, len(rec.prims) // 4))
            return out
        setattr(owner, attr, g)

    wrap(two_bit, "add_or_shift_la", "add")
    orig_gt = wide.gt_tree

    def gt(C_, A, B, *a, **k):
        CMP_IO.update(A=list(A), B=list(B), pos=len(rec.prims) // 4)
        return orig_gt(C_, A, B, *a, **k)
    wide.gt_tree = gt
    wrap(wide, "gt_tree", "cmp")
    orig_w = wide.PrefixTree.write

    def pw(self, *a, **k):
        st = len(rec.prims) // 4
        out = orig_w(self, *a, **k)
        PHASES.append((f"write{len(self.bits)}", st, len(rec.prims) // 4))
        return out
    wide.PrefixTree.write = pw
    wrap(wide.LowAdd, "undo", "low1undo")
    wrap(wide.LowAdd, "write", "low2write")
    orig_init = wide.PrefixTree.__init__

    def pinit(self, C, bits, *a, **k):
        s = len(rec.prims) // 4
        orig_init(self, C, bits, *a, **k)
        PHASES.append((f"Ptree{len(bits)}", s, len(rec.prims) // 4))
    wide.PrefixTree.__init__ = pinit
    orig_low = wide.LowAdd.__init__

    def linit(self, *a, **k):
        s = len(rec.prims) // 4
        orig_low(self, *a, **k)
        PHASES.append(("lowinit", s, len(rec.prims) // 4))
    wide.LowAdd.__init__ = linit


def phase_of(op):
    best = None
    for name, s, e in PHASES:
        if s <= op < e and (best is None or e - s < best[2] - best[1]):
            best = (name, s, e)
    return best[0] if best else "glue"


def build(N):
    mod = TwoBitMod(Q)
    n = mod.n
    C = Circuit()
    rec = Circuit()
    if "--crit" in sys.argv or "--arrive" in sys.argv:
        _instrument(rec)

    def tick(prev, g, nxt, R, S, pin, pout):
        mod.tick(rec, prev or None, g, nxt[0] if nxt else None, R, S, pin=pin[0] if pin else None,
                 pout=pout[0] if pout else None)

    block = rec.record_shape(tick, [2, 2, 1, n, n, 1, 1])
    G = C.alloc(2 * N + 4)
    R, S = C.alloc(n), C.alloc(n)
    hand = C.alloc(1)
    bounds = []
    for i in reversed(range(1, N + 1)):
        pout = C.alloc(1)
        bounds.append(len(C.prims) // 4)
        C.apply(block, [G[2 * i + 2:2 * i + 4], G[2 * i:2 * i + 2], G[2 * i - 2:2 * i - 1], R, S,
                        hand, pout])
        C.free(hand)
        hand = pout
    bounds.append(len(C.prims) // 4)
    return C, block, bounds


def post_pass(C, block, bounds):
    """gate_pass's cancel and control order on the tick body as one definition, then the chain
    re-expanded flat."""
    from qec import gate_pass
    names = {0: "X", 1: "CX", 2: "CCX", 3: "SWAP"}
    local = sorted(block.wires)
    idx = {w: i for i, w in enumerate(local)}
    body = []
    ops = block.ops
    for g in range(len(ops) // 4):
        k, a, b, c = ops[4 * g:4 * g + 4]
        if k == 0:
            body.append({"op": "X", "target": idx[a], "orig": g})
        elif k == 1:
            body.append({"op": "CX", "control": idx[a], "target": idx[b], "orig": g})
        elif k == 2:
            body.append({"op": "CCX", "control1": idx[a], "control2": idx[b], "target": idx[c], "orig": g})
        else:
            body.append({"op": "SWAP", "wire1": idx[a], "wire2": idx[b], "orig": g})
    # each call's wire map: read off the expanded ops of that call
    calls = []
    for j in range(len(bounds) - 1):
        m = [None] * len(local)
        seg = C.prims[4 * bounds[j]:4 * bounds[j + 1]]
        for g in range(len(ops) // 4):
            for s in (1, 2, 3):
                m[idx[ops[4 * g + s]]] = seg[4 * g + s]
        calls.append({"call": 0, "wires": m})
    circ = {"num_wires": C.num_wires, "definitions": [{"operations": body, "arity": len(local)}],
            "operations": calls}
    circ["definitions"][0]["operations"], _ = gate_pass.cancel_pairs(body)
    gate_pass.order_controls(circ)
    out = Circuit()
    out.num_wires = C.num_wires
    nb = [0]
    for call in calls:
        m = call["wires"]
        for r in circ["definitions"][0]["operations"]:
            o = r["op"]
            if o == "X":
                out.x(m[r["target"]])
            elif o == "CX":
                out.cx(m[r["control"]], m[r["target"]])
            elif o == "CCX":
                out.prims.extend((2, m[r["control1"]], m[r["control2"]], m[r["target"]]))
            else:
                out.swap(m[r["wire1"]], m[r["wire2"]])
        nb.append(len(out.prims) // 4)
    nops = len(circ["definitions"][0]["operations"])
    ORIG[:] = [r["orig"] for r in circ["definitions"][0]["operations"]]
    return out, nb, nops


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    N = int(args[0]) if args else 3
    C, block, bounds = build(N)
    if "--gp" in sys.argv:
        C, bounds, _ = post_pass(C, block, bounds)
    end, pred = schedule(C.prims, C.num_wires)
    nops = bounds[1] - bounds[0]
    ccx = sum(1 for k in C.prims[4 * bounds[0]:4 * bounds[1]:4] if k == 2)
    total = max(end)
    # each tick's finish: max end over its ops
    fins = [max(end[bounds[j]:bounds[j + 1]]) for j in range(N)]
    print(f"N {N} wires {C.num_wires} ops/tick {nops} ccx/tick {ccx} total {total} "
          f"finishes {fins} steady {(fins[-1] - fins[0]) / max(N - 1, 1):.1f}")
    if "--arrive" in sys.argv:
        ops = block.ops
        m = {}
        seg = C.prims[4 * bounds[0]:4 * bounds[1]]
        for g in range(len(ops) // 4):
            for k in (1, 2, 3):
                m[ops[4 * g + k]] = seg[4 * g + k]
        # ready of every wire just before the compare's first op
        pos = bounds[0] + CMP_IO["pos"]
        end2, _ = schedule(C.prims[:4 * pos], C.num_wires)
        ready = [0] * C.num_wires
        for g in range(pos):
            for k in (1, 2, 3):
                w = C.prims[4 * g + k]
                ready[w] = max(ready[w], end2[g])
        print("cmp starts after op", pos, "max end so far", max(end2))
        print("A ready:", [ready[m[w]] for w in CMP_IO["A"]])
        print("B ready:", [ready[m[w]] for w in CMP_IO["B"]])
    if "--crit" in sys.argv:
        g = max(range(len(end)), key=end.__getitem__)
        path = []
        while g >= 0:
            path.append(g)
            g = pred[g]
        path.reverse()
        names = {0: "X", 1: "CX", 2: "CCX", 3: "SWAP"}
        prev_end = 0
        verbose = "--v" in sys.argv
        run = None
        for g in path:
            j = max(i for i in range(N) if bounds[i] <= g)
            k, a, b, c = C.prims[4 * g:4 * g + 4]
            ph = phase_of(ORIG[g - bounds[j]] if ORIG else g - bounds[j])
            if verbose:
                print(f"tick{j} op{g - bounds[j]:6d} {ph:10s} {names[k]:4s} {a:5d} {b:5d} {c:5d} end {end[g]:6d} +{end[g] - prev_end}")
            key = (j, ph)
            if run is None or run[0] != key:
                if run:
                    print(f"tick{run[0][0]} {run[0][1]:10s} ops {run[3]:6d}-{run[4]:6d} n {run[5]:4d} ticks {run[2] - run[1]:5d} end {run[2]}")
                run = [key, prev_end, end[g], g - bounds[j], g - bounds[j], 1]
            else:
                run[2] = end[g]; run[4] = g - bounds[j]; run[5] += 1
            prev_end = end[g]
        print(f"tick{run[0][0]} {run[0][1]:10s} ops {run[3]:6d}-{run[4]:6d} n {run[5]:4d} ticks {run[2] - run[1]:5d} end {run[2]}")


if __name__ == "__main__":
    main()
