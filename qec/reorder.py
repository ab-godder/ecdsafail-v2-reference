"""Reorder commuting gates inside each definition body to cut board ticks. Function-preserving.

    python3 -I qec/reorder.py <in.json> <out.json> [--defs 1,2,...] [--orient] [--passes N]

Two X-type gates (X/CX/CCX: target flipped under controls) commute when neither's target is
among the other's controls. On one wire, gates that use it only as a control commute with
each other, and gates that use it only as a target commute with each other, so a body's
legal orders are the topological orders of runs on each wire. SWAP, Z-type gates and calls
are barriers on their wires. A list scheduler under the board's tick model picks an order,
using as input ready times the profile seen at the body's first execution (a body first run
inverted is scheduled reversed). Every flat body with at least --min-ops records is rewritten.
With --orient, each CCX placed also gets the control order that ends its target soonest; CCX
is symmetric in its controls, so this is the same gate.
"""
import heapq
ORIENT = False
import json
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from tick_cost import Coster, ccx_ticks, load


def roles(r):
    """(target, controls) for X-type, else None (barrier) with its wires."""
    o = r.get("op")
    if o == "X":
        return r["target"], ()
    if o == "CX":
        return r["target"], (r["control"],)
    if o == "CCX":
        return r["target"], (r["control1"], r["control2"])
    return None


def wires_of(r):
    if "call" in r:
        return list(r["wires"])
    return [r[k] for k in ("control1", "control2", "control", "target", "wire1", "wire2") if k in r]


def finish(r, ready):
    """New ready times of r's wires if started now (dict wire -> time)."""
    if "call" in r:
        return None
    o = r["op"]
    if o in ("X", "Z"):
        t = r["target"]
        return {t: ready[t] + 1}
    if o in ("CX", "CZ"):
        a, t = r["control"], r["target"]
        x = max(ready[a], ready[t]) + 1
        return {a: x, t: x}
    if o == "SWAP":
        a, t = r["wire1"], r["wire2"]
        x = max(ready[a], ready[t]) + 3
        return {a: x, t: x}
    if o in ("CCX", "CCZ"):
        a, b, t = r["control1"], r["control2"], r["target"]
        ra, rb, rt = ccx_ticks(ready[a], ready[b], ready[t])
        return {a: ra, b: rb, t: rt}
    if o == "CONTROLLED_NEG":
        a = r["control"]
        return {a: ready[a] + 1}
    return {}


def build_dag(ops):
    """preds count and successor lists; virtual barrier nodes keep it linear."""
    n = len(ops)
    succ = [[] for _ in range(n)]
    npred = [0] * n
    # per wire: (run_type, members of current run, barrier node of previous run or None)
    state = {}
    virt = []  # virtual nodes: index n + k

    def add_edge(u, v):
        if u >= n:
            vsucc[u - n].append(v)
        else:
            succ[u].append(v)
        npred[v] += 1

    vsucc = []
    for i, r in enumerate(ops):
        ro = None if "call" in r else roles(r)
        if ro is None:
            uses = [(w, "B") for w in wires_of(r)]
        else:
            t, cs = ro
            uses = [(t, "T")] + [(c, "C") for c in cs]
        for w, typ in uses:
            st = state.get(w)
            if st is None:
                state[w] = [typ, [i], None]
                continue
            rtyp, members, prev = st
            if typ == rtyp and typ != "B":
                members.append(i)
                if prev is not None:
                    add_edge(prev, i)
            else:
                # close current run: a virtual node after all members
                if len(members) == 1:
                    bar = members[0]
                else:
                    bar = n + len(vsucc)
                    vsucc.append([])
                    npred.append(0)
                    for m in members:
                        succ[m].append(bar)
                        npred[bar] += 1
                add_edge(bar, i)
                state[w] = [typ, [i], bar]
    return succ, vsucc, npred


def schedule_body(ops, ready_in):
    """List-schedule ops. ready_in: dict wire -> ready time at body start."""
    n = len(ops)
    succ, vsucc, npred = build_dag(ops)
    total = n + len(vsucc)
    # tail length (unit gate estimate: CCX 10, others 1) for priority
    dur = [10 if ops[i].get("op") in ("CCX", "CCZ") else (3 if ops[i].get("op") == "SWAP" else 1) for i in range(n)]
    tail = [0] * total
    # reverse topological via Kahn on reversed graph: compute order first
    order, indeg = [], list(npred)
    stack = [i for i in range(total) if indeg[i] == 0]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in (succ[u] if u < n else vsucc[u - n]):
            indeg[v] -= 1
            if indeg[v] == 0:
                stack.append(v)
    assert len(order) == total
    for u in reversed(order):
        best = 0
        for v in (succ[u] if u < n else vsucc[u - n]):
            if tail[v] > best:
                best = tail[v]
        tail[u] = best + (dur[u] if u < n else 0)
    ready = dict(ready_in)
    for r in ops:
        for w in wires_of(r):
            ready.setdefault(w, 0)
    indeg = list(npred)
    heap = []
    out = []

    def est(i):
        r = ops[i]
        ws = wires_of(r)
        return max(ready[w] for w in ws)

    def release(u):
        for v in (succ[u] if u < n else vsucc[u - n]):
            indeg[v] -= 1
            if indeg[v] == 0:
                push(v)

    def push(v):
        if v >= n:
            release(v)
        else:
            heapq.heappush(heap, (est(v), -tail[v], v))

    for i in range(total):
        if npred[i] == 0:
            push(i)
    while heap:
        s, nt, i = heapq.heappop(heap)
        e = est(i)
        if e > s:
            heapq.heappush(heap, (e, nt, i))
            continue
        r = ops[i]
        if ORIENT and r.get("op") == "CCX":
            a, b, t = r["control1"], r["control2"], r["target"]
            k1 = ccx_ticks(ready[a], ready[b], ready[t])
            k2 = ccx_ticks(ready[b], ready[a], ready[t])
            if (k2[2], max(k2)) < (k1[2], max(k1)):
                r = dict(r, control1=b, control2=a)
        f = finish(r, ready)
        if f is None:  # never reached: bodies with calls are skipped
            m = max(ready[w] for w in r["wires"]) + 1
            f = {w: m for w in r["wires"]}
        ready.update(f)
        out.append(r)
        release(i)
    assert len(out) == n, (len(out), n)
    return out


def body_profiles(circuit):
    """Ready-time vector (local wires) at each definition's first execution."""
    prof = {}
    c = Coster(circuit)
    bodies = c.bodies
    orig_run = c.run

    def run(ops, wmap, inverse, body_id=None):
        if body_id is not None and body_id not in prof and wmap is not None:
            base = min(c.ready[w] for w in wmap) if wmap else 0
            prof[body_id] = ({k: c.ready[w] - base for k, w in enumerate(wmap)}, inverse)
        return orig_run(ops, wmap, inverse, body_id)

    c.run = run
    c.total()
    return prof


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--defs", default="")
    ap.add_argument("--min-ops", type=int, default=30)
    ap.add_argument("--orient", action="store_true", help="also pick each CCX's control order")
    ap.add_argument("--passes", type=int, default=1)
    a = ap.parse_args()
    global ORIENT
    ORIENT = a.orient
    with open(a.src) as fh:
        doc = json.load(fh)
    c = doc.get("circuit", doc)
    base = Coster(c).total()["ticks"]
    print("ticks in", base, flush=True)
    defs = c["definitions"]
    for _ in range(a.passes):
        one_pass(c, defs, a)
        print("ticks after pass", Coster(c).total()["ticks"], flush=True)
    after = Coster(c).total()["ticks"]
    print("ticks out", after, f"{100*(after-base)/base:+.2f}%")
    with open(a.dst, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))


def one_pass(c, defs, a):
    prof = body_profiles(c)
    targets = [int(x) for x in a.defs.split(",") if x] or [i for i in range(len(defs)) if len(defs[i]["operations"]) >= a.min_ops]
    for d in targets:
        ops = defs[d]["operations"]
        if any("call" in r for r in ops):
            continue
        rin, inverse = prof.get(d, ({}, False))
        if inverse:
            # schedule the reversed body, then reverse back
            new = list(reversed(schedule_body(list(reversed(ops)), rin)))
        else:
            new = schedule_body(ops, rin)
        defs[d]["operations"] = new


if __name__ == "__main__":
    main()
