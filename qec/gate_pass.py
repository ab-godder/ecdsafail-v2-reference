"""Gate-level post-pass for a hierarchical circuit: same unitary, fewer ticks.

    python3 -I qec/gate_pass.py <in.json> <out.json>

<in.json> is a circuit, or an upload whose circuit sits under "circuit"; the
output keeps that shape and is written compact. Two rewrites, each applied to
every definition body and to the top level, so every call site of a body gets
the same change:

1. Cancel pairs. Two identical self-inverse gates (X, CX, CCX, SWAP) cancel when
   every gate between them that shares a wire commutes with them. Gates of the
   form "X on a target under controls" commute when neither's target is among
   the other's controls. A call, or a SWAP, commutes with nothing that shares
   a wire, so it ends the search.
2. Order CCX controls. CCX is symmetric in its controls, so swapping them is the
   same gate. The board lowers CCX to a Clifford+T sequence in which control1 is
   needed 2 ticks later than control2 (finish = max(r_c1 + 7, r_c2 + 9, r_t + 10)
   in ticks, where r_w is when wire w is free). The pass schedules the expanded
   circuit, and for each CCX record it puts the control that is ready later in
   control1, by majority over that record's executions. It repeats while the
   tick count falls.

Every rewrite keeps the circuit's function on every input: (1) removes U U = I
across gates that commute with U, and (2) relabels a symmetric gate.
"""
import json
import sys


def wires_of(r):
    if "call" in r:
        return r["wires"]
    return [r[k] for k in ("control1", "control2", "control", "target", "wire1", "wire2") if k in r]


def key_of(r):
    o = r["op"]
    if o == "X":
        return ("X", r["target"])
    if o == "CX":
        return ("CX", r["control"], r["target"])
    if o == "CCX":
        return ("CCX", frozenset((r["control1"], r["control2"])), r["target"])
    if o == "SWAP":
        return ("SWAP", frozenset((r["wire1"], r["wire2"])))
    return None


def x_form(r):
    """(target, controls) for an X under controls, else None."""
    o = r.get("op")
    if o == "X":
        return r["target"], ()
    if o == "CX":
        return r["target"], (r["control"],)
    if o == "CCX":
        return r["target"], (r["control1"], r["control2"])
    return None


def commute(p, q):
    a, b = x_form(p), x_form(q)
    if a is None or b is None:
        return False
    return a[0] not in b[1] and b[0] not in a[1]


def cancel_pairs(ops):
    alive = [True] * len(ops)
    on_wire = {}
    removed = 0
    for i, r in enumerate(ops):
        ws = wires_of(r)
        k = None if "call" in r else key_of(r)
        partner = None
        if k is not None:
            for w in ws:
                found = None
                for j in reversed(on_wire.get(w, ())):
                    if not alive[j]:
                        continue
                    h = ops[j]
                    if "call" not in h and key_of(h) == k:
                        found = j
                        break
                    if not commute(h, r):
                        break
                if found is None or (partner is not None and partner != found):
                    partner = None
                    break
                partner = found
        if partner is not None:
            alive[partner] = alive[i] = False
            removed += 2
            continue
        for w in ws:
            on_wire.setdefault(w, []).append(i)
    return [r for r, a in zip(ops, alive) if a], removed


def schedule(circuit, vote=None):
    """Board tick model. With vote, tally per CCX record whether control2 is ready later."""
    bodies = [d["operations"] for d in circuit.get("definitions", [])]
    ready = [0] * circuit["num_wires"]

    def run(body_id, ops, wmap, inverse):
        for idx in (range(len(ops) - 1, -1, -1) if inverse else range(len(ops))):
            r = ops[idx]
            if "call" in r:
                sub = [wmap[w] for w in r["wires"]] if wmap is not None else r["wires"]
                run(r["call"], bodies[r["call"]], sub, inverse ^ bool(r.get("inverse", False)))
                continue
            o = r["op"]
            m = (lambda w: wmap[w]) if wmap is not None else (lambda w: w)
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
                ra, rb, rt = ready[a], ready[b], ready[t]
                if vote is not None and ra != rb:
                    tally = vote.setdefault((body_id, idx), [0, 0])
                    tally[rb > ra] += 1
                rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
                x = max(ra, rt) + 1; ra = rt = x; rt += 1
                x = max(rb, rt) + 1; rb = rt = x; rt += 1
                x = max(ra, rt) + 1; ra = rt = x
                rb += 1; rt += 2
                x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
                x = max(ra, rb) + 1; ra = rb = x
                ready[a], ready[b], ready[t] = ra, rb, rt
            elif o == "CONTROLLED_NEG":
                ready[m(r["control"])] += 1
            elif o != "NEG":
                raise SystemExit(f"unknown op {o}")

    run(None, circuit["operations"], None, False)
    return max(ready)


def order_controls(circuit, rounds=4):
    bodies = [d["operations"] for d in circuit.get("definitions", [])]
    best = schedule(circuit)
    for _ in range(rounds):
        vote = {}
        schedule(circuit, vote)
        flips = [(b, i) for (b, i), (keep, flip) in vote.items() if flip > keep]
        if not flips:
            break
        for b, i in flips:
            r = (bodies[b] if b is not None else circuit["operations"])[i]
            r["control1"], r["control2"] = r["control2"], r["control1"]
        ticks = schedule(circuit)
        if ticks >= best:
            for b, i in flips:
                r = (bodies[b] if b is not None else circuit["operations"])[i]
                r["control1"], r["control2"] = r["control2"], r["control1"]
            break
        best = ticks
    return best


def main(src, dst):
    with open(src) as fh:
        doc = json.load(fh)
    circuit = doc.get("circuit", doc)
    before = schedule(circuit)
    removed = 0
    for d in circuit.get("definitions", []):
        d["operations"], n = cancel_pairs(d["operations"])
        removed += n
    circuit["operations"], n = cancel_pairs(circuit["operations"])
    removed += n
    after_cancel = schedule(circuit)
    after = order_controls(circuit)
    with open(dst, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(json.dumps({"ticks_in": before, "records_cancelled": removed, "ticks_after_cancel": after_cancel,
                      "ticks_out": after, "wires": circuit["num_wires"]}))


if __name__ == "__main__":
    main(*sys.argv[1:3])
