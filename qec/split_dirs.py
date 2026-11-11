"""Give definitions that run in both directions a separate reversed twin.

    python3 -I qec/split_dirs.py <in.json> <out.json> [default|all|<def,def,...>]

Every gate in the circuit is self-inverse (X, CX, CCX, SWAP), so a definition run
inverted is the same as its operations run in reverse order. For each chosen
definition d this adds a twin dR whose operations are d's in reverse order, with
each call's inverse flag flipped, and then rewrites every call (d, inverse) as a call
(dR, forward). After that, d and dR each run only forward, so later passes (control
orientation, gate reordering) can tune each for the direction it actually runs in
instead of a compromise between the two. default_targets says which definitions
are split by default.

The function of the circuit is unchanged: a call (dR, forward) executes exactly the
gate sequence that (d, inverse) executed. Each twin is placed right after its
original and the definitions are renumbered, so a definition still calls only
earlier ones.
"""
import json
import sys

SELF_INVERSE = {"X", "CX", "CCX", "SWAP"}


def split(circuit, targets):
    defs = circuit["definitions"]
    targets = set(targets)
    for d in targets:
        for r in defs[d]["operations"]:
            if "call" not in r and r["op"] not in SELF_INVERSE:
                raise SystemExit(f"definition {d} holds {r['op']}, which is not self-inverse")
    # New order: each twin sits right after its original, so callees stay earlier.
    pos, twin, order = {}, {}, []
    for d in range(len(defs)):
        pos[d] = len(order); order.append((d, False))
        if d in targets:
            twin[d] = len(order); order.append((d, True))

    def renumber(r):
        """Copy of record r with the call index in the new numbering; an inverted call to a
        split definition goes to its twin, run forward."""
        r = dict(r)
        if "call" in r:
            d = r["call"]
            if r["inverse"] and d in twin:
                r["call"] = twin[d]
                r["inverse"] = False
            else:
                r["call"] = pos[d]
        return r

    def flip(r):
        r = dict(r)
        if "call" in r:
            r["inverse"] = not r["inverse"]
        return r

    out = []
    for d, rev in order:
        new = dict(defs[d])
        ops = [flip(r) for r in reversed(defs[d]["operations"])] if rev else defs[d]["operations"]
        new["operations"] = [renumber(r) for r in ops]
        out.append(new)
    circuit["definitions"] = out
    circuit["operations"] = [renumber(r) for r in circuit["operations"]]
    return {d: (pos[d], twin[d]) for d in sorted(twin)}


def both_directions(circuit):
    """Definitions that some execution runs forward and another runs inverted."""
    defs = circuit["definitions"]
    seen = {}

    def walk(ops, inverse, memo):
        for r in ops:
            if "call" in r:
                d, inv = r["call"], inverse ^ bool(r.get("inverse", False))
                if (d, inv) in memo:
                    continue
                memo.add((d, inv))
                seen.setdefault(d, set()).add(inv)
                walk(defs[d]["operations"], inv, memo)

    walk(circuit["operations"], False, set())
    return sorted(d for d, s in seen.items() if len(s) == 2)


def run_counts(circuit):
    """Executions of each definition in the expanded circuit."""
    defs = circuit["definitions"]
    runs = [0] * len(defs)
    for r in circuit["operations"]:
        if "call" in r:
            runs[r["call"]] += 1
    for d in range(len(defs) - 1, -1, -1):
        for r in defs[d]["operations"]:
            if "call" in r:
                runs[r["call"]] += runs[d]
    return runs


def default_targets(circuit, min_runs=8):
    """Definitions that run in both directions and either hold calls (so their twin can
    route to the callees' twins) or run at least min_runs times. Splitting every
    definition would add the bytes of every body run both ways; this keeps the bodies
    that dominate the run count and leaves bodies that run only a few times shared."""
    runs = run_counts(circuit)
    defs = circuit["definitions"]
    return [d for d in both_directions(circuit)
            if runs[d] >= min_runs or any("call" in r for r in defs[d]["operations"])]


def main(src, dst, spec="default"):
    with open(src) as fh:
        doc = json.load(fh)
    circuit = doc.get("circuit", doc)
    if spec == "default":
        targets = default_targets(circuit)
    elif spec == "all":
        targets = both_directions(circuit)
    else:
        targets = [int(x) for x in spec.split(",") if x]
    twin = split(circuit, targets)
    with open(dst, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(json.dumps({"twins": {str(k): v for k, v in twin.items()}, "definitions": len(circuit["definitions"])}))


if __name__ == "__main__":
    main(*sys.argv[1:4])
