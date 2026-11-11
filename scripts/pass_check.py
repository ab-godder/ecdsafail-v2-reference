"""Certify the depth pass: each emitted definition is its builder body with commuting ops reordered
and CCX controls possibly exchanged, and the top level differs only in CCX control order.

    python3.12 -I -B scripts/pass_check.py [--toy-q 127] [--no-secp]

Shares no code with qec/depth_pass.py. Both texts come from qec/emit.py, with and without the pass.
An op acts on each of its wires as a reader (c: a control), a target (t: XORed with a function of
wires it does not write) or both (e: SWAP wires). A call's role on a wire is its callee's role over
the whole body, e when the body uses the wire in two roles. Two ops commute when every wire they
share is c/c or t/t, so on each wire the passed body must keep the run of ops between role changes
as the same set, in the same order; every non-commuting pair then keeps its order. The two bodies
must also be the same multiset of ops, a CCX's controls taken as a set.

Control: in the passed secp256k1 (or toy) body with the most ops, the first adjacent pair of ops
that share a wire in conflicting roles is exchanged, and the check must reject it.
"""
import argparse
import io
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves
from qec.emit import write_circuit
from qec.point_add import ControlledPointAdd


def texts(curve, A):
    gen = ControlledPointAdd(curve, A)
    C = gen.build(hierarchical=True)
    out = []
    for dp in (False, True):
        buf = io.StringIO()
        write_circuit(buf, C, gen.layout, depth_pass=dp)
        out.append(json.loads(buf.getvalue()))
    return out


def key(r):
    if "call" in r:
        return ("call", r["call"], tuple(r["wires"]), bool(r.get("inverse", False)))
    o = r["op"]
    if o == "CCX":
        return ("CCX", frozenset((r["control1"], r["control2"])), r["target"])
    if o == "CX":
        return ("CX", r["control"], r["target"])
    if o == "X":
        return ("X", r["target"])
    if o == "SWAP":
        return ("SWAP", frozenset((r["wire1"], r["wire2"])))
    raise ValueError(o)


def roles(r, summaries):
    if "call" in r:
        s = summaries[r["call"]]
        return [(w, s[k]) for k, w in enumerate(r["wires"]) if k in s]
    o = r["op"]
    if o == "CCX":
        return [(r["control1"], "c"), (r["control2"], "c"), (r["target"], "t")]
    if o == "CX":
        return [(r["control"], "c"), (r["target"], "t")]
    if o == "X":
        return [(r["target"], "t")]
    return [(r["wire1"], "e"), (r["wire2"], "e")]


def summarize(defs):
    out = []
    for d in defs:
        acc = {}
        for r in d["operations"]:
            for w, role in roles(r, out):
                acc[w] = role if acc.get(w, role) == role else "e"
        out.append(acc)
    return out


def runs(ops, ids, summaries):
    """Per wire, its list of (role, frozenset of op ids) runs, in the order ops gives."""
    seq = defaultdict(list)
    for r, i in zip(ops, ids):
        for w, role in roles(r, summaries):
            lst = seq[w]
            if lst and role != "e" and lst[-1][0] == role:
                lst[-1][1].add(i)
            else:
                lst.append((role, {i}))
    return {w: [(role, frozenset(s)) for role, s in lst] for w, lst in seq.items()}


def certify(before, after, summaries):
    """None if after is a commuting reorder of before (CCX controls as a set), else a reason."""
    if len(before) != len(after):
        return "op counts differ"
    slots = defaultdict(list)
    for i, r in enumerate(before):
        slots[key(r)].append(i)
    for k in slots:
        slots[k].reverse()
    ids = []
    for r in after:
        s = slots.get(key(r))
        if not s:
            return f"op {r} not in the builder body"
        ids.append(s.pop())
    if runs(before, range(len(before)), summaries) != runs(after, ids, summaries):
        return "a wire's role runs changed"
    return None


def conflict_swap(ops, summaries):
    """ops with the first adjacent conflicting pair exchanged, or None."""
    for i in range(len(ops) - 1):
        a, b = dict(roles(ops[i], summaries)), dict(roles(ops[i + 1], summaries))
        if any(w in b and not (a[w] == b[w] != "e") for w in a):
            return ops[:i] + [ops[i + 1], ops[i]] + ops[i + 2:]
    return None


def check(name, curve, A):
    plain, passed = texts(curve, A)
    assert len(plain["definitions"]) == len(passed["definitions"])
    summaries = summarize(plain["definitions"])
    bad, moved = [], 0
    for d, (x, y) in enumerate(zip(plain["definitions"], passed["definitions"])):
        if x["arity"] != y["arity"]:
            bad.append((d, "arity differs"))
            continue
        why = certify(x["operations"], y["operations"], summaries)
        if why:
            bad.append((d, why))
        moved += x["operations"] != y["operations"]
    if [key(r) for r in plain["operations"]] != [key(r) for r in passed["operations"]]:
        bad.append(("top", "top level differs beyond CCX control order"))
    big = max(range(len(passed["definitions"])), key=lambda d: len(passed["definitions"][d]["operations"]))
    mutant = conflict_swap(passed["definitions"][big]["operations"], summaries)
    control = mutant is not None and certify(plain["definitions"][big]["operations"], mutant, summaries) is not None
    row = {"instance": name, "definitions": len(plain["definitions"]), "changed": moved,
           "rejected": len(bad), "first": bad[:3], "control_rejected": control}
    print(json.dumps(row, default=str), flush=True)
    return not bad and control


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy-q", type=int, default=127)
    ap.add_argument("--no-secp", action="store_true")
    a = ap.parse_args()
    toy = next(c for c in toy_secp_curves(a.toy_q.bit_length()) if c.q == a.toy_q)
    ok = check(f"toy{a.toy_q}", toy, toy.points()[5])
    if not a.no_secp:
        ok &= check("secp256k1", SECP256K1, SECP256K1_G)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
