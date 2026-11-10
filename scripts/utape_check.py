"""The zero argument for tape on U's vacated wires, tested at the moment it is used.

The pass hands U[j] to the tape's supply right after the step at which the width drops to j or
below. The sentence tested is L6's consequence as the layout uses it: at that moment U[j] is 0,
and no gate touches it until a step takes it as a cell. The first part reads the wire on every
lane; the second scans the gates.

Controls: at each lending moment the same test on U[1 + k], for the k-th lent wire, which is
inside every window, must find a nonzero lane; and the zero caller, V = 0, keeps U = q, so a lent
wire must be nonzero there. The wire just below a lent one is no control: L6 is a worst-case bound
and real u is usually far smaller, so that wire is often zero on every lane too.
Usage: utape_check.py
"""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod, load_const
from qec.builder import Circuit
from qec.point_add import DialogMultiplier

SECP = 2**256 - 2**32 - 977


def pack(values, wires):
    return {w: int("".join("1" if v >> i & 1 else "0" for v in reversed(values)) or "0", 2) for i, w in enumerate(wires)}


def build(q):
    """The pass with the op index after every step, and the cells each step wrote."""
    C = Circuit()
    m = DialogMultiplier(C, Mod(q))
    U, V = C.alloc(m.mod.n), C.alloc(m.mod.n)
    U0, V0 = list(U), list(V)
    pool = V + C.alloc(m.tape_external)
    load_const(C, U, q)
    steps, apply = [], C.apply

    def record(block, regs, inverse=False):
        apply(block, regs, inverse)
        if len(regs) == 4 and len(regs[3]) == 3:
            steps.append((len(C.ops), list(regs[3][:2])))

    C.apply = record
    m._euclid_pass(U, pool, C.alloc(m.iters))
    C.apply = apply
    return C, m, U0, V0, steps


def case(q, xs):
    C, m, U, V, steps = build(q)
    n, T, lent = m.mod.n, m.iters, U[m.mod.n - m.u_lent:]
    state, start = sim.pack([0] * len(xs), []) | pack(xs, V), 0
    lanes = len(xs)
    nonzero = low_live = untouched_bad = spans = u0_touched = 0
    lent_seen = []
    for i, (end, cell) in enumerate(steps):
        state = dict(enumerate(sim.run(C.ops[start:end], C.num_wires, state, lanes)))
        start = end
        w = m.widths[i]
        nw = m.widths[i + 1] if i + 1 < T else m.pass_tail_width
        for j in range(nw, w):
            if U[j] in lent:
                lent_seen.append((U[j], end))
                nonzero += state[U[j]] != 0
                low_live += state[U[1 + len(lent_seen) - 1]] != 0
    # No gate touches a lent wire between its lending and the step that takes it as a cell.
    ops = C.ops
    first_cell = {}
    for end, cell in steps:
        for c in cell:
            first_cell.setdefault(c, end)
    for wire, at in lent_seen:
        use_end = first_cell.get(wire)
        step_starts = [e for e, _ in steps]
        k = step_starts.index(use_end)
        use_start = steps[k - 1][0] if k else 0
        spans += use_start > at
        for p in range(at, use_start, 4):
            if wire in (ops[p + 1], ops[p + 2], ops[p + 3]):
                untouched_bad += 1
                break
        # Control for the scan: U[0] is in every step, so the same interval must show it.
        u0_touched += any(U[0] in (ops[p + 1], ops[p + 2], ops[p + 3]) for p in range(at, use_start, 4))
    return {"q": q if q < SECP else "secp256k1", "lanes": lanes, "u_lent": m.u_lent, "lent wires seen": len(lent_seen),
            "lent wire nonzero at lending": nonzero, "gates on a lent wire before its cell step": untouched_bad,
            "control, low wire nonzero at lending": low_live, "intervals with gates": spans, "control, intervals touching U[0]": u0_touched}


def main():
    r = random.Random(20261102)
    rows = [case(q, list(range(1, q))) for q in (127, 251, 4093, 65521, 65519, 35533)]
    xs = [r.randrange(1, SECP) for _ in range(500)] + [1, 2, 3, SECP - 1, SECP - 2]
    rows.append(case(SECP, xs))
    zero = case(65521, [0])
    rows.append({"control, zero caller": zero})
    print(json.dumps(rows, indent=1))
    ok = all(r["lent wire nonzero at lending"] == 0 and r["gates on a lent wire before its cell step"] == 0
             and r["lent wires seen"] == r["u_lent"] for r in rows[:-1])
    controls = all(r["control, low wire nonzero at lending"] > 0 and r["control, intervals touching U[0]"] == r["intervals with gates"] > 0
                   for r in rows[:-1]) and zero["lent wire nonzero at lending"] > 0
    return 0 if ok and controls else 1


if __name__ == "__main__":
    sys.exit(main())
