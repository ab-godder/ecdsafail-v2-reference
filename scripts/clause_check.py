"""The clause as a statement about wires, on every state of two consecutive ticks.

Claim: for any U, V, E, a tick that leaves its cell at g0 = 0, b = 0 is followed by a tick
whose cell has g0 = 1. The states here are all of them, reachable or not. The control asks
the same of a cell (0, 1), which must fail.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier

Q = 251


def pair(w, next_w):
    C = Circuit()
    mod = Mod(Q)
    m = DialogMultiplier(C, mod)
    U, V, E = C.alloc(mod.n), C.alloc(w), C.alloc(1)
    a, b = C.alloc(3), C.alloc(3)
    inputs = U[:w] + list(V) + E
    C.apply(m._euclid_step(w), [U[:w], V, E, a])
    Vn = V[:next_w]
    C.apply(m._euclid_step(next_w), [U[:next_w], Vn, E, b])
    lanes = 1 << len(inputs)
    init = sim.pack(list(range(lanes)), inputs)
    state = sim.run(C.ops, C.num_wires, init, lanes)
    g, t, h = sim.unpack(state, [a[0]], lanes), sim.unpack(state, [a[1]], lanes), sim.unpack(state, [b[0]], lanes)
    seen = sum(1 for i in range(lanes) if not g[i] and not t[i])
    bad = sum(1 for i in range(lanes) if not g[i] and not t[i] and not h[i])
    other = sum(1 for i in range(lanes) if not g[i] and t[i] and not h[i])
    return {"w": w, "next": next_w, "states": lanes, "cell_00": seen, "next_g0_zero": bad, "control_01_next_g0_zero": other}


rows = [pair(w, nw) for w in range(2, 9) for nw in (w, w - 1) if nw >= 2]
print(json.dumps(rows))
sys.exit(0 if all(r["next_g0_zero"] == 0 and r["cell_00"] and r["control_01_next_g0_zero"] for r in rows) else 1)
