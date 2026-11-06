"""Per-component op and CCX counts for the secp256k1 instance."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.builder import Circuit
from qec.curve import SECP256K1, SECP256K1_G
from qec.point_add import ControlledPointAdd, DialogMultiplier


def count(fn, sizes):
    C = Circuit()
    regs = [C.alloc(k) for k in sizes]
    fn(C, *regs)
    k = C.ops[0::4]
    return {"ops": len(C), "ccx": k.count(2), "cx": k.count(1), "x": k.count(0), "wires": C.num_wires}


gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
n, mod = gen.n, gen.mod
out = {}


def euclid_only(C, U, V):
    M = DialogMultiplier(C, mod)
    G = C.alloc(2 * M.iters)
    for i in range(M.iters):
        C.apply(M._euclid_step(M.widths[i]), [U, V, G[2 * i:2 * i + 2]])


def bezout_only(C, Y, S, G):
    M = DialogMultiplier(C, mod)
    for i in reversed(range(M.iters)):
        C.apply(M._bezout(), [G[2 * i:2 * i + 2], Y, S])


def mulsub(C, g, X, Y):
    gen.C, gen._modadd_blk = C, None
    gen._mulsub(g[0], Y, X)


out["euclid_pass"] = count(euclid_only, [n + 2, n + 2])
out["bezout_pass"] = count(bezout_only, [n, n, 2 * 512])
out["dialog_multiplier_M"] = count(lambda C, X, Y: DialogMultiplier(C, mod).apply(X, Y), [n, n])
out["lambda_square_subtract"] = count(mulsub, [1, n, n])
full = gen.build()
k = full.ops[0::4]
out["controlled_point_add"] = {"ops": len(full), "ccx": k.count(2), "cx": k.count(1), "x": k.count(0), "wires": full.num_wires}
print(json.dumps(out, indent=1))
