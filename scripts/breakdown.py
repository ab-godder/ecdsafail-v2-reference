"""Per-component op and CCX counts for the secp256k1 instance."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.builder import Circuit
from qec.curve import SECP256K1, SECP256K1_G
from qec.model import tick_bound
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


def euclid_only(C, X, Y):
    M = DialogMultiplier(C, mod)
    from qec.arith import load_const
    U = C.alloc(n)
    load_const(C, U, mod.q)
    pool = list(X) + C.alloc(M.tape_external)
    M._euclid_pass(U, pool, C.alloc(M.held_count))


def bezout_only(C, Y, S, G):
    # The pass reads the coded tape; the B1 sweep now runs in Euclid's tail.
    M = DialogMultiplier(C, mod)
    M._bezout_pass(M.with_free_slots(G), Y, S)


def mulsub(C, g, X, Y):
    gen.C, gen._modadd_blk = C, None
    gen._mulsub(g[0], Y, X)


out["euclid_pass"] = count(euclid_only, [n, n])
tape = DialogMultiplier(Circuit(), mod)
freed = sum(len(tape._block_free(j, tape.swept)) for j in range(tape.blocks))
out["bezout_pass"] = count(bezout_only, [n, n, 2 * tick_bound(mod.q) - freed])
out["dialog_multiplier_M"] = count(lambda C, X, Y: DialogMultiplier(C, mod).apply(X, Y), [n, n])
out["lambda_square_subtract"] = count(mulsub, [1, n, n])
full = gen.build()
k = full.ops[0::4]
out["controlled_point_add"] = {"ops": len(full), "ccx": k.count(2), "cx": k.count(1), "x": k.count(0), "wires": full.num_wires}
print(json.dumps(out, indent=1))
