"""Euclid tick 0 against tick 0 at a reference commit, on every V at toy widths.

    python3.12 -I -B scripts/tick0_check.py [--ref 98a5409]

Tick 0 runs on U = q, E = 0 and a zero cell. For every prime given whose first width is its bit
length, both blocks run on every V in 0 .. 2^n - 1, and at secp256k1 on 4,096 random V and the
edges; U, V, E and the cell must match and every other wire must end zero. Controls: the current
block with one CCX dropped, at a few places, must differ somewhere.
"""
import argparse
import random
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from step_check import load_ref
from qec import sim
from qec.arith import Mod
from qec.builder import CCX, Circuit
from qec.model import SECP256K1_Q
from qec.point_add import DialogMultiplier


def run(block, q, n, vs, ops=None):
    """Outputs (U, V, E, cell) per lane and the OR of every other wire, on U = q, E = cell = 0."""
    wires = block.flat_in
    num = max(max(block.wires), max(wires), max(block.flat_out)) + 1
    U, V = wires[:n], wires[n:2 * n]
    init = sim.pack([q] * len(vs), U) | sim.pack(vs, V)
    state = sim.run(block.ops if ops is None else ops, num, init, len(vs))
    out = block.flat_out[:2 * n + 3]
    dirty = 0
    for w in range(num):
        if w not in set(out):
            dirty |= state[w]
    return [state[w] for w in out], dirty


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="98a5409")
    ap.add_argument("--repo", default=str(ROOT))
    ap.add_argument("--primes", type=int, nargs="*", default=[13, 31, 61, 127, 251, 509, 1021, 4093, 65521])
    a = ap.parse_args()
    ok = True
    caught = controls = 0
    with tempfile.TemporaryDirectory() as tmp:
        ref_pa, ref_builder, ref_arith = load_ref(a.ref, tmp, a.repo)
        rng = random.Random(4748)
        for q in a.primes + [SECP256K1_Q]:
            n = q.bit_length()
            cur = DialogMultiplier(Circuit(), Mod(q))
            if cur.widths[0] != n:
                continue
            new = cur._euclid_tick0()
            old = ref_pa.DialogMultiplier(ref_builder.Circuit(), ref_arith.Mod(q))._euclid_tick0()
            assert new.wires <= set(new.flat_in) and set(new.flat_out) == set(new.flat_in)
            if n <= 16:
                vs = list(range(1 << n))
            else:
                vs = sorted({0, 1, 2, q - 1, q - 2, (1 << n) - 1} | {rng.randrange(1 << n) for _ in range(4096)})
            o_new, d_new = run(new, q, n, vs)
            o_old, d_old = run(old, q, n, vs)
            assert not d_old
            good = o_new == o_old and not d_new
            ok &= good
            ccx = [k for k in range(len(new.ops) // 4) if new.ops[4 * k] == CCX]
            hit = 0
            for frac in (0.1, 0.4, 0.7, 0.95):
                k = ccx[int(len(ccx) * frac)]
                ops = new.ops[:4 * k] + new.ops[4 * k + 4:]
                o_mut, d_mut = run(new, q, n, vs, ops)
                controls += 1
                hit += o_mut != o_old or bool(d_mut)
            caught += hit
            print(f"q={'secp256k1' if q == SECP256K1_Q else q}: {len(vs)} V, "
                  f"CCX {old.ops[0::4].count(CCX)} -> {new.ops[0::4].count(CCX)}, "
                  f"{'match' if good else 'MISMATCH'}, controls {hit} of 4", flush=True)
    print(f"tick0_check: {'PASS' if ok and caught == controls else 'FAIL'}; controls caught {caught} of {controls}")
    sys.exit(0 if ok and caught == controls else 1)


if __name__ == "__main__":
    main()
