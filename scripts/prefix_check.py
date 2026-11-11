"""Prefix(bits) then write(t), against the same at a reference commit, on every input.

    python3.12 -I -B scripts/prefix_check.py [--ref 3cb492c] [--max-bits 14]

For every length h up to --max-bits, with and without release()/restore() between the two, both
versions run on all 2^(h+1) values of (bits, t) with clean scratch. The outputs must match and
every other wire must be zero; the result must also be bits + t mod 2^h. Control: the current
version with one CCX of write dropped must fail.
"""
import argparse
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from qec import sim
from qec.builder import CCX, Circuit
from qec.two_bit import Prefix
from step_check import load_ref


def build(Circ, Pre, h, release, drop=None):
    C = Circ()
    bits, (t,) = C.alloc(h), C.alloc(1)
    P = Pre(C, bits)
    if release:
        P.release()
        P.restore()
    start = len(C.prims) // 4
    P.write(t)
    ops = C.prims
    if drop is not None:
        ccx = [k for k in range(start, len(ops) // 4) if ops[4 * k] == CCX]
        if not ccx:
            return None
        k = ccx[int(len(ccx) * drop)]
        ops = ops[:4 * k] + ops[4 * k + 4:]
    return ops, bits, t, C.num_wires


def run(built, h):
    ops, bits, t, num = built
    lanes = 1 << (h + 1)
    init = sim.pack(range(lanes), bits + [t])
    state = sim.run(ops, num, init, lanes)
    out = sim.unpack(state, bits, lanes)
    dirty = 0
    for w in range(num):
        if w not in bits and w != t:
            dirty |= state[w]
    ok_t = state[t] == init[t]
    return out, dirty, ok_t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="3cb492c")
    ap.add_argument("--max-bits", type=int, default=14)
    ap.add_argument("--repo", default=str(ROOT), help="git checkout that holds --ref")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        load_ref(a.ref, tmp, a.repo)
        import qec_ref.builder as rb
        import qec_ref.two_bit as rt
        bad = cases = caught = controls = 0
        for h in range(1, a.max_bits + 1):
            for release in (False, True):
                new = run(build(Circuit, Prefix, h, release), h)
                old = run(build(rb.Circuit, rt.Prefix, h, release), h)
                want = [(v & ((1 << h) - 1)) + (v >> h) & ((1 << h) - 1) for v in range(1 << (h + 1))]
                cases += 1
                if new != old or new[0] != want or new[1] or not new[2]:
                    bad += 1
                    print("MISMATCH", h, release, flush=True)
                if h >= 3:
                    for frac in (0.0, 0.5, 0.99):
                        built = build(Circuit, Prefix, h, release, frac)
                        if built is None:
                            continue
                        controls += 1
                        mut = run(built, h)
                        caught += mut[0] != want or bool(mut[1]) or not mut[2]
        print(f"{cases} cases, {bad} mismatched; controls caught {caught} of {controls}")
        sys.exit(1 if bad or caught != controls else 0)


if __name__ == "__main__":
    main()
