"""Prefix(bits) then write(t), against the same at a reference commit, on every input.

    python3.12 -I -B scripts/prefix_check.py [--ref 3cb492c] [--max-bits 14]

For every length h up to --max-bits, with and without release()/restore() between the two, both
versions run on all 2^(h+1) values of (bits, t) with clean scratch. The outputs must match and
every other wire must be zero; the result must also be bits + t mod 2^h. With each block cap in
CAPS that splits a block, the current version alone must give bits + t with clean scratch. Control:
the current version with one CCX of write dropped must fail, capped or not. The tick's lengths
(WIDE) run at their caps on sampled lanes: a run of ones of every length from bit 0, which
reaches every block's top, under random high bits, with t = 0 and 1, and 512 random lanes.
"""
import argparse
import random
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


def build(Circ, Pre, h, release, drop=None, cap=None):
    C = Circ()
    bits, (t,) = C.alloc(h), C.alloc(1)
    P = Pre(C, bits) if cap is None else Pre(C, bits, cap=cap)
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


def run(built, h, values=None):
    ops, bits, t, num = built
    values = range(1 << (h + 1)) if values is None else values
    lanes = len(values)
    init = sim.pack(values, bits + [t])
    state = sim.run(ops, num, init, lanes)
    out = sim.unpack(state, bits, lanes)
    dirty = 0
    for w in range(num):
        if w not in bits and w != t:
            dirty |= state[w]
    ok_t = state[t] == init[t]
    return out, dirty, ok_t


CAPS = (1, 2, 3, 5)
WIDE = ((221, None), (221, 37), (221, 23), (222, 37), (190, 37), (64, 13))


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
                # Capped blocks (P_CAP), against bits + t alone: the reference has no cap.
                for cap in CAPS:
                    if cap >= h - 1:
                        continue
                    got = run(build(Circuit, Prefix, h, release, cap=cap), h)
                    cases += 1
                    if got[0] != want or got[1] or not got[2]:
                        bad += 1
                        print("MISMATCH cap", h, release, cap, flush=True)
                    for frac in (0.0, 0.99):
                        built = build(Circuit, Prefix, h, release, frac, cap=cap)
                        if built is not None:
                            controls += 1
                            mut = run(built, h)
                            caught += mut[0] != want or bool(mut[1]) or not mut[2]
                if h >= 3:
                    for frac in (0.0, 0.5, 0.99):
                        built = build(Circuit, Prefix, h, release, frac)
                        if built is None:
                            continue
                        controls += 1
                        mut = run(built, h)
                        caught += mut[0] != want or bool(mut[1]) or not mut[2]
        # The tick's lengths at its caps, on lanes that reach every block's top: runs of ones from
        # bit 0 of every length, with t = 1 and 0, and random words.
        rng = random.Random(7)
        for h, cap in WIDE:
            vals = [((1 << k) - 1) | ((rng.getrandbits(h) << (k + 1)) % (1 << h)) | tb << h
                    for k in range(h + 1) for tb in (0, 1)]
            vals += [rng.getrandbits(h + 1) for _ in range(512)]
            want = [(v & ((1 << h) - 1)) + (v >> h) & ((1 << h) - 1) for v in vals]
            for release in (False, True):
                got = run(build(Circuit, Prefix, h, release, cap=cap), h, vals)
                cases += 1
                if got[0] != want or got[1] or not got[2]:
                    bad += 1
                    print("MISMATCH wide", h, cap, release, flush=True)
                built = build(Circuit, Prefix, h, release, 0.5, cap=cap)
                controls += 1
                mut = run(built, h, vals)
                caught += mut[0] != want or bool(mut[1]) or not mut[2]
        print(f"{cases} cases, {bad} mismatched; controls caught {caught} of {controls}")
        sys.exit(1 if bad or caught != controls else 0)


if __name__ == "__main__":
    main()
