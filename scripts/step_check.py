"""Each Euclid step against the same step at a reference commit, on every input state.

    python3.12 -I -B scripts/step_check.py [--ref 3cb492c] [--max-width 7]

The depth changes rewrite the step's gates but not its map. For every width up to --max-width
and every top_clear, narrow, last, keep, sw and copy-count variant the build allows, both steps
run on all 2^(2w+1) values of (U, V, E), with the cell, kept and copy wires zero as the pass
gives them. On every state of the domain _top_bit's mode names (all states when the step has no
top mode), every output wire must match the reference step's, every copy wire must be zero, and
so must every other wire. This covers states the pass never reaches.

Controls: a copy of the current step with one CCX dropped must fail, for each of a few CCX.
"""
import argparse
import importlib
import io
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec import sim
from qec.arith import Mod
from qec.builder import CCX, Circuit
from qec.point_add import DialogMultiplier


def load_ref(commit, tmp, repo=ROOT):
    data = subprocess.run(["git", "-C", str(repo), "archive", commit, "qec"], check=True,
                          capture_output=True).stdout
    tarfile.open(fileobj=io.BytesIO(data)).extractall(tmp, filter="data")
    (Path(tmp) / "qec").rename(Path(tmp) / "qec_ref")
    sys.path.insert(0, tmp)
    return importlib.import_module("qec_ref.point_add"), importlib.import_module("qec_ref.builder"), \
        importlib.import_module("qec_ref.arith")


def domain(w, top_clear, narrow):
    """Lane mask of the states a top mode's _top_bit is exact on: top_clear, V's top bit is 0 on
    entry; narrow, the step leaves U's top bit 0, U' = min(U, V ^ E) when V ^ E is odd. Lane
    bits are U, then V, then E, as outputs packs them."""
    m = 0
    for lane in range(1 << (2 * w + 1)):
        u, v, e = lane & ((1 << w) - 1), lane >> w & ((1 << w) - 1), lane >> 2 * w
        x = v ^ e
        u_out = min(u, x) if x & 1 else u
        if top_clear and v >> (w - 1):
            continue
        if narrow and u_out >> (w - 1):
            continue
        m |= 1 << lane
    return m


def outputs(block, sizes, ops=None):
    """Run block on every (U, V, E) with the rest zero; return the output wires' lane masks."""
    w = sizes[0]
    lanes = 1 << (2 * w + 1)
    wires = block.flat_in
    num = max(max(block.wires, default=0), max(wires), max(block.flat_out)) + 1
    init = {}
    for k in range(2 * w + 1):
        m = 0
        for lane in range(lanes):
            if lane >> k & 1:
                m |= 1 << lane
        init[wires[k]] = m
    state = sim.run(block.ops if ops is None else ops, num, init, lanes)
    out = [state[x] for x in block.flat_out]
    scratch = [x for x in range(num) if x not in set(block.flat_out)]
    dirty = 0
    for x in scratch:
        dirty |= state[x]
    return out, dirty


def variants(max_w):
    for w in range(2, max_w + 1):
        for top_clear in (False, True):
            if top_clear and w < 3:
                continue
            for narrow in (False, True):
                for last in (False, True):
                    if last and (top_clear or narrow):
                        continue
                    kc = DialogMultiplier._kept_count(w, top_clear, narrow, last)
                    chain = (w - 1 if w >= 2 and (top_clear or narrow) else w)
                    keeps = [0] + [k for k in range(1, kc + 1) if 1 <= k <= chain - 1]
                    for keep in keeps:
                        for sw in (False, True):
                            if last and (keep or sw):
                                continue
                            full = 0 if last else max(DialogMultiplier._strip_copies(w, top_clear),
                                                      2 * DialogMultiplier.SELECT_COPIES)
                            for cp in sorted({0, min(1, full), min(3, full), full}):
                                yield w, top_clear, narrow, last, keep, sw, cp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="3cb492c")
    ap.add_argument("--max-width", type=int, default=7)
    ap.add_argument("--repo", default=str(ROOT), help="git checkout that holds --ref")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        ref_pa, ref_builder, ref_arith = load_ref(a.ref, tmp, a.repo)
        q = 251
        cur = DialogMultiplier(Circuit(), Mod(q))
        ref = ref_pa.DialogMultiplier(ref_builder.Circuit(), ref_arith.Mod(q))
        bad = total = 0
        controls_caught = controls = 0
        for v in variants(a.max_width):
            w, top_clear, narrow, last, keep, sw, cp = v
            new = cur._euclid_step(w, top_clear, narrow, last, keep, sw, cp)
            old = ref._euclid_step(w, top_clear, narrow, last, keep, sw)
            sizes = [w, w, 1, 4 if sw else 3]
            dom = domain(w, top_clear, narrow) if w >= 2 and not last else (1 << (2 * w + 1)) - 1
            o_new, dirty = outputs(new, sizes)
            o_old, dirty_old = outputs(old, sizes)
            o_new = [x & dom for x in o_new]
            o_old = [x & dom for x in o_old]
            dirty, dirty_old = dirty & dom, dirty_old & dom
            assert not dirty_old
            n_old = len(o_old)
            copies_clean = not any(o_new[n_old:])
            ok = o_new[:n_old] == o_old and copies_clean and not dirty
            total += 1
            if not ok:
                bad += 1
                print("MISMATCH", v, flush=True)
            # Controls: drop one gate at a few positions of the new step; each must differ.
            if cp and keep == 0 and not sw and w >= 4:
                ccx = [k for k in range(len(new.ops) // 4) if new.ops[4 * k] == CCX]
                for frac in (0.3, 0.6, 0.9):
                    k = ccx[int(len(ccx) * frac)]
                    ops = new.ops[:4 * k] + new.ops[4 * k + 4:]
                    o_mut, dirty_mut = outputs(new, sizes, ops)
                    o_mut, dirty_mut = [x & dom for x in o_mut], dirty_mut & dom
                    controls += 1
                    if o_mut[:n_old] != o_old or any(o_mut[n_old:]) or dirty_mut:
                        controls_caught += 1
        print(f"{total} variants, {bad} mismatched; controls caught {controls_caught} of {controls}")
        sys.exit(1 if bad or controls_caught != controls else 0)


if __name__ == "__main__":
    main()
