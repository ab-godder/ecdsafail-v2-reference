"""The two-wire-cell Euclid step against the same step at a reference commit, on every input state,
to widths where the strip's copy fan has more than two copies.

    python3.12 -I -B scripts/strip_check.py [--ref 74f114a] [--min-width 2] [--max-width 9]

The step's tail was reordered for depth: the cell's swap before the strip, the stripped bit cleared
at the top after the first swap layer, the strip's copies taken from the far end of copies and fanned
from a second copy, and the fans after the unwind run on g1's and g0's copies alone. None of it may
change the map. For every width up to --max-width, top_clear, narrow and last, and copy counts
0, 1, 3 and the pass's full count, the step runs on all 2^(2w+1) values of (U, V, E) with the cell
and copies zero. On every state of the domain _top_bit's mode names, every output wire must match
the reference step's, the copies must be zero, and so must every other wire.

Controls, each of which must fail: copies of the current step with one CCX dropped, at the cell's
swap and at the strip's clear and at a few other places; and a copy with the strip's clear moved
to another wire of V."""
import argparse
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import functools

from step_check import load_ref
from qec import sim
from qec.arith import Mod
from qec.builder import CCX, Circuit
from qec.point_add import DialogMultiplier


@functools.lru_cache(maxsize=None)
def bit_masks(nbits):
    """Per input bit k, the lanes whose index has bit k set, built in linear time (step_check's
    quadratic loop is too slow past width 8)."""
    lanes = 1 << nbits
    out = []
    for k in range(nbits):
        period = 1 << (k + 1)
        block = ((1 << (1 << k)) - 1) << (1 << k)
        out.append(block * (((1 << lanes) - 1) // ((1 << period) - 1)))
    return out


def lanes_where(nbits, pred):
    """Lane mask of the lanes whose index satisfies pred."""
    buf = bytearray((1 << nbits) // 8 or 1)
    for lane in range(1 << nbits):
        if pred(lane):
            buf[lane >> 3] |= 1 << (lane & 7)
    return int.from_bytes(buf, "little")


@functools.lru_cache(maxsize=None)
def domain(w, top_clear, narrow):
    """step_check.domain, the same states, built in linear time."""
    def ok(lane):
        u, v, e = lane & ((1 << w) - 1), lane >> w & ((1 << w) - 1), lane >> 2 * w
        x = v ^ e
        u_out = min(u, x) if x & 1 else u
        return not (top_clear and v >> (w - 1)) and not (narrow and u_out >> (w - 1))
    return lanes_where(2 * w + 1, ok)


def outputs(block, sizes, ops=None):
    """step_check.outputs with linear-time input masks."""
    w = sizes[0]
    wires = block.flat_in
    num = max(max(block.wires, default=0), max(wires), max(block.flat_out)) + 1
    masks = bit_masks(2 * w + 1)
    init = {wires[k]: masks[k] for k in range(2 * w + 1)}
    state = sim.run(block.ops if ops is None else ops, num, init, 1 << (2 * w + 1))
    out = [state[x] for x in block.flat_out]
    live = set(block.flat_out)
    dirty = 0
    for x in range(num):
        if x not in live:
            dirty |= state[x]
    return out, dirty


def variants(min_w, max_w):
    for w in range(min_w, max_w + 1):
        for top_clear in (False, True):
            if top_clear and w < 3:
                continue
            for narrow in (False, True):
                for last in (False, True):
                    if last and (top_clear or narrow):
                        continue
                    full = 0 if last else max(DialogMultiplier._strip_copies(w, top_clear),
                                              2 * DialogMultiplier.SELECT_COPIES)
                    for cp in sorted({0, min(1, full), min(3, full), full}):
                        yield w, top_clear, narrow, last, cp


def tail_ccx(block, wires):
    """Positions of the CCX whose target is one of wires, in op order."""
    ops = block.ops
    return [k for k in range(len(ops) // 4) if ops[4 * k] == CCX and ops[4 * k + 3] in wires]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="74f114a")
    ap.add_argument("--min-width", type=int, default=2)
    ap.add_argument("--max-width", type=int, default=9)
    a = ap.parse_args()
    q = 251
    with tempfile.TemporaryDirectory() as tmp:
        ref_pa, ref_builder, ref_arith = load_ref(a.ref, tmp)
        cur = DialogMultiplier(Circuit(), Mod(q))
        ref = ref_pa.DialogMultiplier(ref_builder.Circuit(), ref_arith.Mod(q))
        cur.HELD = ref.HELD = False
        total = bad = controls = caught = 0
        for v in variants(a.min_width, a.max_width):
            w, top_clear, narrow, last, cp = v
            new = cur._euclid_step(w, top_clear, narrow, last, 0, False, cp)
            old = ref._euclid_step(w, top_clear, narrow, last, 0, False, cp)
            dom = domain(w, top_clear, narrow) if w >= 2 and not last else (1 << (2 * w + 1)) - 1
            o_new, dirty = outputs(new, [w])
            o_old, dirty_old = outputs(old, [w])
            o_new, o_old = [x & dom for x in o_new], [x & dom for x in o_old]
            assert not dirty_old & dom
            n = 2 * w + 3
            ok = o_new[:n] == o_old[:n] and not any(o_new[n:]) and not dirty & dom
            total += 1
            if not ok:
                bad += 1
                print("MISMATCH", v, flush=True)
            if last or w < 4:
                continue
            flat = new.flat_in
            V, E, cell = flat[w:2 * w], flat[2 * w], flat[2 * w + 1:2 * w + 3]
            mutants = []
            # The cell's swap's CCX targets b; the strip's clear targets V's wire that ends at the top.
            for k in tail_ccx(new, {cell[1]})[-1:]:
                mutants.append(("swap", new.ops[:4 * k] + new.ops[4 * k + 4:]))
            clear = [k for k in tail_ccx(new, set(V)) if new.ops[4 * k + 2] == E]
            for k in clear[-1:]:
                mutants.append(("clear", new.ops[:4 * k] + new.ops[4 * k + 4:]))
                moved = list(new.ops)
                moved[4 * k + 3] = next(x for x in V if x != new.ops[4 * k + 3])
                mutants.append(("moved clear", type(new.ops)("i", moved)))
            ccx = [k for k in range(len(new.ops) // 4) if new.ops[4 * k] == CCX]
            for frac in (0.3, 0.9):
                k = ccx[int(len(ccx) * frac)]
                mutants.append((f"drop {frac}", new.ops[:4 * k] + new.ops[4 * k + 4:]))
            for name, ops in mutants:
                o_mut, dirty_mut = outputs(new, [w], ops)
                o_mut = [x & dom for x in o_mut]
                controls += 1
                if o_mut[:n] != o_old[:n] or any(o_mut[n:]) or dirty_mut & dom:
                    caught += 1
                else:
                    print("control not caught", name, v, flush=True)
        print(f"{total} variants, {bad} mismatched; controls caught {caught} of {controls}")
        sys.exit(1 if bad or caught != controls else 0)


if __name__ == "__main__":
    main()
