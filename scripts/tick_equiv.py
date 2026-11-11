"""The full-width Bezout tick under the current two_bit flags against the same tick under REF's flags,
on random basis states, for every tick variant the pass records.

    python3.12 -I -B scripts/tick_equiv.py [lanes] [seed]

Every flag in REF only chooses among exact constructions of the same map (block sizes, tree
shapes, copies, gate order, wire choice), so the two must agree on every input, not just on the
promise domain: each variant's output registers must match lane for lane, and its scratch must
be clean wherever the reference's is. Exits 1 on any difference.
"""
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from qec import sim, two_bit, wide
from qec.builder import Circuit
from qec.two_bit import TwoBitMod

Q = 2**256 - 2**32 - 977

# The flags of b9b464a (PR #18), and the off value of each flag added since.
REF = dict(WIDE_ADD=7, WIDE_ADD_RAMP=2, WIDE_ADD_COPIES=16, WIDE_ADD_BACK=(10, 1), WIDE_WAVE=False,
           P_FAN=8, CMP_SPLIT=0, FINAL_EARLY=0, WAVE_SEG=0, ONE_BEFORE_DROP=False,
           IDLE_FROM_START=False, FINAL_INTERLEAVE=False, WIDE_ASIDE=False, ASIDE_UNDO_FIRST=False, SECOND_FRESH=False, IDLE_SLACK=0, KUPD_FAST=False, WIDE_M=3, WIDE_ADD_FSIZES=None,
           WIDE_ADD_BSIZES=None, LA_TREE_B=0, WIDE_WAVE_BACK=False, ROT_EXTRA=0, TOUCH_CAP=0,
           **{"wide.MID_AFTER_CARRIES": False, "wide.PEBBLE_MID": 0})


def _owner(name):
    return (wide, name[5:]) if name.startswith("wide.") else (two_bit, name)


def record(key, flags):
    saved = {k: getattr(*_owner(k)) for k in flags}
    for k, v in flags.items():
        setattr(*_owner(k), v)
    try:
        first, last, pebble, consume, produce = key
        mod = TwoBitMod(Q)
        n = mod.n
        C = Circuit()

        def tick(prev, g, nxt, R, S, pin, pout):
            mod.tick(C, prev or None, g, nxt[0] if nxt else None, R, S, pebble=pebble,
                     pin=pin[0] if pin else None, pout=pout[0] if pout else None)

        sizes = [0 if first else 2, 2, 0 if last else 1, n, n, int(consume), int(produce)]
        return C.record_shape(tick, sizes)
    finally:
        for k, v in saved.items():
            setattr(*_owner(k), v)


def run(block, init_vals, lanes):
    nw = max(max(block.wires), max(block.flat_in)) + 1
    init = {w: v for w, v in zip(block.flat_in, init_vals)}
    s = sim.run(block.ops, nw, init, lanes)
    out = [s[w] for w in block.flat_out]
    live = set(block.flat_out)
    dirty = 0
    for w in block.wires:
        if w not in live:
            dirty |= s[w]
    return out, dirty


def main():
    lanes = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    cur = {k: getattr(*_owner(k)) for k in REF}
    print("current flags:", {k: v for k, v in cur.items() if v != REF[k]})
    bad = 0
    for key in itertools.product([False, True], repeat=5):
        first, last, pebble, consume, produce = key
        if consume and first:
            continue
        ref, new = record(key, REF), record(key, cur)
        assert len(ref.flat_in) == len(new.flat_in)
        init = [rng.getrandbits(lanes) for _ in ref.flat_in]
        # pin holds the previous cell's B1 mark, NOR of prev: random marks break both alike, so
        # give it its promised value.
        if consume:
            pa, pb = init[0], init[1]
            init[len(init) - 1 - int(produce)] = ~(pa | pb) & ((1 << lanes) - 1)
        if produce:
            init[-1] = 0
        out_r, dirty_r = run(ref, init, lanes)
        out_n, dirty_n = run(new, init, lanes)
        diff = 0
        for a, b in zip(out_r, out_n):
            diff |= a ^ b
        extra = dirty_n & ~dirty_r
        status = "ok" if not diff and not extra else "DIFF"
        bad += status != "ok"
        print(f"{key} ops {len(ref.ops) // 4}/{len(new.ops) // 4} lanes differing "
              f"{bin(diff).count('1')} dirty ref {bin(dirty_r).count('1')} new {bin(dirty_n).count('1')} {status}")
    print("FAIL" if bad else "PASS")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
