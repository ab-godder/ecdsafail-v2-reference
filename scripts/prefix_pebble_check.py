"""wide.PrefixTree.write in every write mode, on every input up to 13 bits and on drawn inputs at the
tick's 221 bits: bits end as bits + t mod 2^h, t and the mid probe are as stated, and every other
wire ends clean. Runs with and without the tick's idle allocator, which the pebbled write's wire
choice reads. Exits 1 on any failure, and on a mutant the check does not see.

    python3.12 -I -B scripts/prefix_pebble_check.py
"""
import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim, wide
from qec.builder import Circuit
from qec.two_bit import _Idle

MODES = [dict(pebble=None), dict(pebble=-1), dict(pebble=-1, fresh=99), dict(pebble=1, fresh=3),
         dict(pebble=99), dict(pebble=99, fresh=5), dict(pebble=99, fresh=99)]


def build(h, block, mode, idle, mutate=None):
    C = Circuit()
    bits, (t, probe) = C.alloc(h), C.alloc(2)
    spare = C.alloc(6)            # clean wires on the stack before the tree, as in the tick
    C.free(spare)
    if idle:
        I = _Idle(C, 0)
        C.alloc, C.alloc_after = I.alloc, I.alloc_after
    P = wide.PrefixTree(C, bits, block, **mode)
    if mutate is not None:
        mutate(C)
    seen = []

    def mid(z):
        seen.append(z)
        C.cx(z, probe)

    P.write(t, mid)
    if idle:
        del C.alloc, C.alloc_after
    assert len(seen) == 1, "mid must run exactly once"
    return C, bits, t, probe


def check(h, block, mode, idle, values, mutate=None):
    """values: list of (bits, t). Returns the number of failing lanes."""
    C, bits, t, probe = build(h, block, mode, idle, mutate)
    lanes = len(values)
    init = sim.pack([v for v, _ in values], bits)
    init.update(sim.pack([c for _, c in values], [t]))
    state = sim.run(C.ops, C.num_wires, init, lanes)
    got = sim.unpack(state, bits, lanes)
    got_t = sim.unpack(state, [t], lanes)
    got_p = sim.unpack(state, [probe], lanes)
    mask = (1 << h) - 1
    fails = 0
    for lane, (v, c) in enumerate(values):
        new = (v + c) & mask
        if got[lane] != new or got_t[lane] != c or got_p[lane] != int((~new & mask) == mask):
            fails += 1
    dirty = [w for w in range(C.num_wires) if w not in set(bits) | {t, probe} and state[w]]
    return fails + (lanes if dirty else 0)


def layouts(h):
    out = [1, 2, 3, 4, 5]
    if h >= 6:
        out += [(2, 3), (4, 3, 2), (5, 4)]
    return out


def main():
    rng = random.Random(1)
    total = cases = 0
    for h in range(1, 14):
        values = [(v, c) for v in range(1 << h) for c in (0, 1)]
        for block, mode, idle in itertools.product(layouts(h), MODES, (False, True)):
            total += check(h, block, mode, idle, values)
            cases += len(values)
    # The tick's own width, with the lowest zero at every position and drawn bits above it.
    h = 221
    values = []
    for z in range(h + 1):
        for _ in range(4):
            hi = rng.getrandbits(h) << (z + 1)
            values.append((((1 << z) - 1 | hi) & ((1 << h) - 1), rng.getrandbits(1)))
    values += [((1 << h) - 1, 1), ((1 << h) - 1, 0), (0, 1), (0, 0)]
    for block, mode, idle in itertools.product((4, 3, (4,) * 54 + (5,)), MODES, (False, True)):
        total += check(h, block, mode, idle, values)
        cases += len(values)
    print(f"prefix write: {cases} lanes, {total} failures")
    # Control: dropping any one CCX of a pebbled write must be seen.
    caught = tried = 0
    values = [(v, c) for v in range(1 << 9) for c in (0, 1)]
    for mode in (dict(pebble=99, fresh=99), dict(pebble=-1)):
        writes = []

        def count(C):
            orig = C.ccx

            def ccx(a, b, c):
                writes.append(1)
                orig(a, b, c)
            C.ccx = ccx
        build(9, 2, mode, True, count)
        for drop_at in range(len(writes)):
            def mutate(C, drop_at=drop_at):
                orig = C.ccx
                seen = [0]

                def ccx(a, b, c):
                    seen[0] += 1
                    if seen[0] != drop_at + 1:
                        orig(a, b, c)
                C.ccx = ccx
            try:
                bad = check(9, 2, mode, True, values, mutate)
            except AssertionError:
                bad = 1
            tried += 1
            caught += bad > 0
    print(f"mutants: {caught} of {tried} dropped CCX seen")
    sys.exit(0 if total == 0 and caught == tried else 1)


if __name__ == "__main__":
    main()
