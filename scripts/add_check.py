"""The tick's add-or-shift, every variant, on every input of its domain.

    python3.12 -I -B scripts/add_check.py [--max-bits 7]

For n from 3 to --max-bits, held_lane on and off, and each SPLIT_ADD / SPLIT_COPIES setting, the
add runs on all (J, a, word) with a < 2^n and word < 2^(n + 2), clean scratch. On the domain (word's
top two bits zero on J, and everywhere without held_lane) it must give word = 2*word on J and
word = (word + (a >> 1)) mod 2^(n + 2) off J, with J and a unchanged and every other wire zero.
Controls, which must fail: the split add with one CCX dropped, for a few CCX, and with its up pass's
carry-line XORs emitted top down, which reads a line after it is written. The one-pass add with one
CCX dropped, with its parity D's last read dropped, and with each line restored only after the sum
reads it (so the sum takes l_i = c_i ^ o_i where it needs c_i).

The hand-over variant (g0 given) runs every add on all (g0, J, a, word) with g0 = 1 only off J. On
a B1 lane (g0 = J = 0, word's top two bits zero) the sum's bits n - 1 and n must leave word for the
returned pair (hn, hb), which is 0 on every other lane; g0 must come back unchanged. Its controls
drop one CCX of the add top's hand-over, which must fail.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.two_bit import TwoBitMod

PRIMES = {3: 5, 4: 13, 5: 29, 6: 61, 7: 113, 8: 251, 9: 509}


def build(mod, held, mutate=None):
    n = mod.n
    C = Circuit()
    (J,), a, word = C.alloc(1), C.alloc(n), C.alloc(n + 2)
    mod._add_or_shift(C, J, a, word, held_lane=held)
    ops = C.prims
    if mutate is not None:
        ops = mutate(ops)
    return ops, J, a, word, C.num_wires


def run(built, n, held):
    """Lanes of the domain where some wire is wrong."""
    ops, J, a, word, num = built
    wires = [J] + a + word
    lanes = 1 << (2 * n + 3)
    state = sim.run(ops, num, sim.pack(range(lanes), wires), lanes)
    want, domain = [], 0
    for v in range(lanes):
        j, av, wv = v & 1, v >> 1 & ((1 << n) - 1), v >> (n + 1)
        if not ((j or not held) and wv >> n):
            domain |= 1 << v
        want_w = 2 * wv if j else (wv + (av >> 1)) % (1 << (n + 2))
        want.append(j | av << 1 | (want_w % (1 << (n + 2))) << (n + 1))
    exp = sim.pack(want, wires)
    wrong = 0
    for w in range(num):
        wrong |= state[w] ^ exp.get(w, 0)
    return bin(wrong & domain).count("1")


def build_hand(mod, held, mutate=None):
    n = mod.n
    C = Circuit()
    (g0, J), a, word = C.alloc(2), C.alloc(n), C.alloc(n + 2)
    hn, hb = mod._add_or_shift(C, J, a, word, held_lane=held, g0=g0)
    ops = C.prims
    if mutate is not None:
        ops = mutate(ops)
    return ops, g0, J, a, word, hn, hb, C.num_wires


def run_hand(built, n, held):
    """Lanes of the hand-over variant's domain where some wire is wrong."""
    ops, g0, J, a, word, hn, hb, num = built
    wires = [g0, J] + a + word
    lanes = 1 << (2 * n + 4)
    state = sim.run(ops, num, sim.pack(range(lanes), wires), lanes)
    want, domain = {}, 0
    for v in range(lanes):
        g, j = v & 1, v >> 1 & 1
        av, wv = v >> 2 & ((1 << n) - 1), v >> (n + 2)
        b1 = not g and not j
        if not (g and j) and not ((j or b1 or not held) and wv >> n):
            domain |= 1 << v
        res = 2 * wv if j else (wv + (av >> 1)) % (1 << (n + 2))
        top = res >> (n - 1) & 3 if b1 else 0
        res ^= top << (n - 1)
        for w, bit in zip([g0, J] + a + word + [hn, hb],
                          [g, j] + [av >> i & 1 for i in range(n)]
                          + [res >> i & 1 for i in range(n + 2)] + [top & 1, top >> 1]):
            if bit:
                want[w] = want.get(w, 0) | 1 << v
    wrong = 0
    for w in range(num):
        wrong |= state[w] ^ want.get(w, 0)
    return bin(wrong & domain).count("1")


def drop_ccx_onto(t, last=False):
    """Drop the first CCX that targets wire t, or the last with last."""
    def mutate(ops):
        hits = [k for k in range(len(ops) // 4) if ops[4 * k] == CCX and ops[4 * k + 3] == t]
        k = hits[-1 if last else 0]
        return ops[:4 * k] + ops[4 * k + 4:]
    return mutate


def drop_cx_onto(t):
    """Drop the last CX that targets wire t."""
    def mutate(ops):
        k = [k for k in range(len(ops) // 4) if ops[4 * k] == CX and ops[4 * k + 2] == t][-1]
        return ops[:4 * k] + ops[4 * k + 4:]
    return mutate


def drop_ccx(frac):
    def mutate(ops):
        ccx = [k for k in range(len(ops) // 4) if ops[4 * k] == CCX]
        k = ccx[int(len(ccx) * frac)]
        return ops[:4 * k] + ops[4 * k + 4:]
    return mutate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-bits", type=int, default=7)
    args = ap.parse_args()
    cases = failed = controls = caught = 0
    for n in range(3, args.max_bits + 1):
        mod = TwoBitMod(PRIMES[n])
        assert mod.n == n
        for held in (False, True):
            for split, copies in ((False, 0), (True, 1), (True, 2), (True, 4), (True, 7)):
                mod.SPLIT_ADD, mod.SPLIT_COPIES = split, copies
                bad = run(build(mod, held), n, held)
                cases += 1
                failed += bad > 0
                if bad:
                    print(f"FAIL n={n} held={held} split={split} copies={copies}: {bad} lanes")
            mod.SPLIT_ADD, mod.SPLIT_COPIES = True, 4
            for frac in (0.0, 0.3, 0.6, 0.99):
                controls += 1
                hit = run(build(mod, held, drop_ccx(frac)), n, held) > 0
                caught += hit
                if not hit:
                    print(f"control not caught: drop {frac}, n={n} held={held}")
            # Top-down carry-line XORs: each c ^= a then reads a line the one above already wrote.
            orig = mod._add_then_rotate

            def top_down(C, J, a, word, held_lane=True, g0=None):
                real_cx = C.cx
                ops_ = a[1:]
                lines = set(zip(ops_, [None] + ops_[:-1]))
                held_back = []

                def cx(c, t):
                    if (c, t) in lines and len(held_back) < len(ops_) - 1:
                        held_back.append((c, t))
                        if len(held_back) == len(ops_) - 1:
                            for c2, t2 in reversed(held_back):
                                real_cx(c2, t2)
                        return
                    real_cx(c, t)
                C.cx = cx
                out = orig(C, J, a, word, held_lane, g0)
                del C.cx
                return out
            if n >= 4:
                mod._add_then_rotate = top_down
                controls += 1
                hit = run(build(mod, held), n, held) > 0
                caught += hit
                if not hit:
                    print(f"control not caught: top-down lines, n={n} held={held}")
                del mod._add_then_rotate
            # The one-pass add: one CCX dropped; the parity D left with one line in it; and the
            # sum reading line i before its restore (l_i, not c_i), the order the restore-first
            # down pass exists to avoid.
            mod.SPLIT_ADD = False
            for frac in (0.1, 0.5, 0.9):
                controls += 1
                hit = run(build(mod, held, drop_ccx(frac)), n, held) > 0
                caught += hit
                if not hit:
                    print(f"control not caught: one-pass drop {frac}, n={n} held={held}")
            D = 1 + n + (n + 2) + 1
            controls += 1
            hit = run(build(mod, held, drop_cx_onto(D)), n, held) > 0
            caught += hit
            if not hit:
                print(f"control not caught: one-pass D read dropped, n={n} held={held}")
            orig_d = mod._add_or_shift_d

            def late_restore(C, J, a, word, held_lane=True, g0=None):
                real_cx, real_ccx = C.cx, C.ccx
                ops_ = a[1:]
                restore = {(ops_[i], ops_[i - 1]): (ops_[i - 1], word[i]) for i in range(1, n - 1)}
                state = {"down": False, "held": {}}

                def ccx(x, y, z):
                    state["down"] |= x == J
                    real_ccx(x, y, z)

                def cx(x, y):
                    if state["down"] and (x, y) in restore:
                        state["held"][restore[(x, y)]] = (x, y)
                        return
                    real_cx(x, y)
                    if (x, y) in state["held"]:
                        real_cx(*state["held"].pop((x, y)))
                C.cx, C.ccx = cx, ccx
                out = orig_d(C, J, a, word, held_lane, g0)
                del C.cx, C.ccx
                return out
            mod._add_or_shift_d = late_restore
            controls += 1
            hit = run(build(mod, held), n, held) > 0
            caught += hit
            if not hit:
                print(f"control not caught: one-pass late restore, n={n} held={held}")
            del mod._add_or_shift_d
            mod.SPLIT_ADD = True
            for split, copies in ((False, 0), (True, 1), (True, 4)):
                mod.SPLIT_ADD, mod.SPLIT_COPIES = split, copies
                bad = run_hand(build_hand(mod, held), n, held)
                cases += 1
                failed += bad > 0
                if bad:
                    print(f"FAIL hand n={n} held={held} split={split} copies={copies}: {bad} lanes")
                # The hand-over's two CCX on its wires: the one that writes hn, and t's clear.
                plain = sum(op == CCX for op in build(mod, held)[0][::4])
                ops, *_, hn, hb, _ = build_hand(mod, held)
                assert sum(op == CCX for op in ops[::4]) == plain + 1
                for mutate in (drop_ccx_onto(hn), drop_ccx_onto(hb, last=True)):
                    controls += 1
                    hit = run_hand(build_hand(mod, held, mutate), n, held) > 0
                    caught += hit
                    if not hit:
                        print(f"control not caught: hand split={split}, n={n} held={held}")
            mod.SPLIT_ADD, mod.SPLIT_COPIES = True, 4
    print(f"{cases} cases, {failed} failed; controls caught {caught} of {controls}")
    sys.exit(1 if failed or caught != controls else 0)


if __name__ == "__main__":
    main()
