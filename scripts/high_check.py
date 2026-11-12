"""Checks of the high-word Euclid step and its stretch plan (DialogMultiplier.high_plan).

    python3.12 -I -B scripts/high_check.py [--max-width 7] [--primes 251 509 ...]

1. Step. For every width w up to --max-width, every high length m (1 .. w - 1), with and without
   top_clear and narrow, init or not: the high step runs on all 2^(2w+2) values of (U, V, E, label).
   On every state where U and V are not both 2^K or more (K = w - m), the label names the one that
   is (either when neither is), and the plain step's top mode is exact, U, V, E and the cell must
   match the plain step, the label must end as label ^ g1 (g1 = g0 AND E after the step), with init
   [U >= 2^K] ^ g1, and every other wire must be zero.
   Controls, each of which must fail: a copy with one CCX dropped, at a few places; the right
   states with the label flipped where the word is nonzero; states where both words are nonzero.
2. Multiplier. The toy primes' own plans are empty, so the plans here are forced through
   HIGH_STRETCHES: every tick its own stretch, the longest stretches, and seeded random cuts.
   M runs forward and inverse on every x against six y; Y must be right and every other wire zero.
   Controls, each of which must fail somewhere: a stretch one bit under its joint threshold, and a
   step that skips the label's update.
3. Plan. At secp256k1 and the toy primes, the shipped plan's label is simulated on the physical
   states of the pass (qec.model.tick_pass) for sampled x: at each stretch tick, U and V must not
   both reach 2^K and the label must name the one that does.
"""
import argparse
import functools
import io
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod
from qec.builder import CCX, CX, Circuit
from qec.emit import write_circuit
from qec.model import SECP256K1_Q, tick_bound, tick_joint, tick_pass
from qec import point_add
from qec.point_add import DialogMultiplier, Stretch


@functools.lru_cache(maxsize=None)
def lanes_of(nbits):
    """Per input bit k, the lane mask of lanes whose index has bit k set."""
    lanes = 1 << nbits
    out = []
    for k in range(nbits):
        m = 0
        for lane in range(lanes):
            if lane >> k & 1:
                m |= 1 << lane
        out.append(m)
    return out


def run_step(block, w, label_in, ops=None):
    """Run a step on every (U, V, E, z); returns (state, lanes). Lane bits: U, V, E, then z."""
    nbits = 2 * w + 1 + (1 if label_in else 0)
    masks = lanes_of(nbits)
    wires = block.flat_in
    num = max(max(block.wires, default=0), max(wires), max(block.flat_out)) + 1
    init = {wires[k]: masks[k] for k in range(2 * w + 1)}
    if label_in:
        init[wires[-1]] = masks[2 * w + 1]
    return sim.run(block.ops if ops is None else ops, num, init, 1 << nbits), 1 << nbits


def domain_masks(w, m, top_clear, narrow, init):
    """(good, flipped, joint_broken) lane masks over (U, V, E, z)."""
    K = w - m
    good = flipped = broken = 0
    for lane in range(1 << (2 * w + 2)):
        u, v, e = lane & ((1 << w) - 1), lane >> w & ((1 << w) - 1), lane >> 2 * w & 1
        z = lane >> (2 * w + 1) & 1
        if init and z:
            continue
        x = v ^ e
        u_out = min(u, x) if x & 1 else u
        if top_clear and v >> (w - 1):
            continue
        if narrow and u_out >> (w - 1):
            continue
        uh, vh = u >> K, v >> K
        if uh and vh:
            broken |= 1 << lane
            continue
        if init:
            good |= 1 << lane
            continue
        right = (uh and z) or (vh and not z) or (not uh and not vh)
        if right:
            good |= 1 << lane
        else:
            flipped |= 1 << lane
    return good, flipped, broken


def step_part(max_w):
    q = 251
    cur = DialogMultiplier(Circuit(), Mod(q))
    cur.HELD = False
    total = bad = 0
    controls = caught = 0
    for w in range(2, max_w + 1):
        for m in range(1, w):
            K = w - m
            cp = max(cur._high_scratch(m), max(cur._strip_copies(w, False), 2 * cur.SELECT_COPIES))
            for top_clear in (False, True):
                if top_clear and w < 3:
                    continue
                for narrow in (False, True):
                    for init in (False, True):
                        new = cur._euclid_step(w, top_clear, narrow, False, 0, False, cp, (m, init))
                        ref = cur._euclid_step(w, top_clear, narrow, False, 0, False, cp)
                        good, flipped, broken = domain_masks(w, m, top_clear, narrow, init)
                        s_new, lanes = run_step(new, w, True)
                        s_ref, _ = run_step(ref, w, False)
                        # The reference ignores the label bit, so its lanes repeat in two halves.
                        half = 1 << (2 * w + 1)
                        ref_out = [s_ref[x] | s_ref[x] << half for x in ref.flat_out]
                        new_out = [s_new[x] for x in new.flat_out]
                        U, V, E, cell = (slice(0, w), slice(w, 2 * w), 2 * w, (2 * w + 1, 2 * w + 2))
                        n_ref = 2 * w + 3
                        g0, e_out = new_out[cell[0]], new_out[E]
                        g1 = g0 & e_out
                        uhigh = 0
                        for i in range(K, w):
                            uhigh |= lanes_of(2 * w + 2)[i]
                        z_in = 0 if init else lanes_of(2 * w + 2)[2 * w + 1]
                        want_z = (uhigh if init else z_in) ^ g1
                        z_out = new_out[-1]
                        live = set(new.flat_out[:n_ref]) | {new.flat_out[-1]}
                        num = len(s_new)
                        dirty = 0
                        for x in range(num):
                            if x not in live:
                                dirty |= s_new[x]

                        def wrong(mask, state_out=new_out, dirty=dirty, z_out=z_out):
                            d = dirty
                            for a, b in zip(state_out[:n_ref], ref_out):
                                d |= a ^ b
                            d |= z_out ^ want_z
                            return d & mask

                        total += 1
                        if wrong(good):
                            bad += 1
                            print("MISMATCH", w, m, top_clear, narrow, init, flush=True)
                        if not init and flipped:
                            controls += 1
                            caught += bool(wrong(flipped))
                            if not wrong(flipped):
                                print("uncaught flipped-label control", w, m, top_clear, narrow)
                        if broken:
                            controls += 1
                            caught += bool(wrong(broken))
                            if not wrong(broken):
                                print("uncaught joint control", w, m, top_clear, narrow, init)
                        # At K = 1 an init label is 0 whenever the word is, and then no borrow
                        # can leave the top: dropping that chain's last link is invisible there.
                        if w >= 4 and m >= 2 and not top_clear and not narrow and (K >= 2 or not init):
                            ccx = [k for k in range(len(new.ops) // 4) if new.ops[4 * k] == CCX]
                            for frac in (0.2, 0.5, 0.8):
                                k = ccx[int(len(ccx) * frac)]
                                ops = new.ops[:4 * k] + new.ops[4 * k + 4:]
                                s_mut, _ = run_step(new, w, True, ops)
                                mut_out = [s_mut[x] for x in new.flat_out]
                                d = 0
                                for x in range(len(s_mut)):
                                    if x not in live:
                                        d |= s_mut[x]
                                controls += 1
                                hit = bool(wrong(good, mut_out, d, mut_out[-1]))
                                caught += hit
                                if not hit:
                                    print("uncaught drop", w, m, init, k, new.ops[4 * k:4 * k + 4].tolist())
    print(f"step: {total} variants, {bad} mismatched; controls caught {caught} of {controls}",
          flush=True)
    return bad == 0 and caught == controls


def m_bad(q, ys, inverse, xs=None):
    pairs = [(x, y) for x in (xs or range(1, q)) for y in ys]
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    DialogMultiplier(C, mod).apply(X, Y, inverse=inverse)
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y}, depth_pass=False)
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    xs, values = zip(*pairs)
    init = sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1])
    state = sim.run(ops, C.num_wires, init, len(pairs))
    expected = [y * (pow(x, -1, q) if inverse else x) % q for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def eligible(m, K_of):
    """Ticks that fit a high word at the threshold K_of(t)."""
    spare = m._spare_counts()
    return [t for t in range(m.iters) if m._high_fits(t, K_of(t), spare)]


def plans(q, rng):
    probe = DialogMultiplier(Circuit(True), Mod(q))
    T, W = probe.iters, probe.widths
    joint, spare = tick_joint(q, T), probe._spare_counts()
    single = [(t, t + 1, joint[t]) for t in eligible(probe, lambda t: joint[t])]
    longest = []
    t = 1
    while t < T:
        if probe._high_fits(t, joint[t], spare):
            K, a = joint[t], t
            while t < T and probe._high_fits(t, K, spare):
                t += 1
            longest.append((a, t, K))
        else:
            t += 1
    out = {"single": single, "longest": longest}
    for r in range(3):
        cut = []
        s = 1
        while s < T:
            e = s
            while e < min(T, s + rng.randint(1, 4)) and probe._high_fits(e, joint[s], spare):
                e += 1
            if e > s:
                cut.append((s, e, joint[s]))
            s = max(e, s + 1)
        out[f"random{r}"] = cut
    under = [(t, t + 1, joint[t] - 1) for t in range(1, T - 1)
             if joint[t] >= 2 and probe._high_fits(t, joint[t] - 1, spare)]
    return out, under


def m_part(primes):
    """The toy pass has few clean wires, so these builds take no copies of g0*label and let a
    tick run one high bit; the copies run in the step part and in the secp256k1 build."""
    ok = True
    rng = random.Random(2026)
    control_bad = {"under": 0, "no_update": 0}
    DialogMultiplier.HIGH_COPIES, DialogMultiplier.HIGH_MIN = 0, 1
    for q in primes:
        ys = sorted({0, 1, q - 1, 3 % q}) if q < 3000 else [1, q - 1]
        named, under = plans(q, rng)
        for name, plan in named.items():
            DialogMultiplier.HIGH_STRETCHES = plan
            bad = m_bad(q, ys, False) + m_bad(q, ys, True)
            ticks = sum(b - a for a, b, _ in plan)
            print(f"q={q} plan={name} stretches={len(plan)} ticks={ticks} bad={bad}", flush=True)
            ok &= bad == 0 and len(plan) > 0
        if under:
            # The pass refuses a K under the joint threshold, so the control lowers its bound too.
            DialogMultiplier.HIGH_STRETCHES = under
            saved_joint = point_add.tick_joint
            point_add.tick_joint = lambda q, T: [k - 1 for k in saved_joint(q, T)]
            bad = m_bad(q, ys, False)
            point_add.tick_joint = saved_joint
            control_bad["under"] += bad
            print(f"q={q} control: K one under the joint threshold, bad={bad}", flush=True)
        DialogMultiplier.HIGH_STRETCHES = named["longest"]
        saved = DialogMultiplier._shared_difference

        def no_update(self, U, V, g0, g1, top=None, kept=None, copies=(), high=None):
            saved(self, U, V, g0, g1, top, kept, copies, high)
            if high is not None:
                self.C.cx(g1, high[1])

        DialogMultiplier._shared_difference = no_update
        bad = m_bad(q, ys, False)
        DialogMultiplier._shared_difference = saved
        control_bad["no_update"] += bad
        print(f"q={q} control: label not updated, bad={bad}", flush=True)
        DialogMultiplier.HIGH_STRETCHES = None
    DialogMultiplier.HIGH_COPIES, DialogMultiplier.HIGH_MIN = 16, 2
    caught = sum(v > 0 for v in control_bad.values())
    print(f"multiplier: {'ok' if ok else 'FAILED'}; controls caught {caught} of {len(control_bad)} "
          f"{control_bad}", flush=True)
    return ok and caught == len(control_bad)


def plan_part(primes, samples):
    ok = True
    for q in primes + [SECP256K1_Q]:
        m = DialogMultiplier(Circuit(), Mod(q))
        plan = m.high_plan()
        T = tick_bound(q)
        joint = tick_joint(q, T)
        rng = random.Random(q)
        xs = range(1, q) if q < 5000 else [rng.randrange(1, q) for _ in range(samples)] + [1, 2, q - 1]
        ticks = labels_ok = 0
        for x in xs:
            states, cells = tick_pass(q, x, T)
            for s in plan:
                assert all(s.K >= joint[t] for t in range(s.start, s.stop))
                U, V, E = states[s.start]
                z = int(U >> s.K != 0)
                for t in range(s.start, s.stop):
                    U, V, E = states[t]
                    vhat = V ^ E
                    uh, vh = U >> s.K, vhat >> s.K
                    good = not (uh and vh) and (not uh or z) and (not vh or not z)
                    ticks += 1
                    labels_ok += good
                    add, _ = cells[t]
                    swap = add and states[t + 1][2]
                    z ^= swap
        print(f"plan q={'secp256k1' if q == SECP256K1_Q else q}: {len(plan)} stretches, "
              f"{ticks} stretch ticks, {ticks - labels_ok} wrong", flush=True)
        ok &= labels_ok == ticks
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-width", type=int, default=6)
    ap.add_argument("--primes", type=int, nargs="*", default=[251, 509, 1021, 2039, 4093])
    ap.add_argument("--samples", type=int, default=3000)
    ap.add_argument("--skip", nargs="*", default=[])
    a = ap.parse_args()
    ok = True
    if "step" not in a.skip:
        ok &= step_part(a.max_width)
    if "m" not in a.skip:
        ok &= m_part(a.primes)
    if "plan" not in a.skip:
        ok &= plan_part(a.primes, a.samples)
    print("high_check:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
