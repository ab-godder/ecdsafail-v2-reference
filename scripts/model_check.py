"""Classical and bit-sliced gate checks of the dialog multiplier and pass.
Small primes cover every x; y is exhaustive below 300 and sampled above."""
import io
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod, load_const
from qec.builder import Circuit
from qec.emit import write_circuit
from qec.model import tick_bound, tick_mul, tick_pass, tick_sweep, tick_tape, tick_widths
from qec.point_add import DialogMultiplier


def gate_mul(q: int, ys, inverse: bool = False, iters: int | None = None) -> int:
    """Check every x against the supplied y values."""
    return gate_pairs(q, [(x, y) for x in range(1, q) for y in ys], inverse, iters)


def gate_pairs(q: int, pairs, inverse: bool = False, iters: int | None = None) -> int:
    """Count bad output or dirty-scratch lanes; malformed builds raise."""
    pairs = list(pairs)
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    m = DialogMultiplier(C, mod)
    if iters is not None:
        m.iters = iters
    m.apply(X, Y, inverse=inverse)
    assert len(set(X + Y)) == 2 * mod.n, "aliased output register"
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y})
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    xs, values = zip(*pairs)
    init = sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1])
    state = sim.run(ops, C.num_wires, init, len(pairs))
    expected = [y if x == 0 else y * (pow(x, -1, q) if inverse else x) % q
                for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def b1_swept(cells):
    """In Bezout order, the reverse of the tape, rotate b left by one in each run of no-subtract
    cells, so each B1 moves to its run's end."""
    out = list(reversed(cells))
    t = 0
    while t < len(out):
        end = t
        while end < len(out) and not out[end][0]:
            end += 1
        if end > t:
            bs = [b for _, b in out[t:end]]
            bs = bs[1:] + bs[:1]
            out[t:end] = [(0, b) for b in bs]
        t = end + 1
    return list(reversed(out))


def _pass_bad(q, xs, iters):
    C = Circuit()
    mod = Mod(q)
    m = DialogMultiplier(C, mod)
    if iters is not None:
        m.iters = iters
    U, V = C.alloc(mod.n), C.alloc(mod.n)
    inputs = list(V)
    pool = V + C.alloc(m.tape_external)
    load_const(C, U, q)
    held = C.alloc(m.iters)
    m._euclid_pass(U, pool, held)
    G = m.coded_tape(pool)
    for j in range(m.blocks):
        m._decode(G, j, C.alloc(len(m._block_free(j, m.swept))), m.swept)
    freed = sum(len(m._block_free(j, m.swept)) for j in range(m.blocks))
    assert len(set(U + pool + G)) == len(set(U + pool)) + freed, "aliased pass register"
    structural(C.prims, C.num_wires)
    state = sim.run(C.prims, C.num_wires, sim.pack(xs, inputs), len(xs))
    # Keep the reference model independent of tick_pass and tick_sweep.
    cells = [tick_tape(q, x, m.iters)[2] for x in xs]
    if m.swept:
        cells = [b1_swept(c) for c in cells]
    tapes = [sum((add + 2 * bit) << (2 * i) for i, (add, bit) in enumerate(c)) for c in cells]
    want = sim.pack([q if x == 0 else 1 for x in xs], U) | sim.pack(tapes, G)
    want |= sim.pack([held_offsets(tick_tape(q, x, m.iters)[2]) for x in xs], held)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def held_offsets(tape):
    """Bit i is a subtract cell's preceding offset: the swap of the last subtract before it."""
    out, e = 0, 0
    for i, (add, bit) in enumerate(tape):
        if add:
            out |= e << i
            e = bit
    return out


def gate_pass(q: int, iters: int | None = None) -> int:
    """Check terminal words, exact tape and all scratch on remapped inputs."""
    return _pass_bad(q, [x or 1 for x in range(q)], iters)


def gate_pass_zero(q: int) -> int:
    """Mismatched bits for a V = 0 caller against the old generic result: q kept, idle tape,
    clean scratch. Late tape takes U's vacated wires, so that caller is outside the pass's domain
    and this must be nonzero; the multiplier never makes the call, by the remap control."""
    return _pass_bad(q, [0], None)


def width_edge_inputs(q: int) -> list[tuple[int, int]]:
    """Schedule-edge powers and neighbours, boundaries, and seed-2026 samples."""
    rng = random.Random(2026)
    xs = {x for k in range(q.bit_length())
          for x in (2**k, 2**k - 1, 2**k + 1) if 0 <= x < q}
    xs |= {0, 1, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2}
    xs |= {rng.randrange(q) for _ in range(64)}
    return [(x, rng.choice([0, 1, q - 1, rng.randrange(q)])) for x in sorted(xs)]


def offset_checks(q, xs):
    """Compare held states and swept cells with the independent semantic tape."""
    T = tick_bound(q)
    offset_ok = tape_ok = True
    for x in xs:
        states, cells = tick_pass(q, x, T)
        if len(states) != T + 1 or len(cells) != T:
            return False, False
        u, v = q, x
        for i, (U, V, E) in enumerate(states):
            offset_ok &= U == u and V + E == v and V >= 0
            if i < T:
                if v & 1:
                    if u > v:
                        u, v = v, u
                    v = (v - u) // 2
                else:
                    v >>= 2 if v % 4 == 0 else 1
        swept, final = tick_sweep(cells, states[-1][2])
        tape_ok &= swept == tick_tape(q, x, T)[2] and states[-1][1:] == (0, 0) and final == 0
    return offset_ok, tape_ok


def trace(q, x):
    u, v, sizes = q, x, []
    while v:
        sizes.append(max(u, v).bit_length())
        if v & 1:
            if u > v:
                u, v = v, u
            v = (v - u) // 2
        else:
            v >>= 2 if v % 4 == 0 else 1
    return sizes


rows = []
for q in [3, 5, 7, 11, 13, 31, 43, 61, 127, 251, 509, 1021, 4093]:
    B, W = tick_bound(q), tick_widths(q, tick_bound(q))
    worst = 0
    width_ok = True
    for x in range(q):
        t = trace(q, x)
        worst = max(worst, len(t))
        width_ok &= all(s <= W[i] for i, s in enumerate(t))
    ys = range(q) if q < 300 else random.Random(q).sample(range(q), 8)
    mul_ok = all(tick_mul(q, x, y) == (x * y % q if x else y) for x in range(q) for y in ys)
    pass_end = all(tick_tape(q, x or 1, B)[:2] == (1, 0) for x in range(q))
    remap_control = [x for x in range(q) if tick_tape(q, x, B)[:2] != (1, 0)] == [0]
    short_rejected = any(tick_tape(q, x or 1, worst - 1)[1] != 0 for x in range(q))
    short_rejected &= gate_mul(q, ys, iters=worst - 1) > 0
    pass_gates = gate_pass(q) == 0
    pass_zero_outside = gate_pass_zero(q) > 0
    forward_bad = gate_mul(q, ys)
    inverse_bad = gate_mul(q, ys, inverse=True)
    pass_offset, pass_tape = offset_checks(q, range(1, q))
    rows.append({"q": q, "worst_iterations": worst, "bound": B,
                 "width_schedule_ok": width_ok, "multiplier_ok": mul_ok,
                 "pass_end_ok": pass_end, "remap_control_ok": remap_control,
                 "pass_short_rejected": short_rejected, "pass_gates_ok": pass_gates,
                 "multiplier_gates_ok": forward_bad == inverse_bad == 0,
                 "forward_bad": forward_bad, "inverse_bad": inverse_bad,
                 "pass_offset_ok": pass_offset, "pass_tape_ok": pass_tape, "pass_zero_outside_domain": pass_zero_outside})
    print(json.dumps(rows[-1]))
q = 2**256 - 2**32 - 977
rng = random.Random(1)
pairs = [(rng.randrange(q), rng.randrange(q)) for _ in range(2000)]
ok = all(tick_mul(q, x, y) == (x * y % q if x else y) for x, y in pairs)
pass_end = all(tick_tape(q, x or 1, tick_bound(q))[:2] == (1, 0)
               for x in [0, 1, q - 1] + [x for x, _ in pairs])
offset_ok = all(offset_checks(q, [1, q - 1] + [x for x, _ in pairs]))
edge_pairs = [(x, y) for x, y in width_edge_inputs(q) if x != 0]
width_edge_ok = gate_pairs(q, edge_pairs) == 0 and gate_pairs(q, edge_pairs, inverse=True) == 0
print(json.dumps({"q": "secp256k1", "bound": tick_bound(q), "random_pairs": 2000,
                  "multiplier_ok": ok, "pass_end_ok": pass_end, "offset_ok": offset_ok,
                  "width_edge_ok": width_edge_ok, "width_edge_lanes": len(edge_pairs)}))
keys = ("width_schedule_ok", "multiplier_ok", "pass_end_ok", "remap_control_ok",
        "pass_short_rejected", "pass_gates_ok", "multiplier_gates_ok",
        "pass_offset_ok", "pass_tape_ok", "pass_zero_outside_domain")
sys.exit(0 if ok and pass_end and offset_ok and width_edge_ok and all(r["worst_iterations"] <= r["bound"]
         and all(r[k] for k in keys) for r in rows) else 1)
