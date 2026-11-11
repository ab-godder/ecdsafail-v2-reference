"""Reorder commuting gates inside definition bodies to cut board ticks. Same unitary.

    python3 -I qec/reorder_pass.py <in.json> <out.json> [--jobs N] [--max-ops N] [--defs LIST]

<in.json> is a circuit or an upload with the circuit under "circuit"; the output
keeps that shape. scripts/post_pass.py runs it after qec/split_dirs.py and
qec/gate_pass.py.

For each definition body that holds only gates (no calls), a local search changes
the order of the body's records and the control order of its CCX records:

- Two gates of the form "X on a target under controls" (X, CX, CCX) commute when
  neither's target is among the other's controls. A move takes a commuting edge
  p -> g of the body's critical path (p before g, sharing a wire) and moves g to just
  before p, or p to just after g, but only when every gate in between that shares a
  wire with the moved gate commutes with it. So the new order is a reordering of
  commuting gates and the body's unitary is unchanged.
- A CCX on the critical path may swap control1 and control2; CCX is symmetric in its
  controls, so that is the same gate.

The search objective is the body's depth under the board tick model with all inputs
ready at once, summed over the directions (forward, inverted) the body runs in,
tie-broken by the sum of output ready times. A move is kept only if it lowers the
objective. The board scores the expanded circuit, where call bodies overlap, so the
pass then puts the new bodies in, one group of bodies with the same callers at a
time, and keeps a group only if the whole circuit's ticks do not rise; a group that
fails is retried one body at a time.
"""
import json
import sys
import time

INF = 10 ** 9


def ccx(ra, rb, rt):
    """Board lowering of CCX(a, b -> t): new ready times of a, b, t."""
    rt += 1
    x = max(rb, rt) + 1; rb = rt = x; rt += 1
    x = max(ra, rt) + 1; ra = rt = x; rt += 1
    x = max(rb, rt) + 1; rb = rt = x; rt += 1
    x = max(ra, rt) + 1; ra = rt = x
    rb += 1; rt += 2
    x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
    x = max(ra, rb) + 1; ra = rb = x
    return ra, rb, rt


LCCX = [list(ccx(*[0 if j == i else -INF for j in range(3)])) for i in range(3)]


# ---------------------------------------------------------------- body model
# A compiled gate is [kind, wires, target, controls]; kind 0 X, 1 CX, 2 SWAP, 3 CCX.
# wires of a CCX are (control1, control2, target).

def compile_ops(ops):
    out = []
    for r in ops:
        o = r["op"]
        if o == "X":
            out.append([0, (r["target"],), r["target"], ()])
        elif o == "CX":
            out.append([1, (r["control"], r["target"]), r["target"], (r["control"],)])
        elif o == "CCX":
            out.append([3, (r["control1"], r["control2"], r["target"]), r["target"],
                        (r["control1"], r["control2"])])
        elif o == "SWAP":
            out.append([2, (r["wire1"], r["wire2"]), None, None])
        else:
            raise ValueError(o)
    return out


def commute(a, b):
    if a[2] is None or b[2] is None:
        return False
    return a[2] not in b[3] and b[2] not in a[3]


SNAP = 64


class Run:
    """One direction's schedule of a body: per position the gate's wires and output ready
    times, the ready vector every SNAP positions, and what the critical-path trace needs."""

    def __init__(self, cops, seq, arity, rin=None):
        ready = list(rin) if rin else [0] * arity
        last = [-1] * arity
        self.seq = seq
        self.outs = outs = []
        self.snaps = snaps = []
        self.rec = rec = []
        for k, gi in enumerate(seq):
            if k % SNAP == 0:
                snaps.append(ready[:])
            g = cops[gi]
            ws = g[1]
            rs = tuple(ready[w] for w in ws)
            rec.append((gi, rs, tuple(last[w] for w in ws)))
            vals = gate_out(g[0], rs)
            for w, v in zip(ws, vals):
                ready[w] = v
                last[w] = k
            outs.append((ws, vals))
        self.ready, self.last = ready, last
        self.kinds = [cops[gi][0] for gi in seq]
        self.value = (max(ready), sum(ready))

    def path(self, cops):
        """Critical path as (pred_gate or -1, gate), last gate first."""
        ready, rec = self.ready, self.rec
        w = max(range(len(ready)), key=ready.__getitem__)
        at = self.last[w]
        if at < 0:
            return []
        pos = cops[rec[at][0]][1].index(w)
        path = []
        while at >= 0:
            gi, rs, preds = rec[at]
            if cops[gi][0] == 3:
                ip = max(range(3), key=lambda i: rs[i] + LCCX[i][pos])
            else:
                ip = max(range(len(rs)), key=rs.__getitem__)
            p = preds[ip]
            path.append((rec[p][0] if p >= 0 else -1, gi))
            if p < 0:
                break
            pos = cops[rec[p][0]][1].index(cops[gi][1][ip])
            at = p
        return path

    def trial(self, cops, seq, lo, hi, exact=False):
        """(max, sum) of the final ready times for seq, which equals self.seq outside
        positions lo..hi; reuses this run's times wherever they cannot have changed.
        Returns (value, is_exact). Unless exact, it stops as soon as no wire is ahead of
        this run and some wire is behind: the schedule is monotone in ready times, so
        every final time is then at least this run's, and value is this run's (a lower
        bound)."""
        outs = self.outs
        k0 = lo - lo % SNAP
        ready = self.snaps[k0 // SNAP][:]
        for k in range(k0, lo):
            ws, vals = outs[k]
            for w, v in zip(ws, vals):
                ready[w] = v
        for k in range(lo, hi + 1):
            g = cops[seq[k]]
            ws = g[1]
            vals = gate_out(g[0], tuple(ready[w] for w in ws))
            for w, v in zip(ws, vals):
                ready[w] = v
        old = {}
        for k in range(hi, lo - 1, -1):
            ws, vals = outs[k]
            for w, v in zip(ws, vals):
                if w not in old:
                    old[w] = v
        ahead = {w for w, v in old.items() if ready[w] < v}
        behind = {w for w, v in old.items() if ready[w] > v}
        if not ahead and (not behind or not exact):
            return self.value, not behind
        kinds = self.kinds
        for k in range(hi + 1, len(outs)):
            ws, vals = outs[k]
            if ahead.isdisjoint(ws) and behind.isdisjoint(ws):
                for w, v in zip(ws, vals):
                    ready[w] = v
                continue
            new = gate_out(kinds[k], tuple(ready[w] for w in ws))
            for w, v, o in zip(ws, new, vals):
                ready[w] = v
                if v < o:
                    ahead.add(w); behind.discard(w)
                elif v > o:
                    behind.add(w); ahead.discard(w)
                else:
                    ahead.discard(w); behind.discard(w)
            if not ahead and not exact:
                return self.value, not behind
        if not ahead and not behind:
            return self.value, True
        return (max(ready), sum(ready)), True


def gate_out(kind, rs):
    if kind == 3:
        return ccx(*rs)
    if kind == 0:
        return (rs[0] + 1,)
    x = max(rs) + (3 if kind == 2 else 1)
    return (x, x)


def runs(cops, order, arity, scenarios):
    """One Run per scenario (inverted, input ready times or None for all at 0)."""
    return [(inv, Run(cops, order[::-1] if inv else order, arity, rin)) for inv, rin in scenarios]


def total(rs):
    return tuple(map(sum, zip(*(r.value for _, r in rs))))


def trial(cops, rs, new, i, j, best):
    """Objective of order new, which differs from the runs' order only at positions i..j,
    or None when it is certainly not below best."""
    n = len(new)
    rev = None
    vals = []
    for inv, r in rs:
        if inv:
            rev = rev or new[::-1]
            vals.append(r.trial(cops, rev, n - 1 - j, n - 1 - i) + (rev, n - 1 - j, n - 1 - i))
        else:
            vals.append(r.trial(cops, new, i, j) + (new, i, j))
    m = sum(v[0][0] for v in vals); s = sum(v[0][1] for v in vals)
    if (m, s) >= best:
        return None
    if all(v[1] for v in vals):
        return m, s
    m = s = 0
    for (inv, r), (val, ok, seq, lo, hi) in zip(rs, vals):
        if not ok:
            val, _ = r.trial(cops, seq, lo, hi, exact=True)
        m += val[0]; s += val[1]
    return (m, s) if (m, s) < best else None


def moved(cops, order, i, j, which):
    """Order with g=order[j] moved before p=order[i] ('g') or p moved after g ('p'), or None
    when an in-between gate on a shared wire does not commute with the moved gate."""
    if which == "g":
        g = cops[order[j]]; ws = set(g[1])
        for q in range(i, j):
            o = cops[order[q]]
            if not ws.isdisjoint(o[1]) and not commute(o, g):
                return None
        return order[:i] + [order[j]] + order[i:j] + order[j + 1:]
    p = cops[order[i]]; ws = set(p[1])
    for q in range(i + 1, j + 1):
        o = cops[order[q]]
        if not ws.isdisjoint(o[1]) and not commute(o, p):
            return None
    return order[:i] + order[i + 1:j + 1] + [order[i]] + order[j + 1:]


def flip(g):
    a, b, t = g[1]
    g[1] = (b, a, t)
    g[3] = (b, a)


def search(cops, arity, scenarios, max_moves=100000, order=None):
    """Local search over orders and CCX control orders, from order (default: as given).
    Returns (order, start, end, moves); CCX flips are made in cops.

    Each sweep takes the critical path of every scenario once and tries its moves in
    turn, keeping each that lowers the objective right away; sweeps repeat until one
    keeps nothing."""
    order = list(range(len(cops))) if order is None else order
    rs = runs(cops, order, arity, scenarios)
    best = start = total(rs)
    moves = 0
    while moves < max_moves:
        kept = 0
        for p, g in [e for _, run in rs for e in run.path(cops)]:
            posof = None
            if cops[g][0] == 3:
                posof = {x: k for k, x in enumerate(order)}
                flip(cops[g])
                q = posof[g]
                ob = trial(cops, rs, order, q, q, best)
                if ob is not None:
                    best = ob
                    rs = runs(cops, order, arity, scenarios)
                    kept += 1
                    continue
                flip(cops[g])
            if p < 0 or not commute(cops[p], cops[g]):
                continue
            if posof is None:
                posof = {x: k for k, x in enumerate(order)}
            i, j = sorted((posof[p], posof[g]))
            for which in ("g", "p"):
                new = moved(cops, order, i, j, which)
                if new is None:
                    continue
                ob = trial(cops, rs, new, i, j, best)
                if ob is not None:
                    best, order = ob, new
                    rs = runs(cops, order, arity, scenarios)
                    kept += 1
                    break
        assert total(rs) == best, (total(rs), best)
        moves += kept
        if not kept:
            break
    return order, start, best, moves


def kicked(cops, order, arity, scenarios, rng, kicks):
    """order after kicks random legal moves on commuting edges of the critical path,
    whether or not they help."""
    for _ in range(kicks):
        run = Run(cops, order[::-1] if scenarios[0][0] else order, arity, scenarios[0][1])
        edges = [(p, g) for p, g in run.path(cops) if p >= 0 and commute(cops[p], cops[g])]
        if not edges:
            break
        p, g = rng.choice(edges)
        posof = {x: k for k, x in enumerate(order)}
        i, j = sorted((posof[p], posof[g]))
        new = moved(cops, order, i, j, rng.choice("gp"))
        if new is not None:
            order = new
    return order


def search_body(args):
    """Local search, then rounds of iterated local search: kick the best order with a few
    random legal moves, search again, keep the result if it is better. The random
    generator is seeded by the definition index, so the result is reproducible."""
    d, ops, arity, scenarios, rounds = args
    cops = compile_ops(ops)
    t0 = time.time()
    order, start, end, moves = search(cops, arity, scenarios)
    import random
    rng = random.Random(d)
    for _ in range(rounds):
        trial_cops = [list(g) for g in cops]
        o2 = kicked(trial_cops, order, arity, scenarios, rng, 4)
        o2, _, e2, n2 = search(trial_cops, arity, scenarios, order=o2)
        if e2 < end:
            order, end, cops, moves = o2, e2, trial_cops, moves + n2
    new = []
    for gi in order:
        r = ops[gi]
        if r["op"] == "CCX" and cops[gi][1][0] != r["control1"]:
            r = dict(r, control1=r["control2"], control2=r["control1"])
        new.append(r)
    return d, new, start, end, moves, time.time() - t0


# ---------------------------------------------------------------- circuit level

def schedule(circuit):
    """Board ticks of the expanded circuit."""
    bodies = [d["operations"] for d in circuit.get("definitions", [])]
    ready = [0] * circuit["num_wires"]

    def run(ops, wmap, inverse):
        for r in (reversed(ops) if inverse else ops):
            if "call" in r:
                run(bodies[r["call"]], [wmap[w] for w in r["wires"]], inverse ^ bool(r.get("inverse", False)))
                continue
            o = r["op"]
            if o == "X":
                ready[wmap[r["target"]]] += 1
            elif o == "CX":
                a, t = wmap[r["control"]], wmap[r["target"]]
                ready[a] = ready[t] = max(ready[a], ready[t]) + 1
            elif o == "CCX":
                a, b, t = wmap[r["control1"]], wmap[r["control2"]], wmap[r["target"]]
                ready[a], ready[b], ready[t] = ccx(ready[a], ready[b], ready[t])
            elif o == "SWAP":
                a, b = wmap[r["wire1"]], wmap[r["wire2"]]
                ready[a] = ready[b] = max(ready[a], ready[b]) + 3
            else:
                raise SystemExit(f"op {o} not handled")

    run(circuit["operations"], list(range(circuit["num_wires"])), False)
    return max(ready)


def directions(circuit):
    defs = circuit["definitions"]
    seen = {}

    def walk(ops, inverse):
        for r in ops:
            if "call" in r:
                d, inv = r["call"], inverse ^ bool(r.get("inverse", False))
                s = seen.setdefault(d, set())
                if inv not in s:
                    s.add(inv)
                    walk(defs[d]["operations"], inv)

    walk(circuit["operations"], False)
    return {d: tuple(sorted(s)) for d, s in seen.items()}


def run_counts(circuit):
    """Executions of each definition in the expanded circuit."""
    defs = circuit["definitions"]
    runs = [0] * len(defs)
    for r in circuit["operations"]:
        if "call" in r:
            runs[r["call"]] += 1
    for d in range(len(defs) - 1, -1, -1):
        for r in defs[d]["operations"]:
            if "call" in r:
                runs[r["call"]] += runs[d]
    return runs


def apply(circuit, jobs=4, min_ops=30, max_ops=50000, only=None, log=True, rounds=0, rounds_min_runs=50):
    """Run the pass on circuit in place. Bodies that run at least rounds_min_runs times
    also get rounds of iterated local search. Returns {"ticks_in", "ticks_out"}."""
    say = print if log else (lambda *a, **k: None)
    defs = circuit["definitions"]
    dirs = directions(circuit)
    nrun = run_counts(circuit)
    work = [(d, defs[d]["operations"], defs[d]["arity"], [(inv, None) for inv in dirs[d]],
             rounds if nrun[d] >= rounds_min_runs else 0) for d in sorted(dirs)
            if min_ops <= len(defs[d]["operations"]) <= max_ops
            and not any("call" in r for r in defs[d]["operations"])
            and (only is None or d in only)]
    work.sort(key=lambda w: -len(w[1]) * len(w[3]) * (1 + w[4]))
    base = schedule(circuit)
    say(json.dumps({"ticks_in": base, "bodies": len(work)}), flush=True)
    found = {}
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(jobs) as pool:
            results = pool.imap_unordered(search_body, work)
            for d, new, s, e, n, dt in results:
                found[d] = new
                say(f"  def {d}: {s} -> {e} in {n} moves, {dt:.0f}s", flush=True)
    else:
        for w in work:
            d, new, s, e, n, dt = search_body(w)
            found[d] = new
            say(f"  def {d}: {s} -> {e} in {n} moves, {dt:.0f}s", flush=True)
    # Accept against the whole circuit's ticks: bodies with the same callers as one group,
    # and a group that makes the circuit slower one body at a time.
    callers = {}
    for c, body in enumerate(defs):
        for r in body["operations"]:
            if "call" in r:
                callers.setdefault(r["call"], set()).add(c)
    groups = {}
    for d in sorted(found):
        groups.setdefault(tuple(sorted(callers.get(d, ()))), []).append(d)
    best = base

    def attempt(grp):
        nonlocal best
        old = {d: defs[d]["operations"] for d in grp}
        for d in grp:
            defs[d]["operations"] = found[d]
        t = schedule(circuit)
        keep = t <= best
        say(f"  defs {grp[:4]}{'...' if len(grp) > 4 else ''} ({len(grp)}): {t} {'kept' if keep else 'dropped'}",
            flush=True)
        if keep:
            best = t
        else:
            for d in grp:
                defs[d]["operations"] = old[d]
        return keep

    for grp in groups.values():
        if not attempt(grp) and len(grp) > 1:
            for d in grp:
                attempt([d])
    return {"ticks_in": base, "ticks_out": best}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--min-ops", type=int, default=30)
    ap.add_argument("--max-ops", type=int, default=50000)
    ap.add_argument("--defs", default="", help="only these definitions (comma list, a-b ranges)")
    ap.add_argument("--rounds", type=int, default=0, help="iterated local search rounds for hot bodies")
    a = ap.parse_args()
    only = None
    if a.defs:
        only = set()
        for part in a.defs.split(","):
            lo, _, hi = part.partition("-")
            only.update(range(int(lo), int(hi or lo) + 1))
    with open(a.src) as fh:
        doc = json.load(fh)
    circuit = doc.get("circuit", doc)
    res = apply(circuit, a.jobs, a.min_ops, a.max_ops, only, rounds=a.rounds)
    with open(a.dst, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print(json.dumps(res))


if __name__ == "__main__":
    main()
