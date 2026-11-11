"""Depth pass over the hierarchical emission: reorder commuting ops inside each definition, then choose
each CCX's control order. Exact for every input; it changes no op, wire, call record or count.

Board model. Ops run as soon as possible per wire in program order. A CCX is the 15-gate Toffoli: its
target is busy at ticks 1-10, control2 at 2-11 and control1 at 4-11; X and CX take 1 tick, SWAP 3.

Commutation. Every primitive is X-type, t ^= f(controls). Per wire an op is 'c' (read only), 't' (written
only as an X-type target, never read) or 'e' (SWAP wires, and call wires a callee both reads and writes).
Two ops commute when every wire they share is c/c or t/t. Keeping each wire's maximal runs of compatible
roles in order therefore keeps the permutation; check_order asserts that for every order the pass keeps.

Per definition, in emission order (callees first): canonical control order (control1 < control2), then
an earliest-start list schedule of the run DAG forward and one backward, keeping the input unless a
schedule has fewer standalone forward plus inverse ticks; then each CCX takes the control order its own
standalone forward and inverse runs prefer. Decisions read definition bodies only, so definitions stay
identical over addends, and a second application returns its input.
"""
import hashlib
import heapq
from collections import Counter

from .builder import CCX, CX, SWAP, X, Call

# Local op encoding: (0, t, t, t) X, (1, c, t, t) CX, (2, a, b, b) SWAP, (3, c1, c2, t) CCX,
# (7, callee, wires, inverse) call.
_ENC = {X: 0, CX: 1, SWAP: 2, CCX: 3}
_DEC = {v: k for k, v in _ENC.items()}
_DUR = {0: 1, 1: 1, 2: 3, 3: 11}
_cache = {}


def _ccx_times(ra, rb, rt):
    s2 = max(rb, rt + 1) + 1
    s4 = max(ra, s2 + 1) + 1
    return s4 + 6, s4 + 7


def _canon(ops):
    return [(3, min(a, b), max(a, b), c) if k == 3 else (k, a, b, c) for k, a, b, c in ops]


def _op_roles(op, summary):
    k, a, b, c = op
    if k == 0:
        return ((a, "t"),)
    if k == 1:
        return ((a, "c"), (b, "t"))
    if k == 3:
        return ((a, "c"), (b, "c"), (c, "t"))
    if k == 2:
        return ((a, "e"), (b, "e"))
    rs = summary[a]
    return tuple((b[lw], rs[lw]) for lw in sorted(rs))


def _summaries(defs):
    """Per definition, its role on each local wire over the whole body, callees included."""
    out = []
    for ops in defs:
        acc = {}
        for op in ops:
            for w, r in _op_roles(op, out):
                p = acc.get(w)
                acc[w] = r if p is None or p == r else "e"
        out.append(acc)
    return out


def _dag(ops, summary):
    runs, where = {}, []
    for i, op in enumerate(ops):
        mine = []
        for w, r in _op_roles(op, summary):
            lst = runs.setdefault(w, [])
            if lst and r != "e" and lst[-1][0] == r:
                lst[-1][1].append(i)
            else:
                lst.append([r, [i]])
            mine.append((w, len(lst) - 1))
        where.append(mine)
    return where, runs


def check_order(ops, order, summary):
    """True when every wire's run index is non-decreasing along order, so every non-commuting pair keeps its order."""
    where, _ = _dag(ops, summary)
    if sorted(order) != list(range(len(ops))):
        return False
    last = {}
    for i in order:
        for w, g in where[i]:
            if g < last.get(w, 0):
                return False
            last[w] = g
    return True


def _run(defs, d, r, inv):
    for k, a, b, c in (reversed(defs[d]) if inv else defs[d]):
        if k == 7:
            sub = [r[w] for w in b]
            _run(defs, a, sub, inv ^ c)
            for j, w in enumerate(b):
                r[w] = sub[j]
        elif k == 3:
            t, cc = _ccx_times(r[a], r[b], r[c])
            r[c] = t
            r[a] = r[b] = cc
        elif k == 0:
            r[a] += 1
        else:
            v = max(r[a], r[b]) + _DUR[k]
            r[a] = r[b] = v


def _standalone(defs, arity, d, inv):
    r = [0] * arity[d]
    _run(defs, d, r, inv)
    return max(r, default=0)


def _schedule(ops, summary, defs, dur_call):
    """Earliest-start list schedule; ties by bottom level, then by op content and run membership."""
    n = len(ops)
    where, runs = _dag(ops, summary)
    pending = [sum(1 for _, g in where[i] if g > 0) for i in range(n)]
    left = {w: [len(m) for _, m in lst] for w, lst in runs.items()}
    bl, best = [0] * n, {}
    for i in range(n - 1, -1, -1):
        k = ops[i][0]
        d = dur_call(ops[i][1]) if k == 7 else _DUR[k]
        bl[i] = d + max((best.get((w, g + 1), 0) for w, g in where[i]), default=0)
        for w, g in where[i]:
            if bl[i] > best.get((w, g), 0):
                best[(w, g)] = bl[i]
    canon = [(op[0], op[1], op[2] if op[0] != 7 else 0, op[3] if op[0] != 7 else 0, tuple(where[i]))
             for i, op in enumerate(ops)]
    width = 1 + max((max(op[2]) if op[0] == 7 else max(op[1:]) for op in ops), default=0)
    r = [0] * width

    def est(i):
        k, a, b, c = ops[i]
        if k == 3:
            return _ccx_times(r[a], r[b], r[c])[1] - 11
        if k == 0:
            return r[a]
        if k == 7:
            return max((r[w] for w in b), default=0)
        return max(r[a], r[b])

    heap = [((est(i), -bl[i], canon[i]), i) for i in range(n) if pending[i] == 0]
    heapq.heapify(heap)
    order = []
    while heap:
        p, i = heapq.heappop(heap)
        q = (est(i), -bl[i], canon[i])
        if q != p:
            heapq.heappush(heap, (q, i))
            continue
        order.append(i)
        k, a, b, c = ops[i]
        if k == 3:
            t, cc = _ccx_times(r[a], r[b], r[c])
            r[c] = t
            r[a] = r[b] = cc
        elif k == 0:
            r[a] += 1
        elif k == 7:
            sub = [r[w] for w in b]
            _run(defs, a, sub, c)
            for j, w in enumerate(b):
                r[w] = sub[j]
        else:
            v = max(r[a], r[b]) + _DUR[k]
            r[a] = r[b] = v
        for w, g in where[i]:
            cnt = left[w]
            cnt[g] -= 1
            if cnt[g] == 0 and g + 1 < len(cnt):
                for j in runs[w][g + 1][1]:
                    pending[j] -= 1
                    if pending[j] == 0:
                        heapq.heappush(heap, ((est(j), -bl[j], canon[j]), j))
    assert len(order) == n
    return order


def _reorder(defs, arity):
    defs = [_canon(d) for d in defs]
    summary = _summaries(defs)
    st = {}

    def dur_call(d):
        if d not in st:
            st[d] = _standalone(defs, arity, d, False)
        return st[d]

    for d in range(len(defs)):
        ops, L = defs[d], len(defs[d])
        back = _schedule(ops[::-1], summary, defs, dur_call)
        cands = [list(range(L)), _schedule(ops, summary, defs, dur_call), [L - 1 - j for j in reversed(back)]]
        best = None
        for order in cands:
            assert check_order(ops, order, summary)
            body = [ops[i] for i in order]
            defs[d] = body
            s = _standalone(defs, arity, d, False) + _standalone(defs, arity, d, True)
            if best is None or s < best[0]:
                best = (s, body)
        defs[d] = best[1]
        st.pop(d, None)
    return defs


def _orient(defs, arity):
    out = [list(d) for d in defs]
    for d in range(len(defs)):
        vote = Counter()
        for inv in (False, True):
            r = [0] * arity[d]

            def run(ops, m, inv2, top):
                for i in (range(len(ops) - 1, -1, -1) if inv2 else range(len(ops))):
                    k, a, b, c = ops[i]
                    if k == 7:
                        run(out[a], [m[w] for w in b] if m is not None else list(b), inv2 ^ c, False)
                        continue
                    if m is not None:
                        a, b, c = m[a], m[b], m[c]
                    if k == 3:
                        if top and r[a] != r[b]:
                            vote[i] += 1 if r[b] > r[a] else -1
                        t, cc = _ccx_times(r[a], r[b], r[c])
                        r[c] = t
                        r[a] = r[b] = cc
                    elif k == 0:
                        r[a] += 1
                    else:
                        v = max(r[a], r[b]) + _DUR[k]
                        r[a] = r[b] = v

            run(defs[d], None, inv, True)
        for i, v in vote.items():
            if v > 0:
                k, a, b, c = out[d][i]
                out[d][i] = (k, b, a, c)
    return out


def _local(block, ids):
    """A block's body in its local numbering (sorted wires), calls inline as (7, callee, wires, inverse)."""
    local = {w: k for k, w in enumerate(sorted(block.wires))}
    ops, body, pos = block.ops, [], 0

    def prims(lo, hi):
        for p in range(4 * lo, 4 * hi, 4):
            op, a, b, c = ops[p:p + 4]
            body.append((_ENC[op], local[a], local[b], local[c]))

    for call in block.calls:
        prims(pos, call.pos)
        body.append((7, ids[call.block], tuple(local[w] for w in call.wires), call.inverse))
        pos = call.pos
    prims(pos, len(ops) // 4)
    return body


def _restore(block, body, order_blocks):
    """Back to the builder's wires: a primitive op array and Calls with recomputed positions."""
    from array import array
    wires = sorted(block.wires)
    ops, calls = array("i"), []
    for k, a, b, c in body:
        if k == 7:
            calls.append(Call(len(ops) // 4, order_blocks[a], [wires[w] for w in b], c))
        else:
            ops.extend((_DEC[k], wires[a], wires[b], wires[c]))
    return ops, calls


def optimize(order, ids):
    """{block: (ops, calls)} for every definition, in the emission's numbering."""
    bodies = [_local(b, ids) for b in order]
    arity = [len(b.wires) for b in order]
    key = hashlib.sha256(repr((bodies, arity)).encode()).hexdigest()
    if key not in _cache:
        _cache.clear()
        _cache[key] = _orient(_reorder(bodies, arity), arity)
    return {b: _restore(b, body, order) for b, body in zip(order, _cache[key])}


def canonical_top(prims):
    """Top-level primitives with control1 < control2 on every CCX; nothing else changes."""
    from array import array
    out = array("i", prims)
    for p in range(0, len(out), 4):
        if out[p] == CCX and out[p + 1] > out[p + 2]:
            out[p + 1], out[p + 2] = out[p + 2], out[p + 1]
    return out
