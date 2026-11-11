"""Exhaustive check behind the secp256k1 width cut; see DESIGN.md, "A width cut at secp256k1".

A tick t that runs one bit under the generic schedule needs max(u, v)^2 4^t < 2^768 there.
W3 reduces any failure at a tick t <= SECP256K1_CUT_LAST to a state (u, 1) with
u >= 2^(384 - SECP256K1_CUT_LAST) reached through Euclid groups of total slack below 2 (W4, W5).
`search` enumerates every such group path backwards from (u, 1), W6, in exact integer
arithmetic. `scripts/width_search.c` is the same search, node for node, in fixed-width integers;
the script compiles it, checks that it matches `search` on small primes and on secp256k1 with a
higher floor on u, and runs it sharded at the cut's floor. Every candidate x is then simulated.
Before that, `search` with no floor on u is compared with brute force over every x on small
primes, and must find planted tight paths on wide moduli; those check the code, not the argument.

    python3 -I scripts/width_check.py            # the C search at the cut, on every core
    python3 -I scripts/width_check.py --python   # the Python search at the cut instead (slow)
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec.model import SECP256K1_CUT_LAST, SECP256K1_Q, tick_bound, tick_widths

P = 32                          # y = Y / 2^P
Y_LO = (17 << P) // 100         # 0.17 < 3 - 2 sqrt 2
Y_HI = -(-(583 << P) // 100)    # 5.83 > 3 + 2 sqrt 2
SHARD_DEPTH = 12
SMALL = (1021, 2039, 4093, 8191, 16381, 65519, 65521, 131063, 131071)


def theta(q: int) -> int:
    """2^768 at secp256k1: the cut width at tick t holds M exactly when M^2 4^t < theta."""
    return 1 << (2 * q * q * (q - 1)).bit_length() - 1


def step(u: int, v: int) -> tuple[int, int]:
    if v & 1:
        return min(u, v), abs(u - v) >> 1
    return u, v >> (2 if v & 3 == 0 else 1)


def breaks(q: int, x: int) -> list[int]:
    """Ticks t with v != 0 and max(u, v)^2 4^t >= theta, by direct simulation."""
    th, u, v, t, out = theta(q), q, x, 0, []
    while v:
        M = max(u, v)
        if M * M << 2 * t >= th:
            out.append(t)
        u, v = step(u, v)
        t += 1
    return out


def inputs(q: int, a1: int, a0: int, b1: int, b0: int, k: int) -> set[int]:
    """The x of a first group from (q, m) with q = a1 u + a0 4^k, u odd, m = b1 u + b0 4^k."""
    r = q - (a0 << 2 * k)
    if r <= 0 or r % a1 or not r // a1 & 1:
        return set()
    m = b1 * (r // a1) + (b0 << 2 * k)
    return {m << j0 for j0 in range(1, q.bit_length(), 2) if 2 * m < q and m << j0 < q}


def search(q: int, log_umin: int, shard: int = 0, shards: int = 1) -> tuple[set[int], int]:
    """A superset of the x with odd 2-adic valuation, not a power of two, whose trajectory
    reaches some (u, 1) with u >= 2^log_umin through groups of total slack below 2.
    With shards > 1, only the subtrees at depth SHARD_DEPTH numbered shard mod shards."""
    cands, nodes, seen = set(), 0, 0
    qP = q << P

    def terminal(a1, a0, b1, b0, IN0, IN1):
        # The first group starts at (q, m) and swaps: a1 u + a0 4^k = q with y = u / 4^k in
        # [Y0, Y1] / 2^P, that is IN0 4^k <= q 2^P <= IN1 4^k.
        k_lo = max(1, (qP.bit_length() - IN1.bit_length() - 1) // 2 - 1)
        k_hi = (qP.bit_length() - IN0.bit_length() + 1) // 2 + 1
        for k in range(k_lo, k_hi + 1):
            if IN0 << 2 * k <= qP <= IN1 << 2 * k:
                cands.update(inputs(q, a1, a0, b1, b0, k))

    def node(a1, a0, b1, b0, T, Y0, Y1, depth):
        # Start state D (a1 y + a0, b1 y + b0) on y in [Y0, Y1] / 2^P, after groups of T - 1
        # ticks besides the last; their slack is phi(y) / (4^T y^2) >= phi(y0) / (4^T y1^2).
        nonlocal nodes, seen
        if depth == SHARD_DEPTH:
            seen += 1
            if seen % shards != shard:
                return
        if depth >= SHARD_DEPTH or shard == 0:
            nodes += 1
        IA, IB = a1 * Y0 + (a0 << P), b1 * Y0 + (b0 << P)
        hi, lo = max(IA, IB), min(IA, IB)
        if hi * hi * lo >= Y1 * Y1 << 2 * T + 1 + P:
            return
        cap = q * Y1               # L(y0) > q y1 / (2 U_min), scaled by 2^P
        if hi << log_umin + 1 > cap:
            return
        if Y1 - Y0 > Y0 >> 3:      # split where part of the interval may prune; exact test
            Ym = (Y0 + Y1) >> 1
            Ma, Mb = a1 * Ym + (a0 << P), b1 * Ym + (b0 << P)
            if Y1 > 3 * Y0 or 10 * max(Ma, Mb) ** 2 * min(Ma, Mb) >= 13 * (Ym * Ym << 2 * T + P):
                node(a1, a0, b1, b0, T, Y0, Ym, depth + 1)
                node(a1, a0, b1, b0, T, Ym, Y1, depth + 1)
                return
        j = 0
        while True:
            f = 2 << j
            IN0 = IA + f * IB
            if IN0 << log_umin > cap:   # grows with j; past it even A = q is out of reach
                break
            T2 = T + 1 + (j + 1) // 2
            if IN0 * IN0 * IA < Y1 * Y1 << 2 * T2 + 1 + P:
                n1, n0 = a1 + f * b1, a0 + f * b0
                terminal(n1, n0, a1, a0, IN0, n1 * Y1 + (n0 << P))
                if IN0 << log_umin + 1 <= cap:   # only the first group's start reaches q / 2
                    node(n1, n0, a1, a0, T2, Y0, Y1, depth + 1)   # it swapped: (A + fB, A)
                    node(a1, a0, n1, n0, T2, Y0, Y1, depth + 1)   # it did not: (A, A + fB)
            j += 1

    # The last group starts at (u + D, u) or (u, u + D); its slack is (1 + y)^2 / (4 y).
    if shard == 0:
        terminal(1, 1, 1, 0, Y_LO + (1 << P), Y_HI + (1 << P))
    node(1, 1, 1, 0, 1, Y_LO, Y_HI, 0)
    node(1, 0, 1, 1, 1, Y_LO, Y_HI, 0)
    return cands, nodes


def c_search(exe: str, q: int, log_umin: int, shard: int = 0, shards: int = 1) -> tuple[set[int], int]:
    out = subprocess.run([exe, hex(q), str(log_umin), str(shard), str(shards)],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    cands, nodes = set(), None
    for line in out:
        if line.startswith("T "):
            cands |= inputs(q, *map(int, line.split()[1:]))
        elif line.startswith("N "):
            nodes = int(line.split()[1])
    assert nodes is not None, "width_search did not finish"
    return cands, nodes


def sharded(run, q: int, log_umin: int, shards: int, pool) -> tuple[set[int], int]:
    with pool(max_workers=min(shards, os.cpu_count() or 1)) as ex:
        parts = list(ex.map(run, [q] * shards, [log_umin] * shards, range(shards), [shards] * shards))
    return set().union(*(c for c, _ in parts)), sum(n for _, n in parts)


def brute_force(q: int) -> set[int]:
    """The x the search must return with no floor on u: odd valuation, not a power of two, and
    some (u, 1), u >= 3, with u^2 4^t > q^2 x, which is total group slack below 2 (W4)."""
    need = set()
    for x in range(1, q):
        if (x & -x).bit_length() % 2 or x & (x - 1) == 0:
            continue
        u, v, t = q, x, 0
        while v:
            if v == 1 and u >= 3 and u * u << 2 * t > q * q * x:
                need.add(x)
                break
            u, v = step(u, v)
            t += 1
    return need


def tight(q: int, x: int, log_umin: int) -> bool:
    """Whether x reaches some (u, 1), u >= max(3, 2^log_umin), with u^2 4^t > q^2 x."""
    u, v, t = q, x, 0
    while v:
        if v == 1 and u >= 3 and u >> log_umin and u * u << 2 * t > q * q * x:
            return True
        u, v = step(u, v)
        t += 1
    return False


def planted(count: int = 12, seed: int = 2026) -> list[tuple[int, int, int]]:
    """Wide (q, x, floor) cases with a known tight path, built backwards from (u, 1) with random
    swapping and non-swapping groups; q is any odd number, as W1 to W6 never use primality."""
    rng, out = random.Random(seed), []
    while len(out) < count:
        bits = rng.choice([64, 128, 256])
        k = rng.randrange(bits // 8, bits // 3)
        u = (4**k * rng.randrange(60, 170) // 100) | 1
        A, B = (u + 4**k, u) if rng.random() < 0.5 else (u, u + 4**k)
        for _ in range(rng.randrange(2, 12)):
            if rng.random() < 0.5 and B > A:
                A, B = A, A + (2 << rng.choice([0, 1])) * B
            else:
                A, B = A + (2 << max(0, (A // B).bit_length() - 1 + rng.choice([-1, 0, 1]))) * B, A
        q, m = A + (2 << max(0, (A // B).bit_length() - 1 + rng.choice([0, 1]))) * B, A
        if 2 * m < q and tight(q, 2 * m, u.bit_length() - 1):
            out.append((q, 2 * m, u.bit_length() - 1))
    return out


def main() -> int:
    use_python = "--python" in sys.argv[1:]
    ok = True
    exe = str(Path(tempfile.mkdtemp()) / "width_search")
    subprocess.run(["cc", "-O2", "-o", exe, str(ROOT / "scripts" / "width_search.c")], check=True)
    for q in SMALL:
        cands, nodes = search(q, 0)
        missing = sorted(brute_force(q) - cands)
        c_same = c_search(exe, q, 0) == (cands, nodes)
        ok &= not missing and c_same
        print(json.dumps({"q": q, "nodes": nodes, "candidates": len(cands), "missing": missing,
                          "c_matches": c_same}), flush=True)
    for q, x, log_umin in planted():
        cands, nodes = search(q, log_umin)
        c_same = c_search(exe, q, log_umin) == (cands, nodes)
        ok &= x in cands and c_same
        print(json.dumps({"planted_bits": q.bit_length(), "u_floor_log2": log_umin, "nodes": nodes,
                          "found": x in cands, "c_matches": c_same}), flush=True)
    q, n = SECP256K1_Q, SECP256K1_Q.bit_length()
    whole = search(q, 236)
    same = {"c": c_search(exe, q, 236) == whole,
            "sharded": sharded(search, q, 236, 7, ProcessPoolExecutor) == whole,
            "c_sharded": sharded(lambda *a: c_search(exe, *a), q, 236, 7, ThreadPoolExecutor) == whole}
    ok &= all(same.values())
    print(json.dumps({"q": "secp256k1", "u_floor_log2": 236, "nodes": whole[1], **same}), flush=True)
    T, cut = tick_bound(q), SECP256K1_CUT_LAST
    log_umin = 384 - cut
    t0 = time.time()
    shards = int(os.environ.get("WIDTH_SHARDS", "256"))
    if use_python:
        cands, nodes = sharded(search, q, log_umin, shards, ProcessPoolExecutor)
    else:
        cands, nodes = sharded(lambda *a: c_search(exe, *a), q, log_umin, shards, ThreadPoolExecutor)
    broken = {x: b for x in cands if (b := breaks(q, x))}
    ok &= not broken
    # The schedule this proves: one bit under the generic one at ticks 129 to cut, else generic.
    generic = [min(n, ((2 * q * q * (q - 1) >> 2 * t).bit_length() + 1) // 2) for t in range(T)]
    want = [min(n, 384 - t) if t <= cut else generic[t] for t in range(T)]
    schedule_ok = tick_widths(q, T) == want and generic[cut + 1] == want[cut]
    ok &= schedule_ok
    print(json.dumps({"q": "secp256k1", "search": "python" if use_python else "c", "cut_last": cut,
                      "u_floor_log2": log_umin, "shards": shards, "nodes": nodes,
                      "candidates": [hex(x) for x in sorted(cands)], "broken": len(broken),
                      "schedule_ok": schedule_ok, "seconds": round(time.time() - t0)}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
