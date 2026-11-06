"""Classical integer model of the exact dialog in-place multiplier.

M(x, y) = (x, x*y mod q) for x != 0 and M(0, y) = (0, y). The circuit in
point_add.py mirrors these loops gate for gate, so this model is the oracle
for the iteration bound and the Bezout identity.
"""


def iteration_bound(q):
    """Proven cap on active Euclid iterations. bitlen(u) + bitlen(v) starts at
    2n or less and drops by at least 1 per step, because |v - u|/2 is at most
    (max(u, v) - 1)/2. It is at least 2 while v != 0, so after 2n - 2 steps
    (u, v) is (1, 1) or v is 0, and one more step clears v."""
    return 2 * q.bit_length() - 1


def euclid(q, x, iters):
    u, v, dialog = q, x, []
    for _ in range(iters):
        b0 = v & 1
        sw = b0 & (u > v)
        dialog.append((b0, sw))
        if sw:
            u, v = v, u
        if b0:
            v -= u
        v >>= 1
    return u, v, dialog


def bezout(q, y, dialog):
    """Inverse Euclid on (r, w) = y*(u, V) with V = 2v + u, from (u, V) = (1, 1)."""
    r, w = y, y
    for b0, sw in reversed(dialog):
        if sw:
            r, w = w, r
        w = (2 * w + (r if b0 else -r)) % q
    return r, w


def inplace_mul(q, x, y, iters=None):
    iters = iters or iteration_bound(q)
    u, v, dialog = euclid(q, x, iters)
    assert v == 0, "euclid did not terminate within the bound"
    r, w = bezout(q, y, dialog)
    assert x == 0 or r == 0
    return (r + w) * (q + 1) // 2 % q


def width_schedule(q, iters):
    """Bits needed for u and v before iteration i while v != 0: the smaller one
    takes at least 1 of the 2n - i bits iteration_bound counts."""
    n = q.bit_length()
    return [min(n, 2 * n - 1 - i) for i in range(iters)]
