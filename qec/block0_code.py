"""Swept block-zero D1 code: first subtract swaps for canonical x < q.

P and HOLD from the Round11 checked handover, SHA256
2c6bae09f85cf2f2eb5bafc3cdeba0f4ca8235d0d2425ee1197347cb300c3ed3.
"""
from . import tape_block

FREE = (10, 13)
P = (
    ('cx', 12, 13),
    ('x', 12),
    ('x', 10),
    ('ccx', 12, 10, 15),
    ('x', 10),
    ('x', 8),
    ('ccx', 15, 8, 16),
    ('x', 8),
    ('x', 6),
    ('ccx', 16, 6, 17),
    ('x', 6),
    ('x', 4),
    ('ccx', 17, 4, 18),
    ('x', 4),
    ('x', 2),
    ('ccx', 18, 2, 19),
    ('x', 2),
    ('x', 0),
    ('ccx', 19, 0, 20),
    ('x', 0),
    ('cx', 15, 12),
    ('cx', 11, 13),
    ('ccx', 12, 13, 11),
    ('cx', 11, 13),
    ('cx', 15, 12),
    ('cx', 16, 15),
    ('cx', 9, 13),
    ('ccx', 15, 13, 9),
    ('cx', 9, 13),
    ('cx', 16, 15),
    ('cx', 17, 16),
    ('cx', 7, 13),
    ('ccx', 16, 13, 7),
    ('cx', 7, 13),
    ('cx', 17, 16),
    ('cx', 18, 17),
    ('cx', 5, 13),
    ('ccx', 17, 13, 5),
    ('cx', 5, 13),
    ('cx', 18, 17),
    ('cx', 19, 18),
    ('cx', 3, 13),
    ('ccx', 18, 13, 3),
    ('cx', 3, 13),
    ('cx', 19, 18),
    ('cx', 20, 19),
    ('cx', 1, 13),
    ('ccx', 19, 13, 1),
    ('cx', 1, 13),
    ('cx', 20, 19),
    ('ccx', 20, 13, 21),
    ('cx', 21, 7),
    ('cx', 21, 13),
    ('x', 7),
    ('ccx', 20, 7, 21),
    ('x', 7),
    ('x', 20),
    ('ccx', 12, 20, 13),
    ('x', 20),
    ('x', 0),
    ('ccx', 19, 0, 20),
    ('x', 0),
    ('x', 2),
    ('ccx', 18, 2, 19),
    ('x', 2),
    ('x', 4),
    ('ccx', 17, 4, 18),
    ('x', 4),
    ('x', 6),
    ('ccx', 16, 6, 17),
    ('x', 6),
    ('x', 8),
    ('ccx', 15, 8, 16),
    ('x', 8),
    ('x', 10),
    ('ccx', 12, 10, 15),
    ('x', 10),
    ('x', 12),
)

HOLD = (
    ('x', 12),
    ('cx', 11, 13),
    ('ccx', 12, 13, 11),
    ('cx', 11, 13),
    ('x', 12),
    ('x', 13),
)


def _apply(C, W, prims):
    for p in prims:
        if p[0] == "x":
            C.x(W[p[1]])
        elif p[0] == "cx":
            C.cx(W[p[1]], W[p[2]])
        else:
            C.ccx(W[p[1]], W[p[2]], W[p[3]])


def encode(C, W):
    ancillas = C.alloc(7)
    _apply(C, list(W) + [None] + ancillas, P)
    C.free(ancillas)
    tape_block.encode(C, W)


def hold(C, W, inverse=False):
    _apply(C, W, reversed(HOLD) if inverse else HOLD)
