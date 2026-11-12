"""The Bezout tick's wire moves, at small primes whose high word H has two or more bits.

    python3.12 -I -B scripts/wires_check.py

The dialog multiplier runs forward and inverse on every x, against every y below 300 and 6
sampled y above, and every output and scratch wire must be exact (runs_check's lane count).
The primes cover e's wait on a borrowed H wire (len(H) >= 2) and the handed-over bit on l0.

Controls, which must fail at some prime: e left on its wire when it moves to H, the B1 bit left
in k0 at the tick's end, and c1's update without its CCX.
"""
import sys
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from runs_check import bad_lanes
from qec.builder import Circuit
from qec.two_bit import TwoBitMod

PRIMES = [31, 61, 127, 241, 251, 509, 1021]
SOURCE = ROOT / "qec" / "two_bit.py"


@contextmanager
def dropping(kind, snippet, after=None):
    """Drop every gate of this kind emitted from the line of two_bit.py that reads snippet, the
    first such line below the one that reads after when after is given."""
    lines = [line.strip() for line in SOURCE.read_text().splitlines()]
    start = lines.index(after) + 1 if after else 0
    target = lines.index(snippet, start) + 1
    real = getattr(Circuit, kind)

    def gate(self, *wires):
        caller = sys._getframe(1)
        if caller.f_code.co_filename == str(SOURCE) and caller.f_lineno == target:
            return
        real(self, *wires)

    setattr(Circuit, kind, gate)
    try:
        yield
    finally:
        setattr(Circuit, kind, real)


def sweep(label):
    total = 0
    for q in PRIMES:
        mod = TwoBitMod(q)
        assert mod.merged, q
        for inverse in (False, True):
            bad = bad_lanes(q, inverse)
            total += bad
            if label == "build":
                print(f"q={q} H={mod.n - mod.window} inverse={inverse}: {bad} bad lanes")
    return total


def main():
    failed = sweep("build") > 0
    controls = [
        ("e kept on its wire", ("ccx", "C.ccx(P.top, u, e)")),
        ("B1 bit left in k0", ("cx", "C.cx(l0, k0)", "# l0 holds only the B1's handed-over bit again; take that back out of k0.")),
        ("c1 update without its CCX", ("ccx", "C.ccx(l0, nb1, c1)")),
    ]
    caught = 0
    for name, args in controls:
        with dropping(*args):
            bad = sweep(name)
        print(f"control {name}: {'caught' if bad else 'MISSED'}")
        caught += bad > 0
    ok = not failed and caught == len(controls)
    print(f"wires_check: {'PASS' if ok else 'FAIL'}; controls caught {caught} of {len(controls)}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
