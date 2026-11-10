"""Sweep two_bit flags on the tick bench: python3.12 -I -B scripts/tick_sweep.py 'NAME=val' ..."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qec import two_bit
import tick_bench

base = None
for spec in sys.argv[1:] or ["BASE=0"]:
    for part in spec.split(";"):
        name, val = part.split("=", 1)
        if name != "BASE":
            setattr(two_bit, name, eval(val))
    C, block, bounds = tick_bench.build(3)
    C, bounds, _ = tick_bench.post_pass(C, block, bounds)
    end, _ = tick_bench.schedule(C.prims, C.num_wires)
    fins = [max(end[bounds[j]:bounds[j + 1]]) for j in range(3)]
    ccx = sum(1 for k in C.prims[4 * bounds[0]:4 * bounds[1]:4] if k == 2)
    steady = (fins[-1] - fins[0]) / 2
    print(f"{spec:50s} wires {C.num_wires} ccx {ccx} first {fins[0]} steady {steady}", flush=True)
