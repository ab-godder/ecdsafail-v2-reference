# Limits check

The limits check confirms the circuit fits the benchmark's interface and caps. It covers the register layout, the 4,096-wire cap, the 104,857,600-byte upload cap on the hierarchical JSON, and the gate set of X, CX, CCX and SWAP with no measurement. `I7.charge` computes the live board's charge on the hierarchical upload at G. `scripts/accept.py` runs these as `I1`, `I3`, `I7` and `I7.charge`.

## Sub-features

- `limits-layout` checks control = wire 0, x = wires 1..n and y = wires n+1..2n, at n = 7 and n = 256, for G and every built addend.
- `limits-wires` checks at most 4,096 wires for G and every built addend.
- `limits-upload` checks the hierarchical JSON for G and for every built addend is under 104,857,600 bytes.
- `limits-gates` checks every emitted primitive is X, CX, CCX or SWAP, in the emitted JSON, the JSON of every built addend, and every expanded 256-bit upload.
- `limits-charge` reports the live board's charge at G against its limit.

## How to get to it (user POV)

- Run `python3 scripts/accept.py --only I1 I3 I7 I7.charge`.
- Run `python3 scripts/accept.py --quick`, which includes all four.

## Driving it with the shell

Preconditions:

- `doctor.sh` exits 0.
- No other `accept.py` runs in this checkout.

- **Run the checks.** Run `mkdir -p verify-out; python3 scripts/accept.py --quick --only I1 I3 I7 I7.charge > verify-out/gate.log 2>&1; echo "exit $?" >> verify-out/gate.log`. It spends about 60 s of CPU, which took 16 s on 10 cores at load average 15. Its summary line starts with `gate PASS: 3 pass, 0 fail, 1 pending` and its last line is `exit 0`.
- **Read the layout.** `I1`'s note reads `control=[0], x=1..n, y=n+1..2n at n=7, n=256 and all 9 built addends; 2311 wires, at most 2311 over addends (reference 2311, cap 4096)`.
- **Read the upload size.** `I7`'s evidence has `hier_json_bytes_G` 43,422,940 and `largest_hier_json_over_addends` under the cap. Run `wc -c verify-out/accept/emit/secp256k1_G_circuit.hier.json` and confirm the same byte count on disk.
- **Read the gate set.** `I3`'s evidence lists `["CCX", "CX", "SWAP", "X"]` for `toy127`, `toy127.hier`, `secp256k1_G.hier` and each built addend from `[1]G` on, and `toy_instances_counted` 36.
- **Read the charge.** `I7.charge` has status `pass`, and its evidence holds the charge, records, listed call wires and call depth.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh limits verify-out/gate.log`. The evidence directory holds `accept.json` and `emitted.sha256` with the byte count of each emitted circuit.

## Gotchas

- The caps here are the ones this repo expects from the next benchmark release. The published benchmark also caps a flat upload at 10,000,000 operations, which the 256-bit circuit exceeds. README "Unresolved gates" describes that.
- `I7.charge` applies the board's published formula to the file. A pass does not mean the board accepted anything.
- The charge per call and per definition may rise before the benchmark ships. A change that lowers ops by adding many small definitions or calls can pass `I7` and still cost more.
- A code path that runs only at n = 256 can break `I7` while every toy check passes. A run of `I7` at 256 bits is the only check that sees it.
