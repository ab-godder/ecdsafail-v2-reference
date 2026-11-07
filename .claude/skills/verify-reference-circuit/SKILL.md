---
name: verify-reference-circuit
description: Drive the ecdsafail-v2 reference circuit generator from the shell and prove a change to qec/ keeps the circuit exact, inside the limits, and independent of the addend in its definitions. Use after any edit under qec/ or scripts/, before handing a change back, or when asked to verify, check, or measure the controlled point-addition circuit or its per-addend cost.
---

# Verify the reference circuit

The product is a Python generator, `qec/`, that emits a reversible circuit for |c>|R>|0> -> |c>|R + c*A>|0> on secp256k1. A user touches it through four commands. `scripts/emit.py` writes the circuit JSON. The checkers in `scripts/` simulate it. `scripts/accept.py` runs every invariant check as one gate. The Rust harness in `harness/` replays a circuit on the benchmark's own executor, which needs the unpublished benchmark checkout.

Everything runs in the shell. There is no server and nothing to keep alive. Read `INVARIANTS.md` at the repo root for why each check exists.

All commands below run from the repo root.

## Launch

Python 3.9 or newer with the standard library is the only dependency. There is nothing to install or build for the Python side.

```sh
python3 --version                                   # 3.9 or newer
.claude/skills/verify-reference-circuit/scripts/doctor.sh
```

The checkout is ready when `doctor.sh` exits 0. A fresh clone has no `verify-out/`, and `scripts/emit.py` does not create directories, so make the scratch directories once:

```sh
mkdir -p verify-out/emit
```

The recipes below and in `features/` also create what they write into, so each one runs on its own. The harness alone needs a build, and only when the benchmark checkout exists at `../ecdsafail-v2` (see `features/benchmark-harness.md`).

Teardown is the Cleanup section. Each command is a short-lived process that exits on its own.

## Doctor

Run this first, and again whenever a result looks wrong:

```sh
.claude/skills/verify-reference-circuit/scripts/doctor.sh
```

It is read-only. It reports the Python version, whether `qec` imports with the standard library only, the git head and the number of changed files under `qec/`, whether `verify-out/` is writable, whether another `accept.py` is running, and whether the benchmark checkout is present. It exits 1 on a missing Python, a broken import, or a missing `scripts/accept.py`.

Don't run two gates in one checkout at the same time. Each run deletes and rewrites `verify-out/accept/`. To run two at once, give the second its own `--out` directory.

## Drive

The gate is the main driver. Every check has an id, and `--only` runs a subset.

```sh
python3 scripts/accept.py --quick                     # about 80 s of CPU; every check except I2.secp-flat, 9 addends, 2 simulated
python3 scripts/accept.py                             # about 600 s of CPU; adds I2.secp-flat, 24 addends, all simulated
python3 scripts/accept.py --quick --only I5 I6 --sample 4 --seed 7
```

Wall time depends on cores and load, because the checkers run in parallel. On 10 cores at load average 6 to 15, quick mode took 17 s and full mode 103 s.

The ids are `I1`, `I2.model`, `I2.toy`, `I2.toy-hier`, `I2.components`, `I2.offset`, `I2.multiplier`, `I2.mutation`, `I2.secp`, `I2.secp-flat`, `I3`, `I5`, `I6`, `I7`, `I7.charge` and `R1`. Each prints a row with `PASS`, `FAIL` or `PENDING` and a one-line note. The last line starts with `gate PASS` or `gate FAIL`. The exit code is 0 only when no check fails. No check is pending now.

The feature files in `features/` give the exact commands for each user-facing feature. Start from `features/README.md`.

## Evidence

A proof is the command, its full output, its exit code, and the files it wrote. Capture it like this after any gate run:

```sh
mkdir -p verify-out; python3 scripts/accept.py --quick > verify-out/gate.log 2>&1; echo "exit $?" >> verify-out/gate.log
.claude/skills/verify-reference-circuit/scripts/evidence.sh <feature-id> verify-out/gate.log
```

`evidence.sh` copies `verify-out/accept.json`, the checker logs under `verify-out/accept/logs/`, the checker receipts, a sha256 and byte count of every emitted circuit, and `git rev-parse HEAD` plus `git status` into `verify-out/evidence/<feature-id>-<UTC time>/`. It prints that path. Evidence lives there and nowhere else. `verify-out/` is git-ignored.

Proof standards:

- Drive the real path a user runs, `scripts/emit.py` and `scripts/accept.py`. Don't call `qec` functions in a one-off script and report that as the proof.
- Capture the action and the resulting state. For emission, that is the meta line it printed and the file it wrote. For a check, that is the row and the evidence in `accept.json`.
- Check the side effects. `R1` in `accept.json` confirms the run changed no file git tracks, in content, type or mode. An emitted file has the byte count and sha256 that `evidence.sh` recorded.
- A pass of `I2.secp` alone proves little. Random 256-bit points cannot tell an exact circuit from an approximate one. Exactness rests on `I2.toy`, `I2.toy-hier`, `I2.components`, `I2.offset`, `I2.multiplier` and `I2.model`. `I2.secp` catches a wrong gate that runs only at 256 bits, on the addends and points it samples.
- `I7.charge` computes the live board's charge from the emitted file. It does not run the board.

## Cleanup

```sh
rm -rf verify-out/accept verify-out/emit verify-out/gate.log verify-out/accept.json
```

This removes the emitted circuits, the logs and summary of the last run, and the scratch emit directory. The 256-bit hierarchical circuit alone is 43 MB. It keeps `verify-out/evidence/`, which holds the proof. Kill only processes you started, by the pid your shell reported, never by name.

## Helpers

- `scripts/doctor.sh` is the read-only health check. Run `.claude/skills/verify-reference-circuit/scripts/doctor.sh`. Set `PY=/path/to/python` to check another interpreter.
- `scripts/evidence.sh <label> [file ...]` copies the last gate run and any extra files into `verify-out/evidence/<label>-<UTC time>/`. Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh limits` to save the last gate run alone.
