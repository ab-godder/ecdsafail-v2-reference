# Exact controlled point addition for ecdsafail-v2

## Start here

This repository holds the reference circuit for the ecdsafail-v2 controlled point addition, an exact and complete |c>|R>|0> -> |c>|R + c*A>|0> on secp256k1. Fork it to make the circuit cheaper. A live board tests uploads, and its benchmark source is not public yet, so `scripts/accept.py` checks the invariants a change must keep.

```sh
python3 scripts/accept.py --quick   # about 80 s of CPU, stdlib Python only
python3 scripts/accept.py           # about 600 s of CPU, run it before you hand a change back
```

The checkers run in parallel, so wall time depends on your cores and load. On an Apple M4 with 10 cores at load average 6 to 15, quick mode took 17 s and full mode 103 s.

- `INVARIANTS.md` says what each check protects, what no script can check, where the CCX are, and how to hand a change back.
- `.claude/skills/verify-reference-circuit/` is a Claude Code skill that drives the generator, the checks, the cost report and the limits, with a feature map under `features/`.

## This branch against the reference

Costs at A = G, from `qcb diagnostic check` on the hierarchical upload. This branch has not been entered on the board, so these are local results on sampled inputs.

| | Reference | This branch |
|---|---|---|
| `lowered_t` | 50,465,933 | 28,810,880 |
| `qubit_time` | 164,135,627,418 | 74,522,959,296 |
| `toffoli` | 7,209,419 | 4,115,840 |
| `allocated_wires` | 2,311 | 1,632 |
| `scheduled_ticks` | 71,023,638 | 45,663,578 |
| `expanded_ops` charge | 20,960,620 | 15,647,747 |
| upload bytes, hierarchical circuit JSON | 43,422,940 | 33,674,606 |

## Result on the live board

The board tested this circuit on 2026-11-06 as `reference-c4`, evaluation `eval-6d612fd4c6156cdb`, in round 2. It ended `TESTED_UNQUALIFIED`, with 0 of 10 boundary cases and 0 of 10,000 seeded draws failing. That is a sampled pass and not a qualification. Its raw costs are the line to beat.

| | Reference | Board limit |
|---|---|---|
| `lowered_t` | 50,465,933 | |
| `qubit_time` | 164,135,627,418 | |
| `toffoli` | 7,209,419 | |
| `allocated_wires` | 2,311 | 4,096 |
| `scheduled_ticks` | 71,023,638 | |
| `expanded_ops` charge | 20,960,620 | 21,696,718 |
| encoded records | 974,699 | 10,000,000 |
| listed call wires | 1,055,748 | 10,000,000 |
| upload bytes | 43,423,435 | 104,857,600 |

The charge leaves 736,098 of room, about 3.4% of the limit, and a call costs about 64 of it.

To start from the circuit on the board, emit the hierarchical file. The flat file is 822,274,918 bytes, which is over the upload limit.

```sh
python3 -I -B scripts/emit.py --hier --out my-circuit.json
```

The SHA-256 of `my-circuit.json` is `01eab537d5cd4c81680c17d94f0e664f5b616fc9e4818dc15937c839cdbfee00`. The upload wraps those exact bytes in a `SubmittedMethod` with `public_key` set to G. The top-level `artifact_hash` in the board's receipt is the SHA-256 of that upload, `666d5595e59e0755f91acf327cac07884e0fc58226d6abe8ee9b0845c7018457`. Another wrapper, such as one with a different `method_id`, changes that hash and leaves the circuit unchanged, so compare `my-circuit.json` to check your starting point.

## What this circuit is

This directory holds a generator for the arithmetic kernel that the `secp256k1-controlled-exact-v1` profile asks for: |c>|R>|0> -> |c>|R + c*A>|0> for a fixed public point A. It covers every R on secp256k1, including infinity and R in {A, -A, -2A}. It is built from the exact core of Schrottenloher's dialog Euclid and Bezout multiplier (arXiv 2606.02235). It is emitted in the challenge's own `ArithmeticCircuit` JSON and checked with the challenge's own executor.

**It does not qualify at 256 bits.** It is exact and complete, and its hierarchical upload passes the live board's sampled test. No route to qualification exists in the repository today. See "Unresolved gates" below.

`DESIGN.md` records the contract reading, the three designs compared, the data shape, and the completeness argument.

## Run

```sh
scripts/run_all.sh                                # about 2 minutes; rewrites everything in receipts/ and artifacts/
ECDSAFAIL_V2=../ecdsafail-v2 scripts/run_all.sh   # also builds and runs the harness against the benchmark checkout
```

Without `ECDSAFAIL_V2`, `run_all.sh` skips the harness steps and says so, and the harness receipts keep their committed contents.

Individual steps:

```sh
PY=python3
$PY -I scripts/model_check.py                                  # classical dialog model: Euclid trace exhaustive at small q, multiplier sampled
$PY -I scripts/check.py --toy-bits 4 6 7 8 9 10 --toy-per-size 6   # gate-level exhaustive check on toy curves
$PY -I scripts/mutation_control.py                             # negative control for the checker
$PY -I scripts/component_check.py                              # 256-bit components on boundary inputs
$PY -I scripts/check.py --secp --secp-random 120               # 256-bit full circuit, sampled
$PY -I scripts/emit.py --out artifacts/secp256k1_G_circuit.json.gz
gzip -dc artifacts/secp256k1_G_circuit.json.gz | harness/target/release/qec-harness --circuit - --cases artifacts/secp256k1_cases.json
```

Python is stdlib only. The harness is a Rust crate whose dependencies are path references into a sibling checkout of the benchmark at `../ecdsafail-v2`, next to this repository. That checkout is not public until the benchmark is published. It builds offline into `harness/target` and writes nothing into the sibling repo. Receipts were produced against sibling `HEAD 9152f4d` with its existing uncommitted diff (sha256 prefix of `git diff HEAD` is `3eaa988f0901d0d1`).

## Layout

- `qec/model.py` is the classical integer model of the multiplier, the proven iteration bound, and the Euclid width schedule.
- `qec/builder.py` holds the circuit, allocator, and recordable blocks that replay forward or reversed under wire remapping.
- `qec/arith.py` holds exact CDKM add, subtract and compare, each optionally controlled, plus modular add, double and negate.
- `qec/point_add.py` holds the dialog multiplier, the gen-controlled kernel, and the four-case completeness wrapper.
- `qec/sim.py` is an independent bit-sliced simulator. `qec/curve.py` is a reference curve law. `qec/emit.py` streams challenge-format JSON.
- `scripts/accept.py` is the acceptance gate. `INVARIANTS.md` explains each of its checks.
- `harness/` runs an emitted circuit through `v2_ir` parsing and validation, `v2_runtime::BitSlicedExecutor`, a `SubmittedMethod` round trip, `QualificationGraph::evaluate`, and `v2_lowering`.
- `artifacts/` holds the emitted circuits, test lanes, and the toy `SubmittedMethod`. `receipts/` holds every measured result quoted below.

## Results

All numbers are measured unless marked as estimates.

| | toy q=127 | secp256k1, A = G | cap |
|---|---|---|---|
| wires | 57 | 1,632 | 4,096 |
| ops | 12,277 | 14,020,218 | 10,000,000 |
| CCX | 3,535 | 4,115,840 | |
| T at 7 per CCX | 24,745 | 28,810,880 | |
| measurements, R or Hmr | 0 | 0 | 0 |
| flat circuit JSON bytes | 473,586 | 591,016,012 | 104,857,600 |
| simulated lanes, `qec/sim.py` | 254 of 254 correct (exhaustive) | 510 of 510 correct (sampled) | |
| dirty scratch or phase lanes | 0 | 0 | |

The harness has not been re-run on this circuit, because it needs the unpublished benchmark checkout. `receipts/harness_*.json` and `artifacts/toy127_submission.json` still describe the reference. They recorded `is_qualified` false.

Where the 256-bit CCX go (`receipts/secp256k1_breakdown.json`). One dialog multiplier is 1,964,330 CCX: two Euclid passes of 391,681 each, one Bezout pass of 1,171,212, and 9,756 to start and finish Bezout and pack the dialog. The kernel uses it twice. The lambda-squared subtraction is 169,569 CCX. Flags and constant adds make up the remaining 17,611.

The paragraphs below record how the reference reached 7,209,419 CCX. Later changes took it to the figures above: offset residues and W form in Bezout, an integer square with folding and one Karatsuba level, a one-chain Euclid step, a tighter iteration bound, reusing the Euclid registers' wires in Bezout, a packed dialog, and control order chosen for depth.

The controlled subtract and the controlled modular add run a controlled CDKM ripple instead of copying `ctrl*a` into an AND register. A controlled subtract costs 3w CCX instead of 4w, and a controlled modular add went from 10n + 2 to 9n + 4. That took the 256-bit point addition from 8,686,084 to 7,967,236 CCX, 8.3% fewer.

Classical constants are no longer copied into a register for a full adder. The carry chain reads the constant's bits, so it skips carries that are always 0 and aliases a carry that equals an input bit. Each modular reduction runs one carry chain of s + 2^(n+1) - q. Its carry out is [s >= q]. It lands in the spare top bit, because the top bit and the carry into it are never both 1 when s < 2q. The same carries then drive the conditional subtraction. The reduction costs 3n - 3 CCX instead of 4n + 2, so a controlled modular add now costs 8n - 1, which is 2,047 at n = 256. A modular double costs 3n - 15 instead of 4n + 4 at secp256k1, because 2s has a zero low bit and the next three bits of 2^(n+1) - q are 0, so the first carries are always 0. Its flag is the low bit of the result, so the flag wire becomes that bit with no gate. Negation adds the sparse 2^n - q - 1 and then complements. That took the point addition from 7,967,236 to 7,209,419 CCX, 9.5% fewer, and from 26,931,713 to 18,995,461 ops. The carry chain still ripples the full width, because exactness needs it.

Exactness evidence, strongest first:

1. **Exhaustive gate-level check.** `receipts/toy_exhaustive.json` covers 36 instances: 18 primes of 4, 6, 7, 8, 9 and 10 bits, 2 addends each, every curve point plus infinity, and both control values. It finds 0 failures, with all scratch returned to zero. The same generator code emits the 256-bit circuit. Only the constants and the bound-derived width schedule differ.
2. **Challenge executor.** `receipts/harness_toy127.json` and `receipts/harness_secp256k1.json` replay the emitted JSON on `v2_runtime::BitSlicedExecutor`. Expected values come from `qec/curve.py`, which shares no code with the generator.
3. **256-bit boundary components.** `receipts/secp256k1_components.json` covers modadd with and without ctrl on edge-by-edge pairs around q and 2^255 and at q - f and 2^256 - 1 mod q, where f = 2^256 - q, plus moddbl on the same values and cneg. It runs cmodadd_const and gt_const_into with ctrl 0 and 1, gt_const_into without ctrl, and add_const under ctrl with a carry out, each on the kernel's constants for G (q - a, q - b, 3a and a) and on sparse and dense edge constants (1, f, 2^255, q - 1, q - 2, 2^255 - 1 and 2^256 - 1 mod q), with inputs on both sides of each constant. It also runs the multiplier on x = 0 and on x values whose Euclid run takes 511 of the proven 512 iterations. All pass.
4. **Proven iteration bound.** bitlen(u) + bitlen(v) drops by at least 1 per active Euclid step, so 2n - 1 iterations always suffice. That is 511 at secp256k1. `receipts/model_check.jsonl` shows the observed worst case never exceeds the bound at any of the 13 tested q. It equals the bound at all 13, so the bound is tight on the tested q.
5. **Checker sensitivity.** A must-fail control (an extra X on the output) fails every lane. Exhaustively, 177 of 200 random CCX-to-CX demotions are caught at q = 127, and the survivors are correct on the whole promise domain. Some can fail off the curve, so they are not equivalent circuits. `scripts/multiplier_check.py` runs its own control on the multiplier at q = 59, and catches 194 of 200. Most CCX in constant adders and idle Euclid tails fire only on rare carries, so random 256-bit sampling cannot separate an exact circuit from an approximate one, which is why items 1, 3 and 4 carry the exactness claim.

The circuit is pure X, CX, CCX and SWAP, so it is a phase-free permutation. Basis exactness with clean scratch on the promise domain therefore gives the exact unitary on superpositions and entangled references (QNT-04). Outside the promise it is still a permutation, which satisfies QNT-06.

## Unresolved gates

These gates are listed in the order they blocked the flat upload on the earlier benchmark. Gates 1 and 2 do not block the hierarchical upload on the live board, which has the limits in "Result on the live board" above. Gates 3 to 5 remain open.

1. **Operation cap.** 14,020,218 ops against a cap of 10,000,000 (`MAX_DECLARED_RECORDS`). One multiplier alone is 6,786,012 ops. The paper reports about 2^24 total gates per addition even for its approximate circuits, so no published point addition fits this cap.
2. **Upload cap.** The service parses plain JSON, at about 43 bytes per op. 100 MiB therefore holds about 2.4M ops, roughly 6x fewer than this circuit needs. This also binds before the operation cap does.
3. **Whole-attack T budget (estimate).** DEF-02's binary wrapper makes 512 controlled additions, which comes to about 2.6e10 T against the 3.25e9 cap. A w = 16 windowed interface (QNT-02) would need about 32 additions plus exact lookups, roughly 2.0e9 T. That would fit, but it means a coherent-addend kernel, which is not built here.
4. **Qualification route.** `QualificationGraph` requires a universal secp256k1 certificate checked by a `SoundTheoremChecker`, a Shor wrapper certificate, a schedule, observations, and a retry model. The production registry holds no sound checker, so every artifact ends `is_qualified: false` with `No sound theorem checker exists`. This candidate supplies none of those artifacts.
5. **Formal proof.** Exactness rests on exhaustive toy verification of a size-parametric generator, full-size boundary tests, and the proven iteration bound. That is strong evidence, not a machine-checked proof.

## Optimization headroom (estimates, not built)

These estimates are against the earlier benchmark's caps, not the live board's limits.

- The dialog garbage compression from the paper's Fig. 1 saves qubits, not ops.
- The exact equivalent of the paper's pseudo-Mersenne reduction still needs a full-width carry, so that saving is smaller here than in the paper.
- A further 2x reduction would bring the ops under the operation cap, but the flat file would stay about 3.9x over the upload cap.
