# Invariants for changes to the circuit


## Current mux contract

The merged path now uses the paired-tick contract D1-D4 in [the mux proof](docs/mux-proof.md), superseding the older unrestricted canonical per-symbol Bézout-cell claims below. Logical T and R are canonical and scratch starts clean, but after a non-final B1 the physical word stores H = 2T + R in S4's encoded form, borrowing the next cell's g0 for its high bit. A semantic subtract follows B1; next.g0 must be 1 before holding; a final B1 requires R = 0. The held tick restores the borrowed symbol. Fallback primes retain the old step.

The additional sweep rotates each no-subtract run so B1 is last in Bézout order.
It runs in Euclid's tail, one cell behind the offset sweep, with which it commutes: the offset sweep's swap acts only on b of g0 = 1 cells, and each B1 step only on b of two g0 = 0 cells.
So Euclid hands Bézout the B1-swept tape, and inverse Euclid unsweeps it before undoing its steps; this supersedes the placement in S2's last sentence and S17's last sentence of the mux proof, not the maps they state.
On the raw tape the sweep moved again: Euclid hands Bézout the unswept tape, each tick applies the sweep step into its own cell through a carried toggle wire (`TwoBitMod._sweep_in`), and `_body` unsweeps after the Bézout pass; see DESIGN.md, "The B1 sweep inside the Bezout ticks". With a swap-bit wire per cell the tick reads s from the tape and the unsweep is `_b1_sweep_sw` reversed; see "Swap bits on the tape". The tape every tick reads, and so D1 to D3, are unchanged.
The tape at rest is block-coded, two wires freed per 15 ticks by `qec/tape_block.py`; on the swept tape every B1 follows a subtract cell in tape order, so the code reads each block's cells last first and relies on the same clause.
Cells read previous-symbol wires and can write the next g0 across a block edge, so Bézout decodes a block one tick before its first tick and codes it again one tick after its last, and two blocks are decoded across each edge for two ticks.
The gate's toy primes have at most one block, so the edge path is covered by `scripts/block_check.py` at 20 bits and above and by the secp256k1 rows, not by a toy instance.
The four-symbol alphabet, 384-tick bound and width schedule are unchanged. M and its inverse use clean states in their respective input/output layouts; top-level R remains on-curve or infinity with canonical encoding and zero scratch.

Euclid's first step starts on U = q, E = 0 and a clean cell.
Its specialised block borrows k = n - 1 - lsb(2^n - q + 1) unwritten tape wires, returning them zero before they become tape.
Both forward and inverse boundaries are clean; the block has no internal wires.
If the width, tape supply or block shape does not fit, the general step is used.
Every later step is unchanged.

M's first merged Bézout iteration starts with S = 0 and terminal cell A0, A1 or J, never B1 under the existing tick bound and sweeps.
It uses S XOR= g0*R with one CCX per bit, unchanged R/cell/next and no scratch; S keeps its wire list.
Only the terminal g0 is read, so its b lend is unchanged.
Every later tick and the fallback path retain their existing implementation; see S16 of the mux proof.

The corrected S13 distinguishes logical carries from physical borrowed first-probe wires: only on P = 1 do those wires equal the carries; on P = 0 their junk is unobserved except through P and reversed. Exported e and the second chain remain zero on the inert lane. Credit arXiv:2609.05625 section VI for the offset representation.

The exact-source independent review covers this pair contract; it does not grant upstream checker acceptance. The local held secp BDD probe stopped at the production 4,194,304-node budget. The earlier cell's cap-transfer proof does not automatically transfer. I7.charge, qualification, and the flat-upload operation cap remain open limits.

### Clean wires in the mux tick

Two changes in `qec/two_bit.py` take the merged tick from 26 clean wires to 19 at secp256k1, at +14 CCX per tick.
Each is exact on its own, and a change to the second chain or to `Prefix` must keep the statements below.

The second chain keeps only some of its low carries, in `_chain_up`, `_chain_down` and `_kept_carries`:

- **C1.** Unless bit j − 1 is perched, `chain[j]`, for 1 ≤ j ≤ a, is built by `_selected_carry` from `chain[j − 1]`, `L[j − 1]` and the selectors, and the same call with the same inputs clears it.
- **C2.** Going up, once kept carry i + 1 is built, the dropped carries between it and the previous kept one are cleared from i down.
  Each `chain[j − 1]` is still live when `chain[j]` is cleared.
  No perched bit rises between a dropped carry's build and its clear, because the carry into a perched bit is kept, so `L` and the selectors still hold what the build read.
- **C3.** Going down, falls run top down, and before fall(i) a missing `chain[i]` is rebuilt upward from the largest kept carry below i.
  No bit between that carry and i is perched or has fallen, and every perch above i has fallen and restored its selector.
  A perch below that carry changes only its own bit, the carry into it, and a selector that no higher bit reads.
  So every rebuilt carry reads what it read going up, and equals the one cleared then.
- **C4.** The carry out of a perched bit sits on its selector at no wire and is always kept.
  So is the carry into a perched bit, because the perch changes it before its fall.

`_kept_carries` picks the kept set that holds the fewest low carries on wires at once, while a segment is built or at the chain's top beside `chain[b]`, then the one that drops fewest, at 2 CCX per dropped carry.
With more than 16 free carries it keeps them all.
At secp256k1 the run is (12, 32), the perches are 11, 32, 33 and 34, and the kept set is {2, 6, 9, 11, 12}, with c12 on bit 11's selector.

The quiet run's prefix M releases its top below the last until its write, in `Prefix.release` and `Prefix.restore`:

- **R1.** Right after M is built, its top below the last, p = a·v, is cleared by `ccx(a, v, p)` and its wire freed, but only when none of a, p and v is a wire of the last block's chain.
  Otherwise one of them could hold that chain rather than its own value.
- **R2.** No gate between `release` and `restore` writes a, v or `M.top`, because by `Prefix`'s contract nothing else touches M's bits or live tops before `M.write`.
  The freed wire is not kept: the allocator hands it out again, and other gates write it before `restore`.
- **R3.** Until `M.write`, nothing outside `Prefix` reads M but `M.top`, in the quiet-run skip of both chains.
  The cleared top is read only by `write`, as the second-to-last block's p and as the last block's a.
- **R4.** `restore` recomputes a·v on a newly allocated zero wire just before `M.write`, and puts that wire in both places.
  `write` asserts that no top is still released.

At secp256k1 M has 20 bits in 3 blocks, and the cleared top is AND(L[12:24]).
P's matching top is not released, because `P.write` runs at the second chain's top, which is the tick's peak.

Small moduli barely reach these changes.
Below 2^14, on the merged path, the second chain drops a carry first at q = 97, and two adjacent carries first at q = 397, which is the first case where C2's order and C3's upward rebuild matter.
The release fires there only at q = 3,583, 12,799, 13,309 and 15,359.
A check of either change needs such moduli or secp256k1.


This page is for you if you change `qec/` to make the circuit cheaper. It says what `scripts/accept.py` protects and why, what no script can check for you, and how to hand a change back.

The benchmark that will judge your change is not published yet. `accept.py` encodes the invariants that the circuit must keep under any scoring rule we expect. A change that passes it is likely to be accepted once the benchmark ships. Nothing here is a promise about the final scoring.

## Run the gate

```sh
python3 scripts/accept.py --quick   # about 80 s of CPU, 9 addends built, the 256-bit upload simulated at 2 of them
python3 scripts/accept.py           # about 600 s of CPU, 24 addends built and simulated, plus the flat 256-bit circuit at G
```

The CPU times are user plus system time over the gate and its checker processes, measured with Python 3.14 on an Apple M4. The checkers run as parallel processes, so wall time depends on the cores and on the load from other work. With 10 cores at load average 6 to 15, quick mode took 17 s and full mode 103 s. A reviewer measured full mode at 298 s with load average 93 on the same 10 cores. On fewer free cores the wall time grows toward the CPU totals. Full mode peaks near 4.5 GB of memory across its processes. Python is stdlib only. The gate prints one row per check and writes the full evidence to `verify-out/accept.json`. It exits 0 only when no check fails. A check marked pending waits on the published benchmark and fails nothing.

Run `--quick` while you iterate and the full mode before you hand a change back. `--only I5 I6` runs a subset. `--sample N` and `--seed S` change the random addends that I1, I3, I5, I6 and I7 build and that `I2.secp` simulates. A row's seconds are the time the gate waited for it, so the first row that needs the addends carries their build time.

## What each check protects

The circuit implements |c>|R>|0> -> |c>|R + c*A>|0> for R on the curve or infinity and a fixed public point A on secp256k1, under the frozen profile `secp256k1-controlled-exact-v1`. Every check below guards one part of that contract or one limit of the benchmark.

**I1 interface.** The control is wire 0, x is wires 1..n, and y is wires n+1..2n, at n = 7 and at n = 256, for G and for every built addend. Every other wire is scratch, and the wire count stays at most 4,096. The benchmark reads the declared registers from these positions, so a reordered layout computes the right sum in the wrong place. The reference uses 2,311 wires.

**I2 exactness.** The profile requires an exact unitary, so the circuit must be right on every promised input R on the curve or infinity, including R in {A, -A, -2A}. I2 runs the existing checkers and parses what they print. It does not reimplement them.

- `I2.model` runs `scripts/model_check.py`.
At 13 small primes it traces the Euclid loop for every x, against the proven iteration bound and the width schedule.
It checks the classical multiplier model on every (x, y) below q = 300 and on 8 sampled y per x above, so the multiplier part is not exhaustive.
It also checks the remapped terminal state (1, 0), the unmapped-zero and shortened-pass controls, the standalone gate pass, and structurally validated hierarchical multiplier gates on x in [1, q) in both directions with clean scratch; every new boolean contributes to the exit status.
The offset model is checked against the independent reference recurrence and swept tape at every small-prime x; the gate pass also checks that exact tape and clean scratch.
A V = 0 caller is outside the pass's domain, because late tape takes U's vacated wires, so the row checks that such a caller does not reproduce the old generic result; the remap control checks the integer model, not the kernel flag guarantee.
The secp256k1 row checks the terminal and offset states on boundary and 2,000 sampled x, plus both multiplier directions on 832 nonzero-x width-edge lanes.
**Scratch pass domain and gate limit.**
The pass requires odd prime q, U = q, canonical 1 <= V < q at entry, and clean scratch.
M's standalone contract is narrowed to x in [1, q); x != 0 at M's entries rests on the kernel's flag code.
The build refuses curves with a finite point at x = 0 and addends with x(A) = 0 or x(2A) = 0.
The kernel's f0 AND f1 indicator remaps zero before M inverse and restores it after M, with gen zero on those states.
The direct evidence is the instrumented observation that the indicator equals [X = 0] before the remap and X != 0 at both M entries, on all 25,800 toy lanes and the secp256k1 special points.
The guarantee's control of record removes both remaps with every flag intact; the existing point gates then fail on 72 lanes, c = 0 with R = inf and c = 1 with R = -A, across all 36 toy instances.
The marker-removal, one-remap and old-gen-fold mutants are sensitivity controls; they do not test the nonzero-X guarantee.
The integer-model remap_control_ok row alone does not gate the kernel guarantee.
The dialog_mul component drops six zero-x lanes, giving 51 lanes and a component census of 50 / 9,615.
The gate_mul row excludes x = 0 in both directions, and the secp256k1 width-edge set has 832 nonzero-x lanes.
The pass_zero_outside_domain negative control confirms loss of the old generic zero-caller behavior, not correctness outside the domain.
Every figure carries both named limits: **I2.model's zero-caller row narrowed from the reference's gate**, and **M's standalone contract narrowed to x in [1, q); x != 0 at M's entries rests on the kernel's flag code**.
Any claim for this build must name both limits and the build refusal.

- `I2.toy` runs `scripts/check.py` on 36 toy instances of 4 to 10 bits, two addends on each of 18 curves. It simulates every point with c = 0 and c = 1 at gate level and requires all scratch back at zero. The gate lists the 36 (curve, addend) pairs itself and requires each exactly once, with 2 x (curve order) lanes. It also loads the flat toy file `scripts/emit.py` wrote, `toy127_circuit.json`, and simulates it the same way at the addend the gate asked for.
- `I2.toy-hier` runs the same check on the hierarchical emission. `check.py` parses the JSON text with `qec/flatten.py` and expands every call through its wire map and inverse flag, as the benchmark would. On each toy instance it also compares that expansion with the text `flatten.operations_json` writes. It then simulates the hierarchical toy file `scripts/emit.py` wrote, `toy127_circuit.hier.json`. The upload is hierarchical, so a broken call map would pass `I2.toy` and still upload a wrong circuit.
- `I2.components` runs `scripts/component_check.py`. It checks the 256-bit adders, doubler, negation and multiplier on boundary inputs near q and 2^256 that random points almost never reach.
- `I2.mutation` runs `scripts/mutation_control.py`. It confirms that the checker still rejects mutated circuits. At its fixed seed the reference catches 163 of 200 CCX-to-CX demotions, and the rest are equivalent circuits. The check fails below 163.
- `I2.secp` simulates the 256-bit hierarchical upload, expanded the same way as `I2.toy-hier`. Each addend runs on infinity, A, -A, -2A, 2A, 3A and random R, with c = 0 and c = 1. Quick mode runs G and one random addend with 8 random R. Full mode runs every built addend with 120 random R. It also requires that the JSON it simulated for G is the bytes `scripts/emit.py` writes.
- `I2.secp-flat` runs only in full mode. It checks the flat 256-bit circuit at A = G on the same points with 120 random R, and requires all 252 lanes.

Each checker must reach the coverage pinned in `COVERAGE` at the top of `scripts/accept.py`: 13 model primes, 36 toy instances and 25,800 lanes, each listed instance exactly once, 49 components and 8,951 lanes, 200 mutation trials, and 2 x (6 + random R) lanes per simulated addend. A receipt that covers less fails, even with zero failures.

**I3 gate set.** Every emitted primitive is X, CX, CCX or SWAP, with no measurement and no reset. The profile disables measurement cleanup, and a circuit of these four gates is a permutation with no phases. That is why basis-state exactness with clean scratch gives the exact unitary on superpositions. A single Z, CZ or measurement breaks that argument even when every basis check passes. I3 reads the emitted JSON at toy and 256-bit sizes and the JSON of every built addend, and counts the gate kinds in every toy instance and every expanded 256-bit upload.

**I5 addend independence of definitions.** The gate builds the 256-bit kernel hierarchically for [1]G to [8]G and the seeded random addends, 1 of them in quick mode and 16 in full mode. Every definition must be identical across them. The gate compares the sha256 of each definition after it parses the JSON and writes it back compactly, so whitespace in the file does not count. Key order does, because the digest keeps the order the file uses. Today all 263 definitions are identical, and A enters only through the top level, in the `cmodadd_const` calls in `_kernel` and the four exception compares and loads in `build` (`qec/point_add.py:110-190`). Keep it that way. A whole attack runs this kernel for many addends, and shared definitions let one upload of the expensive parts serve every addend.

**I6 cost at G and sample statistics.** I6 reports CCX, T = 7 x CCX, ops and wires at G. Over the addends the gate built, [1]G to [8]G and a seeded random sample, it reports the sample maximum and the sample minimum. These are statistics of that sample, not bounds. I6 shows the sample maximum T against the reference, 50,465,989, the largest T any sampled addend reached at commit `8ca9b0e`. I6 never fails on cost. It fails only when its counts disagree with the expanded upload that `I2.secp` simulated, at every simulated addend, or with the flat circuit at toy size and, in full mode, at G.

**I7 limits.** Wires stay at most 4,096. The hierarchical JSON for G, which the gate writes to `verify-out/accept/emit/`, stays under 104,857,600 bytes, and so does the JSON for every built addend. Every built addend also stays at most 4,096 wires. The reference is 43,422,940 bytes.

**I7.charge.** The next benchmark release charges a hierarchical upload with a formula that is not public. This check stays pending until it is.

**R1 run hygiene.** The gate records every tracked file and every untracked file git does not ignore, before and after it runs. For each path it records the file type, the permission bits and a hash of the contents, or of the link target for a symlink, without following links. Every checker writes under `verify-out/`, so a run that changes one of those files is a bug in the gate. It compares these records, not git status, so a second edit to a file you already changed still shows, and so do a lost executable bit and a file replaced by a symlink.

## What no script can check

**One generator for every curve.** The exhaustive toy check says something about 256 bits only because the same code emits both sizes. Only the constants and the bound-derived width schedule differ. A code path that runs only at n = 256, or only for secp256k1's q, escapes every exhaustive check. One mutation used to test this gate added 2.8 million identity gates only when n = 256. It passed all of I2 and was caught only by the upload limit. A wrong gate that runs only at n = 256 fails `I2.secp`, which simulates the 256-bit upload in both modes, but only on the points and addends it samples. Write size-generic code. If a trick needs a property of secp256k1's q, test that property on toy primes that share it, or prove it. The secp256k1 width cut is such a property, proved in DESIGN.md ("A width cut at secp256k1") with `scripts/width_check.py`; it changes only the widths at which the shared step code runs.

**The iteration bound is proven, not sampled.** The Euclid loop runs bitlen(q(q-1)) = 512 iterations at secp256k1 because u*v at least halves per step. Random inputs finish far earlier. The paper's approximate circuits use a shorter count and fail on rare inputs. Shrinking the bound needs a proof that covers every input, not a sample that passed.

**Random 256-bit sampling cannot tell exact from approximate.** Most CCX in the constant adders and the idle Euclid tails fire only on rare carries. A circuit that drops them passes random 256-bit lanes and is still wrong (README, exactness evidence item 5). Exactness rests on the exhaustive toy check, the 256-bit boundary components and the proven bound. A passing `I2.secp` alone proves little.

## Where the CCX are

At A = G the circuit has 7,209,419 CCX.

| Part | CCX | Share |
|---|---|---|
| Two dialog multipliers, 3,140,864 each | 6,281,728 | 87.1% |
| Lambda-squared subtract | 910,080 | 12.6% |
| Flags, exception loads, constant adds | 17,611 | 0.2% |

The multipliers and the lambda-squared subtract hold 99.8% of the CCX, and none of them reads A. Spend your effort there. That work also carries over to a future windowed kernel that takes the addend as a quantum input, because that kernel needs the same multiplier and subtract. Constant-adder tricks do not carry over, since a windowed kernel has no classical constant to exploit. They also act on the 0.2% that depends on A.

## Report cost over a sample of addends

Cost depends on A. [2]G and [3]G each cost 7,209,427 CCX, 8 more than G. Over 1,024 sampled addends at the reference, T ranged from 50,465,177 to 50,465,989. The benchmark may test any addend, and a whole attack runs hundreds of them.

The sample maximum that I6 prints is the largest cost among the addends the gate built. It does not bound the cost at other addends. The true maximum over every A is unproven. A change can raise the cost at an addend outside the sample, and the gate will not see it. So keep the cost of your change independent of A. A change inside the definitions does that, because I5 holds them identical across addends. If your change touches the top level, where A enters, argue its maximum over every A. A change that helps G and hurts another addend is a regression. Quote the I6 sample maximum with the sample it came from, not the cost at G.

## Don't trade ops for many small calls

The next benchmark release charges per call and per definition, and those terms may rise before it ships. Splitting the circuit into many small definitions or many calls can lower the op count and still raise the charge. Prefer changes that remove CCX and ops outright.

## Not settled yet

- The scoring rule and the caps of the next benchmark release.
- Which whole-attack route comes first. One candidate is the binary wrapper with 512 controlled additions of fixed points. The other is a windowed kernel with a coherent addend.

## Hand a change back

Open a pull request against this repository with:

1. The change, with nothing under `receipts/` edited by hand.
2. The output of `python3 scripts/accept.py` in full mode, pasted into the description, and its `verify-out/accept.json` attached.
3. The I6 sample maximum T, the sample it came from, and its delta against 50,465,989. For a change to the top level, add the argument for its maximum over every A.
4. For any change to the iteration bound, the width schedule or a modular reduction, the argument that it stays exact on every input.

## Last uncoded tape cell

On the merged pass, the final swept cell obeys g0 XOR b = 1 on the nonzero entry domain.
After the second reverse-order tick it is no longer read.
When it is uncoded, CX then X zeros b and releases its storage; after Bezout a fresh zero wire restores it with X then CX.
The pool adopts the restored tape mapping before inverse Euclid.
Freed physical storage may be reused and need not remain zero.
See `docs/tail-lend-proof.md` for the bound, fallback guard and evidence limits.

## Seven-tick whole-block capacity codec

The seven-cell block codes W[0..12], whose six inner windows (g_i, b_i, g_(i+1)) exclude 000, while W[13] remains raw.
The existing block reading order reads cells last first on the B1-swept tape.
There are 4,021 coded words and 8,042 inputs including both raw-b values.
The primitive list is pinned by canonical-JSON SHA256 8a43178163d443fa37f9f9e8da3431e876f476d8f1d4c87b0710760106813ed2.
Its 326 primitives include 36 CCX: the five-cell part costs 15, the re-tag costs 11, and the U1 hand pack costs 10.
Two clean ancillas at indices 14 and 15 return zero, and FREE = (10,) names the block wire returned zero.
W[13] is absent from every primitive operand.

The generic encode/decode plumbing and sweep order are unchanged.
Each decode supplies one fresh zero wire and each encode returns one wire to the allocator.
At secp256k1, 54 blocks code 378 cells and leave six uncoded, preserving the last-cell lend.
Euclid lends 56 wires: its two tail zeros plus one zero per coded block.
Two decoded blocks still overlap across each edge.
Item 4's borrowed tape storage must be checked on this layout at its four borrow boundaries; later physical reuse is distinct from a live FREE or codec-ancilla role.

This section supersedes the earlier fifteen-tick codec summaries and their 25-block, two-free-wire and 52-wire lending figures.
The whole-block domain requires all six inner windows, including the boundaries the previous three-part codec did not code.
Checks cover the emitted list exhaustively, both sweep states at every gate prime with coded blocks, M and its inverse with every wire compared, the merged block edges, and directed secp256k1 lanes.
The ordinary toy mutation row is finite-domain evidence; separate codec deletion controls test codec sensitivity.
Every figure retains the limit: **`I2.model`'s zero-caller row narrowed from the reference's gate**.
No checker admission or qualification is implied.

Kernel constant additions use a wrap flag and exclusive selected constants, clearing the flag against the canonical output.
Controlled reflection computes (a-b) mod q on all canonical b, including zero, and replaces the kernel's negation plus final addition.
Both gadgets emit inline gates, restore the control and return their flag/carries clean; cneg and its narrower standalone contract are unchanged.
The reflection row adds 670 lanes to the 49 existing component rows, whose names and populations are unchanged.
The flag trie shares x-bit chains where 0, a and x(2A) agree, with an uncontrolled zero leaf and c joining after the root.

### Kernel flag recode across M

At either M boundary, (c, f0, f1, gen) is one of 0000, 0110, 1000, 1100, 1010, 1110 and 1001.
With c = 0 the word is 0000, or 0110 with the zero marker; with c = 1 it is 1001 on the generic path and 1000, 1100, 1010 or 1110 on an exception.
A change to `_eq_flags`, the exception list or `build`'s flag compression must preserve these seven words, or the recode releases a nonzero wire.
CX(gen, c) and CX(gen, f0) move 1001 to 0101, then gen XOR= NOT c AND (f0 XOR f1) clears gen using one CCX and no temporary wire.
The equivalent XOR form relies on the absence of a c = 0 word with f0 = 0 and f1 = 1.
The other six words are unchanged, and the seven remaining (c, f0, f1) words are distinct.
The kernel releases zero gen before each M apply and allocates a zero wire before reversing the recode, restoring c, f0, f1 and gen before any subsequent reader.
The released physical wire may become M scratch or an output register; the reallocated wire carries the restored logical gen and is returned to build.
Four recodes cost four CCX in total and free one live wire throughout both M applies.
Direct observation checks both boundaries on all 25,800 toy lanes and twelve secp256k1 special-point lanes, including zero at release and reallocation, word restoration and final data/control/scratch.
Removing gen's clearing CCX is the freed-wire control of record; removing CX(gen, f0) is a sensitivity control, each checked over the full toy population.
The inherited narrowed zero-caller row and standalone M domain x in [1, q) remain unchanged.

## Five-cell tail code

The merged Bezout pass codes five raw tail cells when the tail lend applies, at least five cells are uncoded, and the width at tick T - 5 is at most six.
The code applies after tick T - 6 and its block encode check. On the default pass it reverses after block zero is encoded and before the lend restore; on M's held path, where block0_raw skips that encode, it reverses at the same point with block zero still raw.
Its 168 primitives cost 19 CCX per apply, return FREE = (2,) zero, and restore three clean ancillas.
The domain is all 230 L6 words on nine live wires; the exhaustive tail_code component row checks distinct images, zero FREE and ancillas, and exact reverse restoration.
A real-tape corpus reaches 84 of those words; the remaining words are covered by the exhaustive component row, not claimed as reached tape states.
The new row adds one component and 230 lanes, changing the census from 50 / 9,615 to 51 / 9,845 without changing any other row.
Both existing named limits remain: I2.model's narrowed zero-caller row, and standalone M restricted to x in [1, q) with nonzero entries resting on kernel flags and the non-full build refusal.
The tail code additionally requires T >= CELLS + 1 so its forward tick exists; the inverse asserts that the forward apply ran.

## Block-zero code

Harnesses decode each block with its own FREE count and sum those counts for tape widths and alias checks.

For canonical x < q, strips before the first subtract leave U = q > V, so that subtract swaps and the swept block-zero word lies in the 4,022-word D1 subset of the L6 domain.
On the merged path with a block, a 21-CCX pre-map followed by the existing codec frees W10 and W13, with seven pre-map ancillas returned before the codec allocates its two.
The extra free wire joins the clean pool that U borrows at the pass boundary, reducing external tape by one; with_free_slots restores both positions for later passes.
A one-CCX HOLD frees b0 while ticks 7 through 2 run and restores it before tick 1, without touching the g1 next-cell lend.
Nonmerged paths keep the original codec, and no new per-prime domain refusal is needed.
The exhaustive D1 component row checks encode and HOLD both ways, clean FREE and ancillas, and distinct images, increasing the census only to 52 rows and 13,867 lanes.

## Pooled tick schedule

Selected merged ticks use a cap-four reversible pebble search on the second chain's unperched prefix, derived from that shape's length and quiet-run top wire.
The construction applies when the sole low perch is at a-1 and a cap-four schedule exists; otherwise that shape keeps the original kept-carry schedule.
At secp256k1 the 106 block-edge ticks, tail-code edge 378, and ticks 1 and 0 use the schedule, adding six CCX each.
At the other measured merged primes below 20,000 where this schedule runs, CCX changes by zero to minus two per use, with no wire change or scratch increase.
The uncoded terminal cell is never B1, so when the existing tail lend applies, tick T - 2 uses the no-prev shape and the terminal b release moves before that tick.
That released b is unread until the existing end-of-pass X and CX restore it.
These changes leave the nonmerged path and no-tail-lend instances unchanged.

The codec's primitives 4, 128, 145, 162, 188 and 197 use the six checked CX substitutions, giving a 30-CCX reversible code.
Its contract holds on all 8,042 window words and all 4,022 pre-map images; 6,254 window words have changed codewords, and each block uses the same shape for encoding and decoding.
Each of the six wrong-control substitutions fails the domain contract.

The merged Euclid sweep starts at mark(T - 2) and carries one mark equal to NOT g_k AND NOT g_(k-1), updating it before the next block encode and clearing it before block zero is encoded.
Each update toggles the mark by NOT g_(k-1) AND (g_k XOR g_(k-2)), so every g read is still raw.
The original carried-mark construction saved 382 CCX per Euclid pass at 384 ticks.
The endpoint identities and skipped terminal transition now remove four further CCX per pass, sixteen per point addition; nonmerged sweeps use the endpoint identities without the merged terminal skip.
