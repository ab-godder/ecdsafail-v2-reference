# Exact controlled point addition for `secp256k1-controlled-exact-v1`

## What the challenge actually requires

Read from a checkout of the ecdsafail-v2 benchmark repository (read-only), which is not public until the benchmark is published.

- Profile `profiles/secp256k1-controlled-exact-v1.json` sets `arithmetic_contract: EXACT_UNITARY`, `approximate_arithmetic: DISABLED`, `measurement_cleanup: DISABLED`, target `CLIFFORD_PLUS_T`.
- SPEC (`ecdsa specs/synthesized-v2/SPEC.md`) DEF-01 requires complete controlled addition of fixed public points. QNT-03 requires identity on control 0 and exact addition on control 1. QNT-05 requires the default exact profile to handle every reachable case, including infinity and R in {A, -A, -2A}. QNT-04 requires correctness on superpositions and with entangled references. SEM-03 rejects `R` and `Hmr`.
- Encoding (QNT-01) is little-endian affine coordinates with `(0,0)` as infinity.
- A submission is a `SubmittedMethod` JSON (`assurance/qualification/src/bundle.rs`) whose `circuit` is an `ArithmeticCircuit` of `X`, `CX`, `CCX`, `SWAP`, `Z`, `CZ`, `CCZ`, `NEG`, `CONTROLLED_NEG` (`contracts/ir/src/arithmetic.rs`).
- Service caps (`service/src/controller.rs`, `capacity/`, `contracts/api/src/quota.rs`) are 4096 wires, 10,000,000 operations, and a 100 MiB upload. The upload is parsed as plain JSON with no compression.
- `limits.json` caps compiled T at 3.25e9 for the whole attack. `compiler/lowering` lowers each CCX to the 7-T Toffoli.
- Qualification (`assurance/qualification/src/graph.rs`) needs a universal certificate checked by a `SoundTheoremChecker`. The production registry holds none, and `qcb diagnostic check` executes nothing. No artifact can be `is_qualified` in this repository today.

So Schrottenloher's published circuits (arXiv 2606.02235) cannot pass as published. They use MSB-only comparisons, pseudo-Mersenne shortcuts that skip rare carries, a fixed iteration count tuned for random inputs (failure up to 2^-13.3), measurement-based uncomputation, and they ignore doubling and infinity. Each of these violates one of the clauses above. What carries over is the exact structure underneath: the dialog Euclid algorithm and the out-of-band Bezout reconstruction that turns an inverse into an in-place multiplication.

## Designs compared

Costs are Toffoli estimates at n = 256 before implementation. Every design here is exact and uses only X, CX, CCX, SWAP.

**A. Exact dialog multiplier (chosen).** Point addition follows the paper's Algorithm 1. It uses an in-place multiplier M(x, y) = (x, xy) and its reverse, plus one out-of-place lambda-squared subtraction. M runs Euclid forward to turn x into a dialog bit-vector, then runs Bezout on (y, 0), then runs Euclid backward. Making it exact means three changes. The iteration count becomes the proven bound bitlen(q(q-1)) = 2n. Comparisons and modular reductions become full width. Uncomputation uses reverse circuits instead of measurement. Euclid widths still shrink, because u*v at least halves each iteration, so the width can be min(n, bitlen(q(q-1) >> i)). Estimated cost is about 58n^2 per M and about 136n^2, roughly 9M Toffolis, per addition. Qubits are about 9n.

**B. Kaliski or binary-EEA inversion plus schoolbook multiplications (Roetteler et al. 2017, Haner et al. 2020).** This design computes lambda with an explicit inverse, then uncomputes the inverse. That costs two inversion passes plus at least three multiply-and-uncompute pairs. Published exact counts are about 2.5e8 Toffolis per addition for Roetteler. Haner's improvements remain several times above design A. Rejected because it costs more on every axis that binds.

**C. Projective complete formulas (Renes, Costello, Batina).** These need no inversion and are complete. They are rejected because projective coordinates are not canonical. An in-place R -> R + A on the declared affine encoding still needs a normalizing inversion to uncompute the old point. That adds design A's cost on top of the formulas.

## Completeness on the fixed addend

For a fixed public A = (a, b), the generic kernel fails only on R in {inf, A, -A, -2A}. On those inputs x - a = 0 or a - x' = 0, and M is not a bijection there. All four inputs and their sums {A, 2A, inf, -A} are public constants. The circuit therefore:

1. computes flags f_k = c AND [R == E_k] for the four exceptional inputs, and gen = c AND NOT(f1 xor f2 xor f3 xor f4);
2. XORs E_k xor O_k into the accumulator under f_k;
3. runs the generic kernel controlled only through gen, using M^-1 ... M. With gen = 0, every gen-controlled step is the identity, and M^-1 followed by M cancels exactly at the circuit level;
4. uncomputes gen, then uncomputes each f_k from the output as c AND [R' == O_k]. This works because R -> R + A is a bijection.

M itself is total. For x = 0, Euclid emits an all-zero dialog and leaves u = q instead of 1. Bit 1 of u then distinguishes the two cases for free, and a controlled swap puts y back in place.

## Data shape

- **Wire.** An integer index. **Register.** A Python list of wires, least significant bit first. Shifts are relabels, implemented as list rotations with no gates.
- **Circuit.** A flat `array('i')` of (opcode, a, b, c), a free-list wire allocator that requires ancillas to be returned clean, and a high-water mark that becomes `num_wires`.
- **Block.** A recorded op range with input layout, output layout and internal wires. It can be instantiated forward or reversed onto any matching layout. That is how M^-1, Euclid-backward and every compute/uncompute pair are emitted without hand-written inverse code.
- **Simulator.** A bit-sliced simulator in which each wire is a Python int and bit j is test lane j. It shares no code with the generator's arithmetic.

## Why the permutation is enough for QNT-04

The circuit contains only X, CX, CCX and SWAP. It is therefore a permutation matrix with no phases. If every promise basis state |c, R, 0> maps to |c, R + cA, 0>, the same matrix acts exactly on any superposition of promise states and on any extension with a reference system. Basis exactness plus clean scratch is the whole unitary contract. Outside the promise, the circuit is still a permutation, so QNT-06's valid-channel requirement holds by construction.

## Qubit-time on PR #18

The sections below are the changes this branch makes to PR #18 (`b9b464a`) to lower qubit-time, the wire count times the board's scheduled ticks. Each one keeps the circuit's function. The CCX count is a budget to spend, as long as T stays under the reference's 50,465,933. Figures are as measured when each change went in. The final circuit is 1,546 wires × 1,920,902 ticks = 2,969,714,492, against PR #18's 1,555 × 2,404,728 = 3,739,352,040.

## Pebbled high increment

The Bézout tick sets the circuit's wire peak, at `P.write`, the tick's 221-bit high increment in `wide.PrefixTree`. Its block chains (165 wires), its AND tree (55) and its down carries (55) were all live at once, beside the multiplier's registers and the second low add.

`PrefixTree._write_pebbled` emits the gates of `write` in another order and on other wires, so the tree and the carries are never all live at once.

- Going down, a node clears its AND as soon as its carry is made. Only its parent's carry read that AND, and its children still hold theirs.
- Coming up on the complemented new bits, post-order, a node rebuilds its AND from its children's and then clears its carry. This is exact because c·AND(~new lo) = c·AND(old lo): with c = 1 the low half's new bits are all zero exactly when its old bits were all ones.
- P (`P_PEBBLE = 99`) clears its whole tree going down. M (`M_PEBBLE = -1`) keeps its down pass and only rebuilds node by node coming up.
- P may take `P_FRESH = 39` wires the tick has not touched, one fewer (`P_FRESH_EDGE`) on the ticks where the Bézout pass holds a second decoded block. Before taking one, it takes a free wire that is free by the time the gate's inputs are, or one that comes free at most `PEBBLE_WAIT = 10` ticks later.
- `PEBBLE_MID = 1` runs the tick's mid hook, CCX(t, z, e) then a drop of e, after the root's carry clear CCX(t, la, ch). Neither gate's target is among the other's wires, so they commute.
- `TOUCH_CAP = 318`: once a tick has touched 318 scratch wires outside its own registers, `_Idle._untouched` hands out the touched free wire that comes free first instead of a new one. Only the choice among clean free wires changes.

The tick's scratch fell from 327 to 314 wires with the pebbled write. The square's `ROW_TEMPS` and `FOLD_TEMPS` dropped from 30 and 126 to 24 and 120 to stay under it. On PR #18 the pebbled write took the circuit from 1,555 wires × 2,404,728 ticks to 1,542 × 2,383,186 at the same CCX. The fresh-wire tuning gave 1,540 × 2,384,510. On the merged branch the mid hook moved 2,088,411 ticks to 2,080,592, and the touch cap took 1,549 wires to 1,546 for 2,439 more ticks. `scripts/prefix_pebble_check.py` runs every write mode on every input to 13 bits and on drawn 221-bit inputs, with and without the tick's idle allocator, and catches every dropped-CCX mutant it makes.


## Tick latency

The Bézout tick is most of the circuit's ticks. These flags in `qec/two_bit.py` and `qec/wide.py` shorten it. Each one picks among exact constructions of the same map: block sizes, tree shapes, extra copies, wire choice, or a reorder of gates that commute.

- `WIDE_WAVE` with `WAVE_SEG = 64`: the add's first CX wave through U runs as a log-depth network in 64-bit segments, and p fans into the word through `P_FAN = 32` copies. The forward sweep's blocks then start together instead of one bit a tick, and the next tick's low bits do not wait for the compare to restore the high ones.
- `WIDE_ADD = 9`, `WIDE_ADD_RAMP = 1` and `WIDE_ADD_FSIZES`: the forward sweep's block sizes, from a local search on the tick alone (`scripts/tick_search.py`). Any split of a carry chain into blocks is exact. `WIDE_ADD_BACK` and `WIDE_ADD_BSIZES` size the backward chain, which only a tick with `chain_back` uses since the backward tree below.
- `WIDE_ADD_COPIES`: more copies of J, so the rotation's swaps do not queue on one copy. It was 24 here and is 26 since the CCX market below.
- `CMP_SPLIT = 14`: the compare's tree is a balanced tree over the bits below 14, joined at the root with a balanced tree over the rest. Any tree over the blocks is exact, and the low bits arrive last.
- `FINAL_EARLY = 1` runs the final swap's bit 0 on sw before sw's fan, since both only read sw. `FINAL_INTERLEAVE` deals the final swap's bits to sw's copies in turn.
- `ONE_BEFORE_DROP` takes one's wire before p's mark is dropped, so one is not on the wire the drop frees late. `IDLE_FROM_START` hands out the tick's wires by when they come free from its first gate. `IDLE_SLACK = -8` asks for wires free a little before the gate's inputs.
- `WIDE_ASIDE` with `ASIDE_UNDO_FIRST`: the trial chain runs on copies beside the second chain, and its undo is emitted first, so its first reverse join reads M's top before the second chain's skip does. `SECOND_FRESH` puts the second chain on wires the tick has not touched, and P's later write takes the trial's freed wires.
- `KUPD_FAST`: the k update is k' = k + e, and the second chain's corner selector is k0'k1' = k0k1 ⊕ e·k1 in the old bits. `three` already holds k0k1, so one CCX turns it into the selector and its wire becomes one. That saves a CCX and a step of the update.
- `MID_AFTER_CARRIES` moves `PrefixTree.write`'s mid hook after the carries' clear. Since P and M now both take the pebbled write, `PEBBLE_MID` governs the emitted circuit instead.

One tick alone, after the post-passes, went from 1,912 to 1,631 ticks. On PR #18 alone the circuit went from 2,404,728 to 2,192,192 ticks (−8.84%) for 31,406 more CCX. `scripts/tick_equiv.py` runs every tick variant the pass records under these flags and under PR #18's, on random basis states at full width. The outputs must agree lane for lane and the scratch must be clean wherever PR #18's is.

## Euclid step latency

The u6 shared difference, which runs the Euclid steps wider than `RIPPLE_MAX_W` whose ring fits it, is reworked in `qec/point_add.py`. Its CCX are unchanged except for the B1 marks.

- `U6_OWN_COPIES = 58`: the step's copies of g0 and g1 are on wires of their own instead of on the ring. The fans then wait for no tree temp, and no tree temp waits for a copy.
- `TWO_EARLY`: the step computes its second-strip flag `two` from V[1] at its start, and fans it onto up to `TWO_COPIES = 200` pass-level copies while the shared difference runs. This is exact because two = 1 implies g0 = 0, and the shared difference is the identity when g0 = 0. The pass alternates two wires for the flag, so a step does not wait for the previous step's uncompute.
- V[1] is read through a copy holding ~(b_1 ⊕ offset). The shared difference clears it after its up-sweep, where V[1] = a_1 ⊕ ~b_1 and U[0] = c_1 ⊕ a_1 with c_1 = ~g0, since a_0 = 1. The offset's XOR into V[0] moves inside the shared difference, after the X and CX on V[0], which it commutes with. A step then waits on neither the CCX that makes `two` nor the previous step's last offset read.
- `U6_KEEP_COPIES`: the own copies stay through the unwind, with no unfan and refan. G0 holds g0 ⊕ g1 during the final XORs, so each bit's V takes one XOR instead of two.
- `U6_RING_SELECT` and `U6_RING_PICK`: the sweep leaves its ring clean before the select, so the select also runs on ring copies of the copies, at most `U6_RING_COPIES = 80` a set. They go on the ring wires the sweep frees first under the tick model.
- `MARK_POOL = 8`: each B1 mark of the pass's down-sweep gets its own wire from a pool, made a few cells ahead from the two g0 bits and cleared after its B1 step. The B1 steps then chain only through the b wires, not through one mark wire too. This costs 1 CCX a cell.
- `ENCODE_SETS = 6`: the back-to-back block encodes rotate over 6 sets of clean wires instead of 3.
- `RIPPLE_MAX_W = 8`: with the early flag, u6 is faster than the in-place ripple from width 9.

On PR #18 alone this took 2,404,728 ticks to 2,302,536 (−4.25%) for 2,136 more CCX. The Euclid steps peak at 1,526 wires, below the Bézout tick.

## Post-passes

`scripts/post_pass.py` runs three passes on the emitted hierarchical circuit, in this order.

1. `qec/split_dirs.py`: every definition that runs both forward and inverted, and either holds calls or runs at least 8 times, gets a reversed twin. Every gate is X, CX, CCX or SWAP, each self-inverse, so the reversed list is the inverse. Inverted calls go to the twin, run forward, so each body runs in one direction only.
2. `qec/gate_pass.py`, unchanged except that it can be called as `apply()`.
3. `qec/reorder_pass.py`: for each body that holds only gates, a local search over legal reorders of commuting X-type gates along the body's critical path, plus CCX control swaps. Two such gates commute when neither's target is among the other's controls, and a move is legal only when every gate it passes that shares a wire with it commutes with it. CCX is symmetric in its controls. The objective is the body's depth under the board's tick model. Bodies run at least 50 times also get seeded rounds of iterated local search, seeded by the definition index, so the result is reproducible. The board scores the expanded circuit, so new bodies go in one caller group at a time, and a group is kept only if the whole circuit's ticks do not rise.

On PR #18 alone the passes cut 4.87% of the ticks. On the final circuit they take 2,000,528 ticks to 1,920,902 (−3.98%) at the same wires and CCX. `scripts/check_passes.py` runs the three passes on every toy circuit that `scripts/check.py` runs, expands the result and checks every lane.

## CCX market

PR #18 sat 13,610 CCX under the reference. This change spends that room where a CCX buys the most ticks and sells CCX where they cost the fewest. The rates below are ticks per CCX on the whole circuit.

- **Buy: the tick add's backward sweep as a tree** (`LA_TREE_B = 1` in `qec/two_bit.py`). The backward chain stepped 11 ticks a block down the word. It is now `arith`'s Brent-Kung lookahead tree in blocks of 1, run in reverse, with the second wave through J's copies. The forward sweep F maps (x, V) to (f, V) with f_(s+1) = x_(s+1) ⊕ V_s·f_s. The W_j writes change V only where f_s = 0, so F's map on the new V still sends x to f. Any exact sweep of that map with clean temps, run backwards, sends f back to x. This saves 1.44 ticks a CCX. `WIDE_ADD_FSIZES` was retuned for it, 42 blocks and 21 CCX fewer a tick.
- **Free: more copies for the rotation** (`ROT_EXTRA = 102`). Once the sweeps' temps are clean, the rotation fans J onto 102 of them and undoes the fan after. With the tree, the rotation no longer hides beside the chain. `U6_RING_PICK` above is also part of this change.
- **Sell: the square's rows as chains** (`qec/square.py`). A chain has fewer CCX and more ticks than a tree. `ROW_TREE_MIN = (126, 126, 16)` makes the halves' rows chains except their top two, and the halves' sum's rows below 16 bits chains. `SPLIT_TREE = (True, False, False)` makes the takes into K and K's add into Q chains. Each is an existing exact add, chosen per add. The halves' rows cost 0.92 ticks a CCX, the takes 0.24, and K's add nothing. Turning every row to a chain (`ROW_TREE = False`) would cost 1.46, more than the tree buys.

Measured after emit and gate_pass, this took 1,549 wires × 2,088,411 ticks = 3,234,948,639 to 1,549 × 2,005,858 = 3,107,074,042 (−3.95%). CCX went from 7,195,809 to 7,208,765, so T is 50,461,355, 4,578 under the reference's 50,465,933. `TwoBitMod.chain_back` and `point_add.CHAIN_EVERY` let some ticks keep the backward chain to set how many CCX the tree takes. They are checked on toys and unused (`CHAIN_EVERY = 0`).
