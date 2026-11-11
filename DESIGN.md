# Exact controlled point addition for `secp256k1-controlled-exact-v1`

## Current mux contract

The merged path now uses the paired-tick contract D1-D4 in [the mux proof](docs/mux-proof.md), superseding the older unrestricted canonical per-symbol Bézout-cell claims below. Logical T and R are canonical and scratch starts clean, but after a non-final B1 the physical word stores H = 2T + R in S4's encoded form, borrowing the next cell's g0 for its high bit. A semantic subtract follows B1; next.g0 must be 1 before holding; a final B1 requires R = 0. The held tick restores the borrowed symbol. Fallback primes retain the old step.

The additional sweep rotates each no-subtract run so B1 is last in Bézout order.
It runs in Euclid's tail, one cell behind the offset sweep, with which it commutes: the offset sweep's swap acts only on b of g0 = 1 cells, and each B1 step only on b of two g0 = 0 cells.
So Euclid hands Bézout the B1-swept tape, and inverse Euclid unsweeps it before undoing its steps; this supersedes the placement in S2's last sentence and S17's last sentence of the mux proof, not the maps they state.
The tape at rest is block-coded, two wires freed per 15 ticks by `qec/tape_block.py`; on the swept tape every B1 follows a subtract cell in tape order, so the code reads each block's cells last first and relies on the same clause.
Cells read previous-symbol wires and can write the next g0 across a block edge, so Bézout decodes a block one tick before its first tick and codes it again one tick after its last, and two blocks are decoded across each edge for two ticks.
The gate's toy primes have at most one block, so the edge path is covered by `scripts/block_check.py` at 20 bits and above and by the secp256k1 rows, not by a toy instance.
The four-symbol alphabet, 384-tick bound and width schedule are unchanged. M and its inverse use clean states in their respective input/output layouts; top-level R remains on-curve or infinity with canonical encoding and zero scratch.

The corrected S13 distinguishes logical carries from physical borrowed first-probe wires: only on P = 1 do those wires equal the carries; on P = 0 their junk is unobserved except through P and reversed. Exported e and the second chain remain zero on the inert lane. Credit arXiv:2609.05625 section VI for the offset representation.

The exact-source independent review covers this pair contract; it does not grant upstream checker acceptance. The local held secp BDD probe stopped at the production 4,194,304-node budget. The earlier cell's cap-transfer proof does not automatically transfer. I7.charge, qualification, and the flat-upload operation cap remain open limits.

## Clean wires in the mux tick

Before these changes the mux tick held the circuit's peak at 1,275 wires, with 26 clean, while Euclid, the code calls and the B1 step peaked at 1,249.
So each clean wire the tick gives up lowers Q by one, until another stage binds.
A wire is worth T / (768 × 7 × Q) CCX per tick, counting 768 tick uses and 7 T per CCX, which is 3.40 on this tree.

Two changes lower the tick to 19 clean wires at secp256k1, for 1,268 wires × 23,151,240 T on this tree, against 1,275 × 23,075,976 before them.
Their statements are C1–C4 and R1–R4 in `INVARIANTS.md`.

- **Kept low carries.** The second chain no longer keeps every low carry through its top.
  It keeps c2, c6, c9 and c11 on wires and c12 on bit 11's selector, clears each other one once the next kept carry is built, and rebuilds them on the way down, at 2 CCX per dropped carry.
  The mux cell already drops c1 to pay for k2, so this saves 6 more wires, at +12 CCX per tick.
  At secp256k1 a pebble search over these carries, whose up pass may clear and rebuild any carry many times, finds this peak at the same cost, one wire lower at +6 CCX more, which is above a wire's worth, and nothing two wires lower.
- **The quiet run's prefix top.** M's top below its last is read only by `M.write`, so it is cleared right after M is built and recomputed on a new wire just before the write.
  That is 1 wire at +2 CCX per tick.

Three other wire trades were priced and set aside:

- `l0` cannot be recomputed across the second chain from g0 and R₀.
  On a non-last B1 it holds H's bit n, parked there by the hide, and that bit is on no other wire.
- Moving `Prefix` tops onto lower wires costs about +6 CCX per wire.
  Each intermediate top is also the next block's `a`, read once when the prefix is built and twice in `write`.
  A conditionally clean copy needs a true control beside it at every read, as in arXiv:2407.17966, so every read pays a gate.
- Releasing P's top below its last saves nothing, because `P.write` runs at the second chain's top, which is the tick's peak.


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
3. runs the generic kernel controlled only through gen, using M^-1 ... M. With gen = 0, every gen-controlled step is the identity, and M^-1 followed by M cancels exactly at the circuit level when M's entries are in [1, q), which the kernel's flag code ensures;
4. uncomputes gen, then uncomputes each f_k from the output as c AND [R' == O_k]. This works because R -> R + A is a bijection.

M's standalone contract is narrowed to canonical x in [1, q); x != 0 at M's entries rests on the kernel's flag code.
The build refuses curves with a finite point at x = 0 and addends with x(A) = 0 or x(2A) = 0.
The kernel codes both the -A exception and the c = 0, R = inf marker as f0 = f1 = 1, while keeping gen zero.
It flips X's low bit under f0 AND f1 before M inverse and restores it after M, holding no extra wire during M.
The direct evidence is the instrumented observation that the indicator equals [X = 0] before the remap and X != 0 at both M entries, on all 25,800 toy lanes and the secp256k1 special points.
The guarantee's control of record removes both remaps with every flag intact; the existing point gates then fail on 72 lanes, c = 0 with R = inf and c = 1 with R = -A, across all 36 toy instances.
The marker-removal, one-remap and old-gen-fold mutants are sensitivity controls; they do not test the nonzero-X guarantee.
This named limit accompanies I2.model's narrowed zero-caller row.

## Reusing the Euclid pass's zero words

These arguments apply to every odd prime q of n >= 2 bits and every canonical x and y in [0, q).
Let P = q(q - 1), T = bitlen(P), and x' = x XOR [x = 0].
Then 1 <= x' < q, including x = 0 and x = 1, whose distinct original values are retained by z.

**F1: the terminal state is (1, 0).**
Initially u = q is positive and odd, v = x' is nonnegative, and gcd(u, v) = 1 by primality.
If v is even, halving v halves the product and preserves the gcd because u is odd.
If v is odd, the conditional swap makes the subtracted operand no larger than v, and (v - u)/2 is nonnegative.
That step preserves the gcd and makes the product at most half its preceding value, with u still positive and odd.
Thus u_i * v_i <= floor(P / 2^i), and at i = T the product is zero.
Since u_T is positive, v_T = 0; the preserved gcd then forces u_T = 1.
This is an argument using the existing iteration bound, not a sampled reduction of that bound.

**F2: the zero remap is reversible.**
The flag z stays allocated while x' is used, and the inverse pass restores x' before the remap is undone.
The resulting x is the original input, so the repeated zero test clears z exactly.
The x' = 1 path multiplies by one and therefore implements the specified identity on y for original x = 0.

**F3: V supplies S.**
V is zero after the pass, even without remapping zero, so its wires can hold the initially zero Bezout word S.
The existing exact Bezout recurrence on canonical R = y and S = 0 ends at R = 0 and S = x' * y mod q.
No arithmetic cell, modular reduction or transcript alphabet changes.

**F4: only U's upper wires are lent.**
With the remap, U = 1, so U[1:] are zero and U[0] must remain allocated.
The allocator can lend those n - 1 zero wires to Bezout scratch; the unchanged exact cells restore that scratch to zero.
Reclaiming n - 1 zero wires afterwards need not return the original wires or their order.
The reclaimed word consists of the still-live U[0] followed by these zero wires, disjoint from the live x and y words.

**F5: restore the terminal state on a new layout.**
X on U[0] zeros the reclaimed word, and swapping it with S moves the product onto the word that becomes out while zeroing V again.
R on y's wires is zero, so X on its low wire recreates one there.
The triple (y wires, V wires, tape) now contains (1, 0, tape), the terminal state required by the inverse Euclid block.
Block replay on that matching layout restores (q, x', zero tape), and loading q again clears y's wires.
The cleared y word and tape are freed before the final equality test, allowing its AND chain to reuse their wires without increasing the peak.
The remap is then undone, z is cleared and freed, and the product word becomes the logical y output.
Every other wire is returned clean; reversing this same block gives M inverse, including its identity at x = 0.

This change relies on the existing proven bound and width schedule and changes neither, nor any modular reduction.
For nonzero remapped inputs, active operands fit the existing schedule, and the idle state (1, 0) also fits every emitted step.
The gate-level small-prime checks and 256-bit boundary components exercise the implementation; they do not replace the all-input argument or constitute acceptance by the separate benchmark proof checker.

## Shared-chain Euclid pass with a held offset

For every odd prime q, put n = bitlen(q), P = q(q - 1), T = bitlen(P), and use the unchanged width_schedule.
Let (u_i, v_i) be the reference Euclid state and (b_i, s_i) its parity and conditional-swap bits.
The following describes the step and sweep introduced with this proof; the existing F1-F5 and C1-C10 arguments remain its external boundaries.

**E1: the physical step, before wire rotation.**
For w >= 1, take a w-bit odd word a, a w-bit word B, prior offset bit e (zero for the first step), two zero tape bits and a zero carry-in.
Let vhat = B XOR e on the low bit, S = a + NOT(vhat) modulo 2^w, and c = [a > vhat], its carry out.
The physical output before interpreting the rotation as halving is:

| Condition | a output | B output | g0 | g1 |
|---|---|---|---|---|
| vhat even | a | vhat | 0 | e |
| vhat odd, c = 0 | a | vhat - a | 1 | 0 |
| vhat odd, c = 1 | vhat | a - vhat - 2 | 1 | 1 |

The upward majority chain stores c_(i+1) on a[i], NOT(vhat_i) XOR a_i on B[i], and c_i XOR a_i on the preceding carry wire.
On descent, one CCX restores a_i and a second adds g0*(c_i XOR a_i) to B[i].
Two CX gates restore the preceding carry and remove a_i from B[i], leaving S_i when g0=1 and NOT(vhat_i) otherwise.
The identity S_i XOR c_i XOR 1 = a_i XOR vhat_i lets an X-conjugated CCX controlled by g1 exchange the selected a bit without a second ripple.
Complementing B gives NOT(S) when g0=1 and vhat otherwise.
A CX from g1 to each upper B bit selects S on those bits in the swap case, while retaining NOT(S_0) on the low bit.
Odd a and odd vhat imply S_0=1, so B's low bit is zero in both odd branches; it is also zero in the even branch.
The final controlled copy adds (NOT(g0))*e to g1; it preserves e, and the carry-in returns to zero.
Thus the low-bit rotation of the full V word is exact halving, and E1's table must not be confused with the logical post-rotation word.
The first variant costs 4w+1 CCX and later variants cost 4w+2, including the offset copy.

**E2: held-state induction.**
For 1 <= x < q, before each step U=u_i and V+e_i=v_i with V>=0, e_0=0, and e_i the preceding cell's second bit for i>=1.
If e=0 the absorb leaves V=v; if e=1 and v is odd it changes the even V=v-1 to v; if e=1 and v is even it changes the odd V=v-1 to v-2.
Therefore vhat has v's parity and is at most v.
For odd v with u>v, the gate outputs after rotation are (U,V,e')=(v,(u-v)/2-1,1); u-v>=2 makes V nonnegative.
For odd v with u<=v they are (u,(v-u)/2,0).
For even v they are (u,v/2,0) if e=0, and (u,v/2-1,1) if e=1.
Each case preserves V+e'=v_(i+1) and U=u_(i+1), with e' the current second tape bit.
During active steps u and v fit the existing window, so vhat<=v and V<=v make this chain exact in that window.
This uses the existing product-halving bound proved in F1; it is not Bernstein-Yang's different divsteps bound.

**E3: zero tail and terminal offset.**
When v_i=0, nonnegative V and e force V=e=0.
Then g0=g1=0 and the upward/downward chain restores U and V for every window, including one narrower than U.
F1 gives u_T=1 and v_T=0 for nonzero canonical x over prime q, hence U=1,V=e=0 at the end.
Primality is used here through F1's gcd argument.
The last odd step has u=v=1, since reaching zero requires equal operands and their gcd is 1, so that step leaves no offset.

**E4: descending tape sweep.**
Write z_i for the unswept second bit: z_i=s_i on odd cells, z_i=z_(i-1) on later even cells, and z_0=s_0.
For i=T-1 down to 1, XOR (NOT(b_i))*z_(i-1) into z_i.
Descending order reads the predecessor before it is changed, so even cells become zero and odd cells keep s_i.
The swept tape is exactly the reference dialog, allowing the unchanged Bezout cell, reduction proof and F1-F5 cleanup to apply.
The first and last offset exchanges now use CX-only endpoint identities, and the merged sweep skips the terminal B1 step and its mark transition.
The sweep remains inline in the pass body so tape blocks can be decoded and coded between steps.

**E5: zero inputs and generic callers.**
The multiplier remaps original x=0 to x'=1, covered by E2; the retained flag distinguishes it from original x=1.
The pass's domain is therefore V != 0, which E7 needs.
A generic caller with V=0 would keep U=q, whose high wires are not zero, so E7's lent wires would not start clean; the pass no longer supports that caller.

**E7: late tape on U's vacated wires.**
For V != 0, the proven bound gives max(u_t, v_t) < 2^(w_t) at every boundary t, or (u_t, v_t) = (1, 0), so U < 2^(w_t).
A step at width w touches only U[:w], and widths never grow, so U's wires from w_t up are zero from boundary t to the end of the pass.
When the width drops, the pass adds those wires to the tape's supply, the highest first, as it already does with V's.
It lends as many as it ends with zero wires: V's tail and two free wires per coded block, 52 at secp256k1.
At the end U takes those zero wires back as its upper wires, so U is again (1, 0, ..., 0) and S reuses the whole word.
The boundary is a permutation of U and the pool together, so replaying the pass block in reverse gives inverse Euclid the mirror layout with no extra gate.
Euclid's pool is smaller by the lent wires, which puts Euclid's peak below Bezout's, so no block is coded during the pass any more.

**E6: reversal and unchanged surfaces.**
Reversing the pass reverses its sweep and then its step calls, including their layout permutations.
F5 restores precisely the terminal words and swept tape needed by that inverse, so no additional garbage or cleanup assumption is introduced.
The external pass relation, iteration bound, width schedule and modular reduction are unchanged; EuclidStep's physical relation and internal second-bit meaning change as E1 states.

At G, sum(W)=98432 and T=512, so each pass costs 4*98432+2*512-1+511=395262 CCX instead of 591104.
Four passes save 783368 CCX, or 5483576 T; these are formulas to check against emitted counts.
This change leaves _bezout and qec/arith.py unchanged.
Finite gate checks and the separate leaf lift probe do not establish acceptance by the benchmark checker or compose this proof automatically.

## Exact merged Bézout reduction

For every odd q >= 3 and canonical R, S, `Mod.dbl_add` computes (2S + g0 R) mod q with clean scratch. The Bézout step then applies the existing g1-controlled swap; its contract and [2, n, n] layout are unchanged.

q is odd with n = bitlen(q), and f = 2^n − q, so f is odd and 2^n < 2q.
m is any integer with bitlen(2f) ≤ m ≤ n, so 2f < 2^m, and the generator takes the smallest.
That range is never empty, because 2f < 2^n.
At secp256k1, m = 34 and h = n − m = 222.

R, S < q and g0 is a bit.
s = 2S + g0·R, so s ≤ 3q − 3.
t = s >> n, low = s mod 2^n, L = low mod 2^m, H = low >> m, and A = [H = 2^h − 1], with A = 1 when h = 0.
For j in {0, 1, 2}, c_j is the carry out of L + j·f on m bits, so c_0 = 0 and c_1 ≤ c_2 ≤ 1.
e = A·c_(t+1) for t ≤ 1, e = 0 for t = 2, and k = t + e.

### Gate order

Every wire below is a named position in `s`'s bits: bit 0 is a fresh wire `l0`, bits 1 to n are S's wires, and bit n + 1 is a fresh wire `top`.
So `low` is `[l0] + S[:n-1]`, t0 is `S[n-1]`, and t1 is `top`.

| # | Action | CCX at n = 256 |
|---:|---|---:|
| 1 | `add(R[1:] + [l0], S, top, ctrl=g0)`, with `l0` still zero as a's top bit | 3n + 1 = 769 |
| 2 | l0 ^= g0·R[0] | 1 |
| 3 | Prefix-AND chain over H, the carries of H + 1, giving A_i and A_h = A | h − 1 = 221 |
| 4 | Set selectors in place: X t1, CX t0 → t1, so t1's wire is [t = 0] and t0's wire is [t = 1] | 0 |
| 5 | Carry chain up on L + [t = 0]·f + [t = 1]·2f, giving c at the top | m = 34 |
| 6 | Fresh `e` ^= A·c, a CX when h = 0 | 1 |
| 7 | The same chain down, then undo step 4 | m = 34 |
| 8 | t += e in place: CCX(t0, e, t1), CX e → t0. The wires now hold k0 and k1 | 1 |
| 9 | Carry chain up on L + k0·f + k1·2f, giving c_k at the top | m = 34 |
| 10 | e ^= A·c_k, which clears `e`, a CX when h = 0 | 1 |
| 11 | For i from h − 1 down to 1: clear A_(i+1) from A_i and the old H_i, then H_i ^= c_k·A_i. Then H_0 ^= c_k, a CX | 2(h − 1) = 442 |
| 12 | For i from m − 1 down to 1: clear c_(i+1) from the old L_i and c_i, then L_i ^= addend_i xor c_i. Then clear c_1 from `l0` and the bit-0 addend, leaving `l0` alone | m = 34 |
| 13 | Fresh z = g0·¬k0 | 1 |
| 14 | `gt_into(R[1:], low[1:], t1, ctrl=z)`, which clears k1 | 2(n − 1) + 1 = 511 |
| 15 | Clear z | 1 |
| 16 | CX l0 → t0, so t0's wire holds S'[0], then l0 ^= g0·R[0], which clears `l0` | 1 |
| 17 | Free `l0`, `top`, `z`, then rotate `S` by one wire, as `moddbl` does | 0 |
| 18 | `cswap(g1, R, S)`, in `_bezout` | n = 256 |

A low chain has one AND per position, with the addend bit at position i being one of 0, the first selector, the second selector, or their XOR.
That is MAJ(L_i, addend, carry), and CX conjugation gives it with one CCX.
`_carries(C, H, 1)` already builds step 3, up to A_(h−1), and A_h is one more carry on top of it, so reuse what fits.
`_reduce`'s down loop has step 11's order: it uncarries position i + 1 before it writes position i.

**Cleanup order in steps 11 and 12.** A stored wire reads the old bit below it, so each one clears before that bit is written.
A_h goes first, right after step 10, because it is used only by steps 6 and 10 and writing H_(h−1) first would leave it dirty. With h = 2, H = 11 and c_k = 1, A_2 would stay 1.
A_1 is an alias of H_0 and has no wire to clear.
With h = 1, A is H_0, step 3 is empty, steps 6 and 10 read H_0, and step 11 is the single CX.
With h = 0, A is the constant 1, step 3 is empty, steps 6 and 10 are CX from c and c_k, and step 11 is empty.
Step 7 runs step 5 in reverse with nothing written, so it has no such order to get wrong.

The cell is inline, with no `C.apply` or `C.record` inside, because `lift` refuses a leaf with calls.

### Correctness and cleanup

**C1.** t is 0, 1 or 2, because s < 3q < 3·2^n.

**C2.** For j ≤ 2, low + j·f ≥ 2^n exactly when A·c_j = 1.
This holds because low + j·f = (H + c_j)·2^m + ((L + j·f) mod 2^m), and H + c_j reaches 2^h only from H = 2^h − 1 with c_j = 1.
The carry is one bit because 2f < 2^m.

**C3.** k = floor(s/q), and S' = s − kq = low + k·f − e·2^n.
With e = 0 and t ≤ 1, C2 gives low + (t + 1)f < 2^n, so 0 ≤ low + t·f < q.
With e = 1, C2 gives low + (t + 1)f ≥ 2^n, so S' ≥ 0, and S' < q is the same as s < (t + 2)q, which holds for t = 0 because s < 2^n < 2q, and for t = 1 because s < 3q.
With t = 2, S' = low + 2f < q is the same as s < 3q.

**C4.** Step 5's addend is (t + 1)·f below 2^m for t ≤ 1 and 0 for t = 2, so step 6 computes the e above.
Step 7 finds the selectors as step 5 left them, because step 6 touches only `e`.

**C5.** e = A·c_k, on every reachable state.
If e = 1 then c_k = c_(t+1) and A = 1.
If e = 0 and t ≤ 1 then c_k = c_t ≤ c_(t+1), so A·c_k ≤ A·c_(t+1) = 0.
If t = 2 then low + 2f < q < 2^n, so A·c_2 = 0 by C2.
So step 10 clears `e` while A and c_k are still live.
This is the sketch's trap: `e` is the overflow out of H as well as k − t, so it needs no clear of its own from the final state.

**C6.** S' = ((H + c_k) mod 2^h)·2^m + ((L + k·f) mod 2^m), which follows from C3, and S' < 2^n.
Steps 11 and 12 write exactly that, and k0, k1 are exclusive because k ≤ 2.
Every stored A_(i+1) and c_(i+1) reads 0 after its clear, because each clear runs on the same inputs that set the wire: the old bit, and a lower chain wire that has not been cleared yet.

**C7.** k1 = 1 needs g0 = 1, because g0 = 0 gives s = 2S < 2q.

**C8.** Step 14 clears k1.
If z = 1 then k is even, so S' and R have the same parity and [R > S'] = [R >> 1 > S' >> 1].
With k = 0, S' = 2S + R ≥ R.
With k = 2, S' = R − 2(q − S) ≤ R − 2.
So the compare equals k1.
If z = 0 then k1 = 0, by C7 or because k0 = 1.
k = 1 can have R > S', so the compare needs z and not g0 alone.

**C9.** S'[0] = k0 xor g0·R[0], because 2S is even and q is odd.
So step 16 leaves S' on the rotated wires and clears `l0`.

**C10.** For each g, S → 2S + g0·R mod q is a bijection on [0, q), because q is odd.
The circuit is a permutation that takes canonical states with zero scratch onto canonical states with zero scratch, so the reversed cell that M⁻¹ runs is exact on canonical inputs too.

Nothing here uses primality, n = 256, or a run of ones.
The run of ones only makes m small.

## Data shape

- **Wire.** An integer index. **Register.** A Python list of wires, least significant bit first. Shifts are relabels, implemented as list rotations with no gates.
- **Circuit.** A flat `array('i')` of (opcode, a, b, c), a free-list wire allocator that requires ancillas to be returned clean, and a high-water mark that becomes `num_wires`.
- **Block.** A recorded op range with input layout, output layout and internal wires. It can be instantiated forward or reversed onto any matching layout. That is how M^-1, Euclid-backward and every compute/uncompute pair are emitted without hand-written inverse code.
- **Simulator.** A bit-sliced simulator in which each wire is a Python int and bit j is test lane j. It shares no code with the generator's arithmetic.

## Why the permutation is enough for QNT-04

The circuit contains only X, CX, CCX and SWAP. It is therefore a permutation matrix with no phases. The point-addition promise is R on the curve or infinity, with canonical encoding and zero scratch. On this domain, x and y parity identify each finite point; the x=0 guard retains a full y comparison when required to distinguish infinity. Both flag computations use this same domain before and after the controlled translation. If every promise basis state |c, R, 0> maps to |c, R + cA, 0>, the same matrix acts exactly on any superposition of promise states and on any extension with a reference system. Basis exactness plus clean scratch is the whole unitary contract. Outside the promise, the circuit is still a permutation, so QNT-06's valid-channel requirement holds by construction.

## Last-cell lend

The merged pass lends the uncoded final tape cell's b after its last read, using its proved complement relation with g0.
Restoration allocates fresh zero storage and updates the tape boundary before inverse Euclid.
The codec and tick schedule do not change.
See `docs/tail-lend-proof.md`.

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

### Euclid tick 0 on its initial boundary

The initial U is q, V is the nonzero-remapped x, E is zero, and the cell and remaining tape are clean.
The first step alone uses `_euclid_tick0` when the width is n and k = n - 1 - lsb(2^n - q + 1) unwritten tape wires are available.
For odd v its output is U=v, V=(q-v)/2-1, E=1, cell=(1,0).
For v=2 mod4 it leaves U=q, E=0 and cell=(0,0), with V=v/2; for v=0 mod4, V=v/4 and cell=(0,1).
The V list rotates exactly as in the general step.
Controlled addition of 2^n-q+1 uses the borrowed wires as carries; it uncomputes them before returning.
Its first carry uses the preceding select identity through the optional low callback rather than another controlled gate.
The recorded block's input/output wire sets agree and contain every used wire.
If width, supply or that shape contract fails, the first step falls back to `_euclid_step`; later steps always use the general block.
The borrowed wires are zero at tick 0's entry/exit and the mirror points in the inverse pass.
They may hold nonzero tape at E's exit, so they are not generally clean scratch there.

### Inline kernel constants and reflection

For cmodadd_const, k=c*[b>=q-v].
The control temporarily holds c XOR k, so selectors (c XOR k,k) exclusively select v or v+2^n-q.
The selected carry chain starts at the lowest set bit of either mask and uncomputes each carry before overwriting its old input bit.
After restoring c, [b'<v] clears k.
For creflect, k=c*[b>a].
Complement b under c and add a+1 or (a+1+q) mod2^n under the same exclusive selectors.
The output remains canonical and lies above a exactly when k, so CX(c,k) followed by comparison a+1>b' clears the flag.
This reflection is valid at b=0 as well.
It replaces kernel cneg and the final constant addition; the intervening subtraction acts on Yr and commutes.
All gates are inline: addend-dependent constants introduce no new recorded definitions.
Existing standalone cneg behavior and its tests remain.

### Kernel flag recode across M

At either M boundary, (c, f0, f1, gen) is one of 0000, 0110, 1000, 1100, 1010, 1110 and 1001.
With c = 0 the word is 0000, or 0110 with the zero marker; with c = 1 it is 1001 on the generic path and 1000, 1100, 1010 or 1110 on an exception.
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

For canonical x < q, strips before the first subtract leave U = q > V, so that subtract swaps and the swept block-zero word lies in the 4,022-word D1 subset of the L6 domain.
On the merged path with a block, a 21-CCX pre-map followed by the existing codec frees W10 and W13, with seven pre-map ancillas returned before the codec allocates its two.
The extra free wire joins the clean pool that U borrows at the pass boundary, reducing external tape by one; with_free_slots restores both positions for later passes.
A one-CCX HOLD frees b0 while ticks 7 through 2 run and restores it before tick 1, without touching the g1 next-cell lend.
Nonmerged paths keep the original codec, and no new per-prime domain refusal is needed.
The exhaustive D1 component row checks encode and HOLD both ways, clean FREE and ancillas, and distinct images, increasing the census only to 52 rows and 13,867 lanes.

### Holding block zero across the inverse pass

The merged multiplier records the Euclid pass with block zero raw and codes that block immediately after the recorded call.
The default pass remains the raw variant followed by the same block-zero encode, with identical operations and output wire lists.
The Bezout pass leaves block zero decoded at its end only when the multiplier requests this variant.
Before the inverse Euclid pass, the multiplier places the block's two free-slot wires in the raw pass's output positions and releases the zero wires previously occupying them.
The default pass and every existing gate-script caller retain their encoded boundary.
This removes the adjacent encode and inverse decode in both multiplier calls, saving 204 CCX and one listed call without increasing the 1,246-wire peak.
The existing breakdown script still measures default passes; its rows include the removed encode/decode work and therefore do not sum to the optimized multiplier's cost.

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

## Uncontrolled top add bit

Without carry output or control, the top majority and unmajority reduce to b XOR a XOR the incoming carry.
Two CX gates replace that pair while every lower bit retains the original chain.
Controlled and carry-output calls keep their existing implementation.
The full circuit loses 44 CCX.

## Bounded square digit adds

Each digit bounds the accumulator before and after its addition using the corners of the caller's input box.
A subtract bounds the complemented value, -v - 1.
Above the first bit containing every bounded value, the remaining bits hold the sign.
A known zero sign allows truncation; a variable sign uses the add's carry output and copies it above that bit before and after the add.
A positive digit added to an initially zero accumulator with no carry input copies its source.
The first fold accepts arbitrary n-bit source and accumulator values, and the second accepts a source at most f and a zero accumulator.
This saves 984 CCX after the top-bit edit; 20 CCX overlap that edit's standalone saving.

## Repeated-width strip boundary

A Euclid tick leaves V below the top bit of that tick's width.
When the preceding tick used the same width and the current width is at least three, the second strip's final swap exchanges two zero bits.
The pass marks precisely those calls with top_clear, and recorded shapes cache both width and that flag.
Unflagged callers retain the full swap sequence.
The full circuit loses 516 CCX.

## Sweep endpoint identities

The first cell has U = q greater than x, so a subtract swaps and g0 implies E XOR b equals one.
At the last cell, no subtract implies E XOR b equals one, giving the other CX-only offset exchange.
For a merged sweep with at least three cells, the last two no-subtract cells have equal b bits.
Their B1 transition never swaps, so the sweep starts its carried mark at T - 2 and skips that terminal B1 step and mark transition.
Together the endpoint and terminal changes remove four CCX per Euclid pass and sixteen per point addition.

## First carry in tick zero

The tick-zero select leaves U[j] = q_j XOR g0 times (v_j XOR q_j).
For the lowest set bit j of the following constant, the first carry g0 times v_j is U[j], complemented with g0 and one when q_j is one.
An optional low callback supplies that carry during computation and uncomputation; every other adder caller retains the default carry primitive.
The full circuit loses eight CCX and keeps the same clean-wire contract.

## Tail-code control substitutions

On the checked tail-code domain, primitives 103 and 109 use their single remaining controls.
The two CCX-to-CX substitutions preserve the code's complete forward and inverse domain contract.
The primitive list retains its length, and its header digest records the substituted list.
The full circuit loses eight CCX.

## Single-Toffoli kernel recode

After the initial flag CX gates, the c = 0 words are 0000, 0110 and 0101.
None has f0 = 0 and f1 = 1, so f0 XOR f1 equals f0 AND NOT f1 on that domain.
The recode toggles gen with NOT c AND (f0 XOR f1), temporarily forming and restoring the XOR.
This uses one CCX without a temporary, saving eight CCX across the four recodes.

## Kept bottom carry for bit one

This branch builds on ab-godder's `mike/comp-ab-godder-pr9-286e36a3ba9f` (`lowt-mike` at `db9ae43`) and adds one edit, in `qec/two_bit.py`.
`_chain_down` now takes `keep = 1` from the add's second chain, so that chain's bottom carry c1 = l0 AND k0 stays built past the compare instead of being cleared.
On a held lane l0 = R0, so c1 is the k0 AND R0 term of bit one, and bit one reads it from c1 instead of computing and uncomputing it on a fresh wire around its use.
c1 is cleared with one CCX, l0 AND k0, after the B1 drop and before k0 changes, and its wire is freed.
Every tick call with a previous cell saves 2 CCX, 764 of the 766 calls at secp256k1, so the circuit at G loses 1,528 CCX at the same 3,135 wires.

## Held offsets beside the tape

Each Euclid cell now has a third wire, held, zero on entry; the pass takes them as a register of `iters` wires and only the inverse pass clears them.
A tick computes two = NOT g0 AND NOT (V0 XOR E) onto held, clears V0 under it and strips on it as before.
Then b ^= held makes b = s XOR two, which is s on a subtract cell and two otherwise, since s implies g0 and two implies NOT g0.
held ^= b leaves held = s, and cswap(g0, E, held) moves the old offset e onto held and s onto E on a subtract cell; on a no-subtract cell s = 0, so held is 0 and E keeps e.
So E follows E2's recurrence unchanged, held ends as g0 AND e, and each subtract cell already carries its own swap: the tape after the ticks is the offset-swept tape of E4, and the offset sweep is gone; the B1 sweep is unchanged, because it reads only g0 and the b of no-subtract cells, which are the same.
E ends zero by E3, as before.
Tick 0 has e = 0 and s = g0 (U = q > v), so it uses no held wire and ends with b ^= g0.
The last tick maps (1, 0, 0), (1, 1, 0), (1, 0, 1) to cell (V0 XOR E, NOT g0), held = E, and zero V and E, by CX.
The inverse pass is the recorded block reversed on the same held register, so held returns to zero; M frees it after that pass.
Per tick this drops the CCX that cleared the old two wire, and per pass the 382 offset-sweep CCX; at G the circuit loses 3,056 CCX for 384 held wires, which live inside the multiplier, below the square's 3,135-wire peak.
`scripts/model_check.py` compares held with the swap of the last subtract before each subtract cell, from the semantic tape.

## Carries kept for the inverse pass

A Euclid step's chain may build each carry c_(i+1) = a ^ (a ^ b)(a ^ c) on a fresh zero wire with one CCX, restoring a, b and c, instead of the in-place majority.
U[i] then still holds a_i on the way down, so each bit's descent is the g0 select of S_i and the g1 exchange, two CCX instead of three; bit 0 and the top bit are as before, the top bit reading c_(w-1) from its kept wire.
The step's U, V, E and cell outputs are the same functions of its input as before; the kept wires end holding the carries of a + NOT(vhat) on that input, so the forward pass's extra output is a function of its input.
The multiplier allocates these wires zero beside held, the Bezout pass never touches them, and the inverse pass, the recorded block reversed on the same wires, clears them, so M returns them zero.
Each kept carry saves one CCX in each of the four passes.
`KEEP` bounds the kept wires; the plan keeps every chain carry of narrow ticks, from the last tick back, while the total fits, and never at full width.
At G, KEEP = 1480 keeps 1,431 carries, ticks of width 55 down to 3, and the circuit loses 5,724 CCX at the same 3,135 wires: the multiplier's peak stays below the square's.
Once the square sits at the 4,096-wire cap (see the square wire trades), the multiplier may grow to the same peak, which is 1,697 wires plus its kept carries.
KEEP = 2400 then keeps the whole ticks of width 70 down to 3, 2,346 carries, for a multiplier peak of 4,043 wires: 3,187,581 CCX at G.
Toy instances keep every narrow tick, so the exhaustive gate checks cover the kept step; `model_check.py`'s pass rows cover the unkept one.

## Handed-over held bits and a shared bottom carry

A non-last B1 hides H's bits n and n + 1 on two fresh wires instead of on l0 and the next cell's g0, and the tick hands that pair to the next tick as its word's top two wires.
They are 0 unless the previous cell is B1, so the held tick needs no cswap to move bit n off S's wire 0 or bit n + 1 off g0, and l0 stays 0 on a B1, so the B1 clear of l0 at the tick's end goes.
That is 3 CCX per tick with a previous and a next cell, for +2 clean wires at the tick's peak and 2 wires held between ticks.
The first chain's bottom carry l0·(one ⊕ K0) sits on a clean wire and survives its fall; after k = K + e, `cx(l0, c1)` and `ccx(l0, e, c1)` make it the second chain's l0·k0, since l0 = 1 implies one = 1.
That is 1 CCX per tick for +1 clean wire during the first chain.
Together they take G from 3,206,027 to 3,202,969 CCX at 3,135 wires.

## Square wire trades

At secp256k1 the square, not the multiplier, sets the peak: M peaks at 1,312 wires, and the square's junk sits beside the 516 live wires.
Four changes in `qec/square.py` spend the square's room up to the 4,096-wire cap.

- **Lent product ancillas.** E, T, t, c, k and the increment carries are first written after the Karatsuba product, so its ancillas take them first. The root's add onto its upper half pads V with 127 zero wires. The peak falls from 3,135 to 3,048 at the same gates.
- **Kept increments.** The fold's two increments above its low window keep their carries on n - m - 1 and n - b - 1 junk wires, at 1 CCX a bit in place of a Prefix's 3. The square's reverse clears them like any other junk fixed by Y. That is 412 wires for 1,648 CCX, 4 CCX a wire.
- **Zero product tops.** R(x) <= 2^(2ws - 2) and T(x) < 2^(2ws - 2), and nothing writes above V after the S child, so the top one or two wires of each S child's output hold 0 from then to the square's end. The fold and the gen-controlled modadd borrow those 14 wires, and the modadd's fresh scratch falls from 41 to 27.
- **Kept leaf carries.** KEEP = 650 junk wires hold the leaves' inner row carries, 2 CCX a wire across the square and its reverse, capped at the carries the leaves have. One KEEP wire is one peak wire, so KEEP is the knob that takes up any slack below the cap.

At G: 3,206,027 -> 3,203,079 CCX (-2,948), 3,135 -> 4,096 wires.
`scripts/square_check.py` checks X -= gen*Y^2 mod q, with gen and Y kept and every other wire zero, through the square, the modadd and the reverse. It runs every (gen, Y, X) at n = 4 to 8, every Y at 10 to 14, and samples at 24 to 48 and at secp256k1, over KEEP values from none to all.

## A partly kept boundary tick

`kept_plan` used to stop at the first tick whose whole chain did not fit in `KEEP`; that tick now keeps its lowest `KEEP - total` carries and runs the rest of its chain in place.
In `_kept_difference` with m kept wires, bits 1 to m build c_(i+1) on kept[i - 1] as before, and bits m + 1 and up run `_shared_difference`'s in-place majority, the first of them taking c_(m+1) from kept[m - 1].
That majority leaves c_(m+1) XOR a_(m+1) on kept[m - 1] during the chain, and its descent's CX restores c_(m+1) before the kept bits' descent, which reads only carry[i] for i <= m.
Each bit's descent produces the same U, V and g1 values as both chains do, so the step's outputs are unchanged, and the kept wires end holding c_2 to c_(m+1) of a + NOT(vhat), a function of the step's input that the inverse pass clears.
`KEEP` drops from 2400 to 2399, the most the multiplier holds without passing the 4,096-wire cap (KEEP = 2400 builds 4,097 wires); at G that is 2,399 kept carries, 53 of them on tick 314 (width 71).
The circuit loses 212 CCX, 4 for each of those 53 carries, at 3,187,369 CCX and 4,096 wires.
The required toy checks never build a partial tick, because toy plans fit in `KEEP`, so `scripts/keep_check.py` sets `KEEP` to every value from 0 to each prime's full plan (every odd prime below 600; 1,881 builds, 1,231 with a partial tick) and gate-checks every x with six y values, forward and inverse, at zero bad lanes.
With the restoring CCX of bit m + 1 removed, that run reports 26,268 bad lanes.

## The B1 sweep inside the Bezout ticks

On the merged path with a raw tape, the Euclid pass no longer runs the B1 sweep; the Bezout pass sweeps the tape as it reads it, and `_body` unsweeps it after the pass, before the inverse Euclid pass, with the old sweep's gates reversed.
This supersedes the placement in "Current mux contract": Euclid hands Bezout the unswept tape, and the map Bezout applies to it is unchanged.

Sweep step k swaps b_k and b_(k-1) when g0_k = g0_(k-1) = 0, for k from T - 2 down to 1, so it runs in the Bezout pass's order.
By S1, no B1 has a no-subtract cell above it in Euclid order, so when both g0 are 0, b_(k-1) is 1 when step k runs, and step k toggles b_k and b_(k-1) by t_k = NOT g0_k AND NOT g0_(k-1) AND NOT x_k, where x_k = b_k XOR t_(k+1) is b_k after step k + 1.
If t_(k+1) = 1 then x_k = 0, and if NOT g0_(k+1) AND NOT g0_k then b_k = 1, so t_(k+1) = NOT g0_(k+1) AND B1_k with B1_k = NOT g0_k AND NOT x_k.
So t_k = t_(k+1) XOR B1_k AND (g0_(k-1) XOR g0_(k+1)).

A wire tau carries t_(k+1) into tick k (`TwoBitMod._sweep_in`).
The tick XORs tau into b, making x_k; `_begin`'s CCX makes sw = g0 AND x_k; sw XOR g0 XOR b XOR 1 is B1_k, so one CCX with g0_(k+1) XOR g0_(k-1) on the next cell's g0 wire moves tau to t_k; then tau XORs into b again, which leaves the swept bit x_k XOR t_k there.
sw is still g0 AND b, because t_k implies NOT g0_k.
The rest of the tick reads the swept cell as before, and step k's toggle of b_(k-1) is the next tick's first XOR.
The pass's first tick has no previous cell and tau = 0 there; it takes g0_(k+1) as 1, so the CCX gives t_k by its definition.
Cell 0 has g0 implies b, so B1_0 = NOT x_0, and one CCX, tau ^= NOT g0_1 AND NOT x_0, clears tau, since t_1 = NOT g0_1 AND B1_0 and t_0 = 0.
The tape the pass leaves is the swept one, cell by cell, so `_body` applies the sweep's inverse, which restores the tape the inverse Euclid pass expects.

Per M application this removes the Euclid pass's sweep, 765 CCX, and adds one CCX to each of the 381 tick calls, one at cell 0, and one wire live through the Bezout pass, which takes `KEEP` from 2399 to 2398.
At G: 3,187,369 -> 3,186,609 CCX (-760), 22,311,583 -> 22,306,263 T, 4,096 wires.
Removing tau's CCX in `_sweep_in`, or ignoring tau in the tick, gives 41,988 bad lanes over the merged odd primes below 600 in `keep_check.py`'s multiplier check; both builds as committed give zero.
The Euclid pass keeps a `sweep` flag, on by default, so `model_check.py` and `utape_check.py` still check the swept pass.

## Swap bits on the tape

Each cell i from 1 to T - 2 now has a fourth wire, sw, zero when M starts, that holds the cell's swap bit s = g0 AND b from the Euclid step to the inverse pass.
The step writes it by CX at the point where held holds s, just before held takes the old offset, so the inverse pass, the same block reversed, clears it with that CX.
Cell T - 1 needs none, because its cell is (1, 0) or (0, 1) and s = 0 there; cell 0 needs none, because tick 0 always swaps, so s = g0 there, which the cell-0 tick already uses.

In the Bezout pass, tick k reads sw_k where it used to compute sw = g0 AND x_k onto a fresh wire with a CCX.
The two are equal: x_k = b_k XOR t_(k+1), and t_(k+1) implies NOT g0_k, so g0 AND x_k = g0 AND b_k = s_k; the sweep never changes a subtract cell's bit.
The tick then forms [cell k is B1] = 1 XOR g0 XOR b XOR sw on sw_k with the same X and CX gates as before, now in place, and tick k - 1 reads it as held.
Where tick k - 1 used to clear that mark with CCX(NOT g0_k, NOT b_k) and free the wire, it now applies X, CX from g0_k and CX from b_k, which return sw_k to s_k.
Tick k - 1 writes neither g0_k nor b_k, and b_k is the swept bit both times, so the mark is cleared exactly; tick 0 restores cell 1's.
So when the pass ends every sw_k is s_k again and the tape is swept, as before.

The unsweep after the pass is now the inverse of `_b1_sweep_sw`, which applies to (G, tau) the gates the ticks apply to them: b_k ^= tau, `_sweep_in` with sw_k for B1_k, and cell 0's clearing CCX.
The ticks' other gates leave every g0, b and sw wire as they found it, so the tape map the pass applies is this sequence's, and its inverse restores the raw tape for the inverse Euclid pass; tau is zero before and after both.
With s on a wire, B1_k = 1 XOR g0 XOR x XOR s is linear, so this sweep costs one CCX a cell, T - 1 in all, where `_b1_sweep`'s cswap and mark transition cost two.

Per M application this removes 382 sw CCX and 382 mark-clearing CCX from the ticks and 382 CCX from the unsweep, 1,146 in all, or 2,292 per point addition.
It adds T - 2 = 382 tape wires at the Bezout peak and removes the tick's own sw wire and the handed-over mark wire, so the peak grows by 380 and `KEEP` falls from 2398 to 2018, giving back 1,520 CCX.
At G: 3,186,609 -> 3,185,837 CCX (-772), 22,306,263 -> 22,300,859 T, 4,096 wires.
Dropping the X of the restore gives 15,084 failing lanes in `check.py`'s 36 toy instances; the build as committed gives zero there, in `keep_check.py` and in `model_check.py`'s gate rows, including the 832 secp256k1 width-edge lanes.
`DialogMultiplier.SWAP_WIRES = False` restores the previous construction; the fallback, coded-tape and short-pass paths never use it.

## Depth changes

The changes below cut scheduled ticks under the board's tick model, not CCX: CCX, T and the 4,096-wire peak do not move.
Several are ideas from iwnlcern's PR #21 (`v2/junk-depth`), re-derived and reimplemented on this tree's Euclid step and Bézout tick, which differ from PR #14's.
`scripts/step_check.py` runs every Euclid step variant at widths 2 to 7 against the step at `3cb492c` on every state of (U, V, E), restricted for a top mode to the domain `_top_bit` names, and needs the same outputs and clean scratch and copies; a dropped CCX must fail it.

### Parallel strip rotation

The strip's serial chain of w - 2 (or w - 3 under top_clear) controlled swaps V[i] <-> V[i + 1] is a left rotation of V[:L], L = w - 1 (w - 2), under held.
It is now the same permutation as two layers of disjoint swaps, all of V[:L] reversed and then all of it but the last wire, still L - 1 controlled swaps.
held is first copied by a doubling CX tree onto up to L // 2 - 1 clean wires, which control the swaps round-robin, and the reversed tree clears them.
The copies are the held wires of the pass's last cells, which no step writes before its own, so they are zero when step i runs, and the inverse pass, undoing step i only after every later step, finds them zero too.
The step takes them as one more register, so no wire is allocated.
At G: 35,620,290 -> 32,524,514 ticks, 3,185,837 CCX, 4,096 wires.

### Euclid selects outside the carry chain

In `_shared_difference` each in-place bit i >= 1 of the chain used to select inside the unwind: restore a_i, then select S_i and apply the swap's correction to U[i], all on the carry chain.
With S = carries[i] = c ^ a, B = V[i] = b ^ a and A = U[i] at the top of the chain, the selects now go in first, top bit down: S ^= g1*B, then B ^= g0*NOT(S).
The unwind is then A ^= S*B and S ^= A alone, one CCX and one CX a bit: A becomes a ^ g1*(b ^ a), since S*NOT(S) = 0 leaves the restore's product S*B as it was apart from g1*B, and S becomes c.
Three CX a bit after the chain, B ^= A, A ^= g1 and B ^= g0, with the existing B ^= g1, give the old outputs.
The pre-select needs g1 => g0, which holds wherever g1 is made: CCX(g0, ...) with no top mode, and H = g0*MAJ(a, b, c) from `_top_bit` on its mode's domain.
Bit 0 loses its empty X pair: c ^= a; v ^= a; v ^= c is v ^= c; c ^= a, then v ^= c.
The CCX are the same three a bit. `_kept_difference`, bit 0's carry and the top bit are unchanged.
This is the identity of PR #21's commits 3 and 6, derived again here for this step, whose in-place branch is the same; `scripts/step_check.py` checks it on every state, within each top mode's domain.
At G: 32,524,514 -> 26,647,800 ticks.

### Prefix bits incremented before the uncompute

`Prefix.write` used to clear each block bit's AND-chain wire and then add the block's carry p into the bit, alternating, so p and the chain shared one serial path.
Each bit now takes x ^= p*NOT(below) before its chain wire is cleared, and x ^= p after the block's chain is down.
The uncompute CCX(below, x, u) reads below*x, and below*(x ^ p*NOT(below)) = below*x, so it is unchanged; the two increments sum to p*below, the old one.
No later gate of the block reads x[i] or writes p, so the deferred CX can wait; the first two bits take the same form with t.
This is PR #21's commit 5; our `Prefix.write` is the one it rewrote, and it serves the tick's P and M and the square's carry.
`scripts/prefix_check.py` runs Prefix and write at every length to 14 bits, with and without release and restore, against `3cb492c` on every input and checks bits + t, with dropped-CCX controls; `square_check.py` passes.
At G: 26,647,800 -> 24,933,111 ticks.

### Copies of the coefficient swap's control

The tick's last step, the controlled swap of R and S under sw, had all 256 CCX on sw, one serial path.
sw is now copied by a star of CX onto free wires, each copy controls an equal share of the swaps, and the reversed star clears the copies.
A wire on the free stack is zero, as every wire is when it is freed, so each copy equals sw and the swaps are the same.
Only wires the tick has touched are used, so the tick's arity and peak do not move; wires touched since the compare began are left out, since their copies would wait on the compare, and in the inverse tick hold its start.
This is PR #21's commits 2 and 7.
At G: 24,933,111 -> 23,458,995 ticks.

### A depth pass at emission

`qec/depth_pass.py` is upstream's (`origin/main`), unchanged. `write_circuit` runs it on the hierarchical branch by default: per definition it reorders commuting ops by list scheduling under the board's CCX lowering, keeps the input unless a schedule is shorter, and then orients each CCX's controls; the top level only gets canonical control order.
Its decisions read definition bodies alone, so definitions stay identical over addends.
The flat branch is never changed, so `check_hier.py` emits with `depth_pass=False` against the flat pin, as do `keep_check.py` and `model_check.py`, which test the builder's M; `check.py --hier` and accept's emissions simulate the passed circuit.
`scripts/pass_check.py`, which shares no code with the pass, certifies each definition: the same multiset of ops (a CCX's controls as a set), and on every wire the same runs of compatible roles in the same order, so every pair of non-commuting ops keeps its order; an exchanged conflicting pair must be rejected.
At G: 23,458,995 -> 19,796,016 ticks.
