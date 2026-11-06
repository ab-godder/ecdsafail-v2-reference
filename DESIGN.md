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

M itself is total. Euclid runs on (u, V) with V = 2v + u, so each step is one uncontrolled add between two b0-controlled complements. For x = 0 a flag z = [x == 0] loads (u, V) = (1, 1) instead of (q, q), so every step idles on registers that fit the shrinking widths and the dialog is all zero. z then drives the controlled swap that puts y back in place, and is uncomputed from the restored x.

## Data shape

- **Wire.** An integer index. **Register.** A Python list of wires, least significant bit first. Shifts are relabels, implemented as list rotations with no gates.
- **Circuit.** A flat `array('i')` of (opcode, a, b, c), a free-list wire allocator that requires ancillas to be returned clean, and a high-water mark that becomes `num_wires`.
- **Block.** A recorded op range with input layout, output layout and internal wires. It can be instantiated forward or reversed onto any matching layout. That is how M^-1, Euclid-backward and every compute/uncompute pair are emitted without hand-written inverse code.
- **Simulator.** A bit-sliced simulator in which each wire is a Python int and bit j is test lane j. It shares no code with the generator's arithmetic.

## Why the permutation is enough for QNT-04

The circuit contains only X, CX, CCX and SWAP. It is therefore a permutation matrix with no phases. If every promise basis state |c, R, 0> maps to |c, R + cA, 0>, the same matrix acts exactly on any superposition of promise states and on any extension with a reference system. Basis exactness plus clean scratch is the whole unitary contract. Outside the promise, the circuit is still a permutation, so QNT-06's valid-channel requirement holds by construction.
