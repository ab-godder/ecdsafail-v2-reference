# The mux merge on a held offset

The cell is in `qec/two_bit.py`, and the multiplier change is in `qec/point_add.py`.
Each statement is numbered so that it can be tested as worded.

## Credit

The representation T = S + R between a doubling and the following subtract is section VI of arXiv:2609.05625.
No novelty is claimed for the offset.
That paper's complement form leaves (c − 1)/p of inputs off the canonical range, so it is not used here.
The subtract below is a true integer map on the domain stated in D1 to D4.

## Notation

q is an odd prime of n bits and f = 2ⁿ − q, which is odd.
A prime takes the new path when `merged` holds: 4f < 2ⁿ and corner ≤ 4f, where corner = f ⊕ 2f ⊕ 3f.
Only 4f < 2ⁿ is used below; the second clause keeps the set of primes on the path where it was.
w = bitlen(4f).
Every other prime keeps today's `step`, with no sweep.

A cell is (g0, b), and the symbols are A0 = (1, 0), A1 = (1, 1), J = (0, 1), B1 = (0, 0).
A subtract cell is one with g0 = 1.
Ticks run in Bézout order, the reverse of the tape, so tick t reads cell T − 1 − t.

| Lane | Condition | Map |
|---|---|---|
| Plain add | g0 = 1, previous cell not B1 | S = 2T + R mod q, then swap on b |
| Held | g0 = 1, previous cell B1 | S = 4T + R mod q, then swap on b |
| J | J | S = 4T mod q |
| B1 | B1, not the last tick | no reduction, H = 2T + R is held |
| Last B1 | B1, last tick | S = 2T mod q |

## Domain

Every statement of exactness below names this domain.

**Top level.** R is on the curve or infinity, in the canonical encoding, with zero scratch.

**Tick level.** T and R are canonical, scratch is zero, and:

- **D1.** A cell after a B1 is a subtract cell.
- **D2.** A B1 that is not the last tick has next g0 = 1 on entry.
- **D3.** A B1 on the last tick has R = 0.
- **D4.** On a held tick, the input is the stored form of S4, not a canonical word.

Today's cell was exact on every symbol and canonical R, S.
This tick is not, so a per-cell all-input check does not apply to it.
D1 to D3 are properties of a swept Euclid tape, by S1 to S3.

## The tape

**S1. Run shape.**
In Bézout order, each maximal run of cells with g0 = 0 holds at most one B1, and it is the run's first cell.
`tick_tape` writes B1 exactly when v ≡ 2 mod 4, which leaves v odd, so the next tape cell is a subtract or the tape ends.
Reversing the tape gives the statement.

**S2. The sweep.**
`_b1_to_run_end` walks t = 0 to T − 2 and swaps b between cells t and t + 1 when both have g0 = 0.
It costs 3 CCX per adjacent pair and one clean wire.
It reads g0 and never writes it.
In each run the b values rotate left by one, so a B1 moves to the run's end and the Js shift down.

After the sweep every B1 is followed by a subtract cell or is the last cell, which gives D1 and D2 on the pass's own tape.
The number of ticks, the alphabet, the symbol counts, `tick_bound` and `tick_widths` are unchanged, and `model.py` is not edited.
B1 and J only double S, R is unchanged across a run, and doublings commute, so S and R after each run equal their unswept values.
The pass's result is therefore `tick_bezout`'s.
`_body` applies the sweep before the first tick and its reverse after the last, so the tape handed to the inverse Euclid block is the original.

**S3. The last tick.**
If the last cell is a B1 then R = 0 on it.
The pass ends with R = 0 by the existing proof of the recurrence, and a B1 does not change R.

## The integers

**S4. Stored form.**
Let H = 2T + R on a B1 that is not the last tick.
Then H < 3q < 2ⁿ⁺², and H ≡ R mod 2.
After the tick, in the rotated S register, wires 1 to n − 1 hold H's bits 1 to n − 1 and wire 0 holds 0. The tick hands the next tick a pair of wires, `hand`, that holds H's bits n and n + 1; on every other lane both are 0.
H's bit 0 is not stored, because it equals R₀, which stands in R.
The previous cell's wires are unchanged, and they still read B1.
So a held tick can recover H from its own inputs.

**S5. The held tick's word.**
Let W = H − R₀ − (R >> 1), which is 2T + (R >> 1).
Then 0 ≤ W < 2ⁿ⁺².
H − R₀ is H with bit 0 cleared, so the subtract needs no carry in.
Let V = 2W + R₀.
Then V = 4T + R, and 0 ≤ V < 5q.
Nothing on this tick is negative, which is the reason for holding.

**S6. The reduction.**
Let V ≥ 0, V < 5q, low = V mod 2ⁿ and K = V >> n, so K is at most 4.
Let probe be f, 2f, 3f, 4f, f for K = 0, 1, 2, 3, 4, and e = [low + probe ≥ 2ⁿ] and k = K + e.
Then k = ⌊V/q⌋, k ≤ 4, [low + kf ≥ 2ⁿ] = e, and S = low + kf − e2ⁿ = V − kq.

For K ≤ 3, V ≥ (K + 1)q exactly when low + (K + 1)f ≥ 2ⁿ.
For K ≤ 2, V < (K + 1)2ⁿ ≤ (K + 2)q because (K + 2)f ≤ 4f < 2ⁿ.
For K = 3, V < 5q bounds k by 4.
For K = 4, low < 2ⁿ − 5f, so low + 4f < 2ⁿ, e = 0 and k = 4.
For K ≤ 3, if e = 0 the write has no carry, since kf ≤ probe, and if e = 1 the carry is the probe's, with kf = probe.

K = 3 needs a true 4f, where the earlier low-chain cell, which is not in this repository, used corner as a don't-care.
A separate check, also not in this repository, shows that corner fails on the held lane at secp256k1.

**S7. The clear.**
Write k = k₀ + 2k₁ + 4k₂, on wires k0, k1, k2.
S₀ = V₀ ⊕ k₀ on every lane, since f is odd, so the cell does not write bit 0 and `cx(l0, k0)` turns k0 into S₀.

- **Plain add.** V = 2T + R < 3q, so k ≤ 2, and k₂ = 0. Odd k is k = 1. For even k, k = 2 exactly when R >> 1 > S >> 1, because S and R agree on bit 0 and k = 0 gives S ≥ R while k = 2 gives S = R − 2(q − T) < R.
- **Held.** S ≡ R + kf mod 4, because 4T ≡ 0 and q ≡ −f. So k₁ = S₁ ⊕ R₁ ⊕ k₀R₀ ⊕ (f₁ · k₀), where f₁ is bit 1 of f. If k₀ = k₁ = 0 then k is 0 or 4, and k = 4 exactly when R >> 1 > S >> 1, by the same argument with S = R − 4(q − T) < R. There is no offset of f.
- **J.** S ≡ kf mod 4, as today.

On a held lane `cswap(p, k1, k2)` first exchanges the two wires, so the compare's target carries k₂, and wire k2 carries k₁.
`mark_even` makes z = (g0 ⊕ wire k2) · ¬k0, which on a held lane is 1 exactly when k is 0 or 4.
On other lanes wire k2 is 0, so z is today's g0 · ¬k0.

## The add and the shift

**S8. Complementary controls.**
Every tick adds or shifts: off J it adds R >> 1 to the word, and on J it shifts the word left.
On a held lane the word is complemented by `cx(p, x)` over all n + 2 wires before and after, and ~(~x + a) = x − a mod 2ⁿ⁺², which is S5's W.
A B1 adds like a plain add, with l0 left at 0 by g0 = 0.

**S9. The rolling XOR.**
Let a_i be `ops[i]`, c_i the carry into bit i of `tgt`, and x_i = a_i ⊕ c_i, so an add makes word[i] ⊕ x_i for i ≤ n − 2.
After the up pass D = x_0 ⊕ … ⊕ x_(n−2), and it is built by CX only.
Going down, pair i takes word[i] and word[i + 1], and both receive one value: D on the add lanes, and word[i] ⊕ word[i + 1] on J.
One CCX does both, as `cx(v, u); cx(D, v); cx(u, D); ccx(J, D, v); cx(u, D); cx(v, u)`.
On J that pair is the ladder's swap, and pairs run top down as the ladder does.
Off J, word[i] receives D_i from pair i and D_(i−1) from pair i − 1, whose sum is x_i, and word[0] receives only D_0 = x_0.
Then `cx(c, D)` removes x_i, so D is 0 after pair 0 on every lane.
Each bit costs 1 CCX up, 1 to restore a_i and 1 for the pair, against 4 today.

**S10. The top of the word.**
The pairs are i = 0 to n − 2, so there are n − 1 of them, and the gross saving is 255, not 256.
Above them, four things happen in this order, on `J`, `carry` and `t`:

1. `carry` = ¬J · c_(n−1).
2. t = carry · word[n − 1], then word[n + 1] ⊕= t · word[n], then word[n] ⊕= t, then word[n − 1] ⊕= carry, then t is cleared from the new word[n − 1] under an X. That is the ripple through the two top bits, 2 CCX over today's half adder.
3. `carry` is cleared, then word[n − 1] ⊕= ¬J · D. Pair n − 2 gives word[n − 1] the value D_(n−2), and no pair above cancels it, so this correction is needed once. 1 CCX.
4. `cswap(J, word[n − 1], word[n])`, J's top pair, 1 CCX. It stays a shift because the top wires hold a carry on add lanes.

Folding the correction into the top pair only moves the junk to `spill`.
On J, c_(n−1) is ignored, D is unused until pair n − 2, and the top pair runs first, as the ladder's does.

## The lanes

**S11. B1 costs no CCX to read inside a bracket.**
Between `_begin` and `_end`, b holds J and g0 is true, so g0 ⊕ b is 0 on a B1 and 1 on the other three symbols.
The first bracket stays open through the first chain, so `one` = g0 ⊕ J is two CXs, up and down.

**S12. The previous cell.**
[previous cell is B1] needs a real AND, so `_after_b1` costs 1 CCX and `_drop_after_b1` 1 more, read twice per tick, 4 CCX.
Its wire does not stay through the reduction, because the reduction is the peak.

**S13. B1 is inert in the reduction.**
On a B1 that is not last, before either chain, k0 and k1 move to two fresh wires, the pair the tick hands over.
The next cell's g0 is not touched, so it stays 1 by D2, and l0 = 0 there by g0 = 0 for the whole tick.
Then every selector is 0, so every logical carry of both chains is 0, and neither chain reads L[0] into anything.
The first chain's carries sit on H's wires, which equal the logical carries only when P = 1.
When P = 0 they may hold junk, which nothing reads except through P and which the first chain's reverse removes.
So e = 0, and every carry of the second chain is 0, both on its clean wires and on the selectors it perches on, which are 0 here.
`P.write(0)` and `M.write(0)` restore their words.
At the end `cx(l0, k0)` leaves wire k0 at 0, and `ccx(g0, R0, l0)` leaves l0 at 0.
That is S4's form.
The first chain's bottom carry l0·(one ⊕ K0) is on a clean wire and is not cleared at its fall; `cx(l0, c1)` and `ccx(l0, e, c1)` turn it into the second chain's l0·k0 = l0·(K0 ⊕ e), because l0 = 1 implies g0 = 1 and so one = 1. On a B1 both are 0.

**S14. Entering a held tick.**
The handed-over pair is the tick's spill and high, word[n] and word[n + 1], so on a held lane word holds H with word[0] = 0, which is H − R₀ by S5, and on every other lane they are 0 as before.
g0 was never borrowed, so the tape is whole.

**S15. The last tick.**
With nxt = None, `one` = 1 and nothing is hidden.
By S3 a B1 there has R = 0, so its add is 0 and V = 2T < 2q, and the reduction runs as on any other lane.
z = 0 on it, and k₁ = 0, so the clear is S₀ alone.

**S16. The first tick.**
With prev = None there is no held lane.
Two gates are then inert on the whole domain: `ccx(t, word[n], word[n + 1])`, since word[n] = 0 on every lane, and `ccx(three, e, k2)`, since K = 3 with e = 1 needs V ≥ 4q, which a tick with no held lane cannot reach.
Those gates remain in the general tick, but M no longer calls that tick first.
M starts with S = 0 and its terminal tape cell is A0, A1 or J, never B1:
the final active cell is a subtract, and every padded cell after it is J;
the offset sweep changes b only on subtracts and the B1 sweep cannot introduce
B1 into that all-J tail. On A0 and A1 the first tick gives (R, R), and on J it
gives (R, 0). Thus M uses S XOR= g0*R, one CCX per bit, with R, cell and next
unchanged and no scratch. S keeps its wire list for the next tick. This uses
the existing tick bound, including the proved 384 bound at secp256k1, and reads
only the terminal cell's g0, so the tail b lend is unaffected. The inverse
uses the same XOR on the image of this boundary. Fallback steps are unchanged.

## The whole circuit

**S17. Both directions.**
`_body` is one recorded block, so M⁻¹ is its reverse, gate for gate.
Forward maps the domain (x, y) with canonical y and clean scratch to (x, x′y mod q), which is every canonical y again because x′ is invertible, so the reverse runs on exactly the forward image.
`_kernel` applies M⁻¹ first and then M.
With gen = 0 the two cancel at circuit level, as before.
The reverse of the pass applies the sweep, then the ticks backward from the last, then the sweep's reverse, so a reversed tick reads the swept tape and the stored form of S4.

**S18. Other moduli.**
The generator runs at every q in the gate.
Primes in the gate on the new path: 7, 13, 31, 61, 97, 127, 251, 397, 433, 509, 823, 907, 967, 1021, 4093, secp256k1.
Primes in the gate on today's path: 3, 5, 11, 43, 67, 79, 163, 211, 349, 547, 577, 733.
`f₁ = 1`, the `f >> 1 & 1` branches, is reached only at toy primes, because secp256k1 has f ≡ 1 mod 4.
Empty H and one-bit H are reached only at toy primes.
Both prefixes at their real lengths, H = 221 bits and the quiet run of 20, are reached only at secp256k1, and `Prefix` is byte-identical to the earlier low-chain cell's.

**S19. Peak.**
The cell needs 26 clean wires, the same as the earlier low-chain cell, at the top of the second chain.
`high` costs one more there, and clearing c₁ once c₂ exists pays it back at +2 CCX, by the same carry rebuilt on the way down.
c₁₂ still sits on `one`, because a perch can now read through other selectors, here k2.
The two held bits need no clean wire: S4 keeps them on wires the tick already owns.
D is live only during the add, which is far below the peak.

**S20. Limits of this proof.**
It proves maps, and the wires by the arguments above and by test, not by a formal method.
`lift` on the new cell has not been tried, and its result belongs in the limit list either way.
Nothing here passes the checker, and `is_qualified` stays false.
