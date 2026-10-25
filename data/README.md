# Data: `zeros_1_701_midpoints.txt`

Plain text: two `#` header lines, then 701 lines; line i holds the decimal midpoint of γ_i, the ordinate of the i-th nontrivial zero of ζ(s) in the upper half-plane (19 724 significant digits each).

- **Computation.** FLINT/Arb, 65 536-bit precision, from rigorously isolated zeros.
- **What the verifier uses.** Only the midpoints. It does not use any stored radius: for i ≤ 700, `verify/verify_li.py zeros` proves |γ_i − mid_i| < 2·10⁻⁹⁶⁶⁵ directly.
  - It shows a sign change of Hardy's Z-function at mid_i ∓ 1.999·10⁻⁹⁶⁶⁵, at 32 600 bits.
  - It checks that the window [mid_i ∓ 2·10⁻⁹⁶⁶⁵] lies strictly between the Turing-indexed isolating balls of zeros i − 1 and i + 1 from `acb.zeta_zeros`, and meets the ball of zero i.
  - So the zero in the window is γ_i, the i-th zero.
- **Shipped proofs.** The proof records are in `results/zero_proofs_shipped.jsonl`. Each record carries the SHA-1 of the midpoint string it proves.
- **Line 701.** It is included for completeness; the verifier uses the isolating ball of zero 701 from Arb, not this midpoint.
