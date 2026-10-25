# Validation guide (for an independent validator)

The claim (`THEOREM.md`) rests on three things:
1. an elementary lemma, with one textbook ingredient (the Gram–Schmidt lower bound for the dual lattice);
2. the zero data;
3. a finite, exact or ball-arithmetic computation.

## 1. Integrity
- Download the bases from the Zenodo record into `bases/`, then run `sha256sum -c SHA256SUMS.txt`.

## 2. The argument
- Check the lemma and its proof in `THEOREM.md`.
- Compare with Best–Trudgian 2015, Theorem 3 and Lemma 3, and Lagarias–Lenstra–Schnorr 1990.
- Points to confirm:
  - the lattice L = {(v, z) : v_j ≡ z Y_j (mod 2^ν)} has the dual vector u = w/2^ν for every relation;
  - the bound ‖w‖² ≤ M²(N + E1²) holds for |m_j| ≤ M;
  - the statement excludes every nontrivial relation in the box, with m₀ free.

## 3. The zeros
- `python verify/verify_li.py zeros 0 1` re-proves |γ_i − mid_i| < 2·10⁻⁹⁶⁶⁵ for i ≤ 700. Measured cost: a few seconds per zero.
  - It can be split: `zeros k K` handles the indices i ≡ k + 1 (mod K).
  - When all 700 proofs exist in `verify/output/`, `cert` and `table` use them instead of the shipped records.
- Optional, strongest: compare the first digits of the midpoints with any published table of zeros, e.g. A. Odlyzko's tables or the LMFDB.

## 4. The certificates
- `python verify/verify_li.py table` checks all six bases against `results/LI_TABLE.json`. Each output line reports:
  - membership (exact);
  - full rank (ball LDL plus rank mod 2⁵⁹ − 55);
  - the largest certified M;
  - whether M is at least the table value.
- Optional, strongest: re-implement the check yourself, e.g. with an exact rational LDL or with exact integer determinants of the leading Gram minors (slow for n = 700, fine for n ≤ 200).
  - The basis format: one vector (v_1 … v_N, z) per line, signed hexadecimal.

## 5. What is deliberately not included
The procedure that produced the reduced bases. It is not needed: any basis of L that passes the checks certifies the stated M.
