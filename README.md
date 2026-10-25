# Certified m-independence of the first n zeros of the Riemann zeta function (n ≤ 700)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23652325.svg)](https://doi.org/10.5281/zenodo.23652325)

Let γ₁ < γ₂ < … be the ordinates of the nontrivial zeros of ζ(s) in the upper half-plane. This repository proves, for the (n, M) below, that **there is no integer relation**

> m₀γ₁ + m₁γ₂ + … + m_{n−1}γ_n = 0   with   (m₁, …, m_{n−1}) ≠ 0   and   max_{j≥1} |m_j| ≤ M,

with m₀ unrestricted. In particular the first n zeros are M-independent: every nontrivial combination of them with integer coefficients of absolute value at most M is nonzero.

| n | log₁₀ M (this work) | Best–Trudgian 2015, Table 2 |
|---|---|---|
| 100 | 94.42 | 86.26 |
| 200 | 45.41 | 39.06 |
| 300 | 28.64 | 21.79 |
| 400 | 19.97 | 12.37 |
| 500 | 12.89 | 5.74 |
| 700 | 5.41 (M = 257 742) | not given |

The exact values of M are in [THEOREM.md](THEOREM.md) and [results/LI_TABLE.json](results/LI_TABLE.json). Best and Trudgian's table (Math. Comp. 2015) goes up to n = 500. The bounds here improve it by 6 to 8 orders of magnitude, and extend it to n = 700. Our statement is also slightly stronger, because m₀ is not bounded.

## How the proof works

For each n, a reduced basis of a lattice built from the ratios γ_{j+1}/γ₁ certifies, through a lower bound for the dual lattice, that no short relation exists. This is the method of Best and Trudgian, with a sharper bound. The lemma and its proof are in [THEOREM.md](THEOREM.md). The check is a finite computation:
- the zeros are given as decimal midpoints, and their accuracy (2·10⁻⁹⁶⁶⁵) and indexing are proved with Hardy's Z-function and Turing's method;
- every basis vector is checked to lie in the lattice, exactly;
- full rank and an upper bound for every Gram–Schmidt length are certified in ball arithmetic;
- M is computed in exact rational arithmetic.

The certificate does not depend on how the bases were reduced.

## Check it yourself

The bases (4.5 GB in total, gzip text) are in the Zenodo record [doi:10.5281/zenodo.23652325](https://doi.org/10.5281/zenodo.23652325). Put them in `bases/`. Then:

```bash
pip install python-flint                 # FLINT/Arb ball arithmetic
cd verify
python verify_li.py table                # all six bases; "all claims verified" at the end
python verify_li.py zeros 0 1            # optional: re-prove the 700 zero enclosures (otherwise the shipped proofs are re-checked by hash)
```

Measured on a desktop computer (AMD Ryzen 7 9700X, one thread per job): the whole table took 6.6 minutes (n = 700: 200 s), with peak memory 4.8 GB. Re-proving the 700 zeros took 1448 s as two parallel jobs (`zeros 0 2` and `zeros 1 2`, about 4 s per zero per job); the result was identical to the shipped proofs. A single basis: `python verify_li.py cert 700 ../bases/LI_n700_nu32000.txt.gz 257742`.

## Contents

```
README.md, THEOREM.md, VALIDATION_GUIDE.md
data/zeros_1_701_midpoints.txt   decimal midpoints of γ₁ … γ₇₀₁ (65 536-bit computation), and data/README.md
results/LI_TABLE.json            the certified (n, M), with the file name and SHA-256 of each basis
results/zero_proofs_shipped.jsonl  the Hardy-Z proofs of |γ_i − mid_i| < 2·10⁻⁹⁶⁶⁵, i ≤ 700
results/certificates.jsonl       the output of verify_li.py table on the published bases
verify/verify_li.py              the verifier
bases/                           (Zenodo only) LI_n<n>_nu32000.txt.gz, one lattice vector per line in signed hexadecimal
SHA256SUMS.txt                   checksums of every file, including the bases
```

## What is not included

The procedure that reduced the bases is not part of this repository. It is not needed to verify the result.

## How this was made

This work is directed by Benjamin Weisz, who sets its goals, decides on publication and is responsible for its content. Nearly all of the code, computations and documentation were produced by an AI system (Claude, by Anthropic) working under that direction. The result does not rest on trusting the AI:
- the certificate is checked by exact integer arithmetic and ball arithmetic;
- the table values were first computed by one program and then confirmed by the separately written verifier published here, which uses a different arithmetic path;
- the zero data are proved by the verifier itself (Hardy-Z sign changes with Turing indexing).

**Review status.** This work has not been peer reviewed. Its proofs have been checked only by the verification programs and by AI agents; no human has reviewed them, and the owner does not count himself as a reviewer.

## Related work

The companion result on the oscillation of the Mertens function, lim sup M(x)/√x > 2.06175 and lim inf M(x)/√x < −2.06471, is at https://github.com/wbeni95/mertens-function-bounds (doi:10.5281/zenodo.23298194).

## Licence and citation

Code (verify/) under the Apache License 2.0; data and documentation under CC BY 4.0: see [NOTICE](NOTICE). Please cite the Zenodo record [doi:10.5281/zenodo.23652325](https://doi.org/10.5281/zenodo.23652325); the metadata is also in [CITATION.cff](CITATION.cff).
