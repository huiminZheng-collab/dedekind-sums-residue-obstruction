# A unipotent residue obstruction for higher-weight generalized Dedekind sums

Preprint and exact-arithmetic reproducibility files, version **1.0.0**.

[Read the paper](paper/obstruction.pdf) · [LaTeX source](paper/obstruction.tex) · [Reproduce the calculations](REPRODUCIBILITY.md)

The paper proves a common coefficient-congruence criterion for period-polynomial cocycles. It shows that four higher-weight generalized Dedekind-sum images are proper subsets of their additive spans, disproving the group-image assertion in the second bullet of Tranbarger's Conjecture 1.15.

| Characters and weight | Group | Additive span |
|---|---|---|
| chi_3, chi_3, weight 4 | Gamma_1(9) | 2 Z |
| chi_3, chi_4, weight 4 | Gamma_1(12) | (2/3) Z |
| chi_3, chi_5, weight 7 | Gamma_1(15) | (4/3) Z |
| chi_3, chi_3, weight 6 | Gamma_1(9) | (10/3) Z |

A separate square-residue law modulo 9 strengthens the second example. The printed weight-7 containment constant 4 is contradicted by a value already in the source's finite sampling set. The other three containment constants in this table remain valid. The paper does not classify the exact images within the surviving residue classes.

Problem source: P. Tranbarger, *Higher Weight Generalized Dedekind Sums*, [arXiv:2512.17139v1](https://arxiv.org/abs/2512.17139v1), Conjecture 1.15. The source is cited for its definitions and conjecture, not as a source of the counterexample proofs.

## Contents and verification

- `paper/`: the frozen 14-page manuscript and its source.
- `src/`, `certificates/`: standard-library Python programs and finite exact inputs/results.
- `AI-USAGE.md`: contribution disclosure.
- `RELEASE-MANIFEST.sha256` and its detached OpenPGP signature: byte-level payload identification.
- `provenance/`: the release public key and DigiCert certificate chain.
- `timestamp/`: RFC 3161 query, response and verification record over the manifest signature.
- `artifacts/reproducibility-v1.0.0.zip`: a self-contained copy of the release, excluding this archive itself.

Run `python verify_release.py --gpg gpg --openssl openssl` to verify the manifest, signature and timestamp. Add `--replay` to execute the mathematical certificates in a temporary copy. See [reproducibility instructions](REPRODUCIBILITY.md) for requirements.

Release-key fingerprint: `490689A7D9DB559372C8620192664F521ACFE6E2`.

The manuscript is content-anonymous. The hosting account and public commit metadata may identify the depositor; author-unlinkable double blindness is not claimed. This is a public preprint release, not a claim of journal acceptance or peer review. The timestamp fixes existence of exact bytes, not mathematical correctness or novelty.

Keywords: generalized Dedekind sums; higher-weight Dedekind sums; Tranbarger Conjecture 1.15; counterexamples; Eisenstein series; period polynomials; cocycles; quadratic characters; unipotent action; Bernoulli polynomials; exact arithmetic.
