# Posterior cones and sharp reconstruction thresholds on finite-mean Galton–Watson trees

Yitzchak Shmalo, October 2026

Paper: [paper/posterior_cones_reconstruction.pdf](paper/posterior_cones_reconstruction.pdf) (16 pages). LaTeX source in `paper/`.

For symmetric ferromagnetic Potts broadcasting with three or four states on an observed Galton–Watson tree with any finite mean offspring d, reconstruction holds exactly when dλ² > 1, for 0 < λ < 1. At equality the root cannot be reconstructed. No assumption on the offspring variance or tail is needed.

The three-state argument uses an explicit invariant cone of posterior laws built from H = y − x²/2, where x and y are the quadratic and cubic centred posterior invariants. Its two global inequalities are certified by exact nonnegative polynomial coefficients. The four-state case combines the same argument with a posterior cone from earlier work, which the paper credits.

## Contents

- `paper/`: the manuscript and its sections.
- `analytic/q3_exact_pair_simple.cpp`: the first three-state coefficient program.
- `computation/q3_independent_coefficients.cpp`: an independent integer implementation. Both produce the same 80,532 coefficients, given in full in `computation/q3_independent_coefficients_A.tsv` and `_B.tsv`.
- `computation/published_q4/`: the four-state verifier, its licence and the output of two complete traversals.
- `SUPPLEMENT.md`: how to compile the paper and reproduce every certificate.
- `MANIFEST.json`: SHA-256 checksums and byte lengths for every file.

All sign decisions use exact integer or rational arithmetic.
