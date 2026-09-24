# Instantiate

`instantiate` generates solver `data` from symbolic template inputs and the
selected Wilson-chain backend. `--diag-seed-only` stops after generating and
diagonalizing the seed; `--wilson-only` generates chain files without a seed.

## Wilson Coefficient Scaling

Full instantiation and `--diag-seed-only` require `rescalexi=false` in `[param]`
(the default). Seed Hamiltonians and solver `data` require unrescaled Wilson
coefficients; the rescaled hopping convention is not supported in these modes.
With `rescalexi=true`, both modes fail before configuration reporting, chain
generation, seed diagonalization, or creating/publishing any output files.
Existing `data`, parameter-section sidecars, and temporary artifacts are left
unchanged, including when `--generate-temporaries` is requested.

`instantiate --wilson-only` still accepts `rescalexi=true` for rescaled chain
output. This restriction is specific to seed-producing `instantiate` modes;
standalone `nrgchain` behavior, output names, and coefficient conventions are
unchanged for both `tridiag_method=lanczos` and `tridiag_method=rkpw`.

The `instantiate1_u1_perchannel_legacy` and `instantiate1_u1_perchannel_rkpw`
regressions cover rejection without modifying caller-owned artifacts,
default/false acceptance in both seed modes, and rescaled Wilson-only output.
After the affected target and test inputs have been rebuilt/staged, run:

```sh
ctest --test-dir build -R '^instantiate1_.*_(legacy|rkpw)$' --output-on-failure --no-tests=error
```
