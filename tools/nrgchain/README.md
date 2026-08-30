# `nrgchain`

`nrgchain` determines the Wilson-chain hopping coefficients `xi`, on-site
energies `zeta`, and hybridization weight `theta`. For a tabulated band it uses
`FSOL.dat` and `FSOLNEG.dat` produced by `adapt`; adaptive meshes additionally
require `GSOL.dat` and `GSOLNEG.dat`.

## Usage

```text
nrgchain [options] [s|l] [parameter_file]
```

- `parameter_file` defaults to `param`.
- `s` calculates and saves the intermediate `de_*` and `du_*` tables without tridiagonalizing.
- `l` loads those tables and `theta.dat`, validates them, and tridiagonalizes without reading the density or `adapt` output.
- `-v` writes the resolved configuration and diagnostics to standard error.
- `-vv` increases verbosity further.
- `-V` or `--version` prints the project version.
- `-h` or `--help` prints the command synopsis.

A command-line table mode overrides `nrgchain_tables_save`,
`nrgchain_tables_load`, and `nrgchain_tridiag` from the parameter file. The
legacy misspelling `nrgchains_tridiag` remains a fallback when
`nrgchain_tridiag` is absent. Without a command-line override,
`nrgchain_tables_save=true` and `nrgchain_tables_load=true` are mutually
exclusive. `instantiate` requires tridiagonalization to be enabled.

## Density interpolation

The parameter

```ini
density_interpolation=linear
```

selects the representation of a tabulated hybridization density. Supported
values are:

- `linear`: compatibility default using piecewise-linear interpolation.
- `steffen`: shape-preserving piecewise-cubic interpolation. It stays within each interval's endpoint range and therefore preserves nonnegativity of nonnegative samples.

`adapt P`, `adapt N`, and `nrgchain` must use the same setting. Changing the
method requires regenerating all `GSOL*` and `FSOL*` files before running
`nrgchain`; these legacy files do not contain interpolation metadata.

Point values and shell weights are derived from the same interpolant. Shell
weights use its analytic interval primitive rather than interpolating a
separately integrated table. Density values must be finite and nonnegative,
energies must be finite and strictly increasing after branch selection, and
Steffen interpolation requires at least three points on each processed branch.
The code also rejects non-positive or non-finite `bandrescale`, total weight,
normalization, and representative energies before tridiagonalization. Loaded
tables must each contain exactly `mMAX+1` finite values, have positive
representative energies and nonnegative amplitudes, and be normalized.

The flat band is interpolation-independent. Interpolation is inactive when
loading saved coefficient tables with mode `l`.

## Outputs

- `theta.dat`: retained hybridization weight.
- `xi.dat`: Wilson-chain hopping coefficients.
- `zeta.dat`: Wilson-chain on-site energies.
- `de_pos.dat`, `de_neg.dat`, `du_pos.dat`, `du_neg.dat`: optional intermediate tables.

See `test/tools/nrgchain/nrgchain1/param` for a complete parameter-file example
and `test/tools/nrgchain/nrgchain19_steffen` for a Steffen pipeline example.
