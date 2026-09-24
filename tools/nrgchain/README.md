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

## Scalar tridiagonalization backend

The opt-in parameter `tridiag_method=rkpw` selects the scalar RKPW
star-to-chain kernel. The default, `tridiag_method=lanczos`, keeps the legacy
GMP calculation. Names are case-sensitive and other values are rejected, even
in save-only mode. `instantiate` uses the same setting for its in-process
Wilson-chain generation; the C++ runtime uses it when `tri=cpp`.

For example, add this to `[param]` and run `nrgchain -v`, or compare the
backends on the same saved star with `nrgchain s` followed by `nrgchain l`
after changing only `tridiag_method`:

```ini
tridiag_method=rkpw
```

`preccpp` defaults to `2000` and specifies **GMP precision in bits**, not decimal
digits. It must be an integer greater than `10` for `lanczos`. With `rkpw`, it
is still parsed as a nonnegative integer but is unused; `0` is allowed and
verbose diagnostics mark it inactive. The tools accept values through the
maximum signed C++ `int`. RKPW uses fixed floating-point arithmetic, not GMP,
and changing `preccpp` does not increase its precision.

Both backends receive the same normalized amplitudes and representative
energies. RKPW preserves shell order, alternating positive and negative
energies from high to low shells. Exactly zero amplitudes are removed and
exactly equal energies are combined. If the resulting finite star has `K`
distinct supported energies, `Nmax+1` must not exceed `K`. At equality, the
last hopping is exactly zero; an overlong request fails rather than padding
or silently shortening the chain. The calculation finishes before RKPW opens
`xi.dat` and `zeta.dat`, so this failure does not truncate existing coefficient
files.

The backend does not change `theta`, normalization, file precision, or
coefficient conventions: both `xi` and `zeta` have `Nmax+1` entries and
`xi[n]` couples sites `n` and `n+1`. `rescalexi=true` divides hoppings only by
`SCALE(n+1)`; on-site energies are never rescaled this way. `bandrescale`
multiplies both coefficient arrays afterward. Runtime `tri=cpp` applies
`bandrescale` but has no `rescalexi` step.
RKPW rejects a nonfinite scaled coefficient or a nonzero coefficient that
underflows to zero during scaling, rather than silently emitting a broken
chain. An exact terminal zero remains valid.

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

The flat band is interpolation-independent. Generated `band=flat` stars require
`adapt=false`; `adapt=true` is explicitly rejected rather than retaining the
old accidental behavior. Interpolation and mesh adaptation are inactive when
loading saved star tables with mode `l`.

### Small shell amplitudes

`TabulatedDensity::normalized_amplitude(lower, upper, root_total_weight)`
takes the **square root** of the normalization weight as its third argument,
not the weight itself. Normal shell masses retain the existing
`sqrt(mass)/root_total_weight` arithmetic. For subnormal or rounded-zero
masses, the fallback integrates local density means in root space and combines
pieces using scaled mantissa/exponent arithmetic. It does not require
multiprecision or extra exponent range from `long double`.

Exactly zero-density intervals stay exactly zero. Representable complete
amplitudes survive even when individual mass products underflow; genuinely
unrepresentable complete amplitudes are errors. A nonpositive computed local
mean on nonzero support is also an error, not empty support. This is not a
universal relative-accuracy guarantee for underresolved polynomial averages or
input energies. If the retained total is subnormal or rounds to zero, branch roots
are combined before the final rounding of `theta`; the physical `theta`
must still be a finite, strictly positive `double`.

## Bounded hard gaps

`hardgap` defaults to `false` and `boundary` to `0`. With `hardgap=true`,
`boundary=b` must be finite and in `[0,1)`. It is in normalized band units:
the physical gap edges are `+/-bandrescale*b`, not `+/-b` in the original
density units. The two branches use separate instances of the shared
`LogMesh`: on the fixed mesh, `eps(x)=1` for `x<=2` and
`eps(x)=b+(1-b)*Lambda^(2-x)` otherwise.

Nonzero-gap star generation requires `adapt=false`. For `band=flat`, the
analytic representative is `E(x)=b+(1-b)*E_flat(x)`, where `E_flat` is the
ungapped flat-band expression. With density `1/2` per branch, the retained
weight is `(1-b)*(1-Lambda^(-z-mMAX))`, including the finite star cutoff.

For a tabulated band, both `FSOL.dat` and `FSOLNEG.dat` must have at least
two `(x, f)` rows, finite strictly increasing abscissas, and positive finite `f`
values. Their **actual** ranges must cover every query from `z+1` through
`z+mMAX+1`; setting a larger `xmax` does not extend a file. Each stored node
is first converted to `E_i=f_i*Lambda^(2-x_i)`, then these **energies** are
interpolated linearly within the table. There is no extrapolation and no
additional affine gap shift. Generate compatible tables with `adapt=false`
and `adapt --integral P` / `adapt --integral N`, or `f_method=integral`;
nonzero-gap ODE and adaptive-mesh calculations are unsupported.

Each generated shell must resolve `b < lo < up <= 1`. Its representative
must satisfy `lo < E < up` for a positive amplitude; endpoints are allowed
only for an exactly zero amplitude. Collapsed shells or out-of-shell energies
fail: reduce `mMAX` or refine/regenerate the tables, rather than extending a
tail artificially. In `adapt`, premature `max_abs` termination is an error
for a nonzero gap. A validated table is published through a checked sibling
temporary and atomic rename, preserving the previous destination on failure;
this is not a crash-durability guarantee. See the [adapt guide](../adapt/README.md#bounded-hardgap).

Load mode instead validates saved energy magnitudes against `b <= E <= 1`
alongside the usual positivity, finiteness, length, and normalization checks.
It does not require density, FSOL, or GSOL files, and `adapt` is inactive;
it does not reconstruct or revalidate the original shell mesh.
`hardgap=false` or `boundary=0` leaves the ungapped mapping unchanged.

## Outputs

- `theta.dat`: retained hybridization weight.
- `xi.dat`: Wilson-chain hopping coefficients.
- `zeta.dat`: Wilson-chain on-site energies.
- `de_pos.dat`, `de_neg.dat`, `du_pos.dat`, `du_neg.dat`: optional intermediate tables.

See `test/tools/nrgchain/nrgchain1/param` for a complete parameter-file example
and `test/tools/nrgchain/nrgchain19_steffen` for a Steffen pipeline example.
`nrgchain20_rkpw_legacy` and `nrgchain20_rkpw_rkpw` independently exercise
calculate/save/load consistency, scaling, and `instantiate` for the selected
backend; finite-support termination checks are RKPW-specific, not cross-backend
comparisons. `nrgchain21_star_ranges_legacy` and `nrgchain21_star_ranges_rkpw`
add analytic/table hard-gap checks and deep, tiny-weight star tests. See the
[focused test commands](../../test/README.md#scalar-star-contracts).
