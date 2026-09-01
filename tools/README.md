# Numerical controls

The tabulated-data tools use the same option names where their numerical
controls overlap. Tool-specific options and input formats are documented in
each tool directory.

## Uniform options

| Option | Tools | Meaning |
| --- | --- | --- |
| `-v` | all tools | Report the resolved configuration and verbose diagnostics on standard error. |
| `-vv` | all tools | Also enable detailed diagnostics where available. |
| `-V`, `--version` | all tools | Print the NRG Ljubljana project version and exit immediately. |
| `--algorithm ALGORITHM` | `hilb`, `kk` | Select adaptive `qag` (default) or the optional analytic transform. `hilb` provides only this long form. |
| `-a ALGORITHM` | `kk` | Short alias for `kk --algorithm`; `hilb -a` instead sets `epsabs`. |
| `-i METHOD`, `--interpolation METHOD` | `hilb`, `kk`, `integ`, `resample` | Select `linear`, `cspline`, `akima`, or `steffen` interpolation. |
| `-j N`, `--jobs N` | `hilb`, `kk` | Select a positive worker count; otherwise use the first value of a fully valid `OMP_NUM_THREADS` list, or one when unset. |
| `--epsabs VALUE` | `hilb`, `kk` QAG, `integ`, `adapt` integral | Set the absolute adaptive-integration tolerance. |
| `--epsrel VALUE` | `hilb`, `kk` QAG, `integ`, `adapt` integral | Set the relative adaptive-integration tolerance. |
| `--workspace-limit N` | `hilb`, `kk` QAG, `integ`, `adapt` integral | Set the GSL workspace capacity; QAG requires `N >= 1` and CQUAD requires `N >= 3`. |
| `--quadrature-rule RULE` | `hilb`, `kk` QAG, `integ` | Select QAG rule `15`, `21`, `31`, `41`, `51`, or `61`. |
| `--gsl-error-policy POLICY` | `hilb`, `kk` QAG, `integ`, `adapt` integral | Select `ignore`, `warn`, or `fail`. |

Configuration reports use `auto -> VALUE` when a value is derived from input
data, another setting, or an environment-dependent default. Verbosity does not
change machine-readable numerical output or numerical output files. The
standard-output exceptions are the `hilb` and `kk` successful file-output
timing reports described below; `-v` adds their performance line.

`hilb --algorithm` is long-only. `hilb` retains `-a` and `-r` as aliases for
`--epsabs` and `--epsrel`, so `-a` does not select its algorithm.
`integ` retains `-w` as an alias for `--gsl-error-policy ignore`; its `-a`
option selects the integral of the absolute value and is not an `epsabs`
alias.

`integ` uses its QAG controls only for Fermi-weighted integration at positive
temperature. The controls are parsed and validated but have no numerical
effect for its analytic modes, including Fermi integration at `T=0`.

## Compatibility defaults

These are the defaults used when the controls are not specified. A dash means
that the control is not available.

| Tool or mode | Interpolation | `epsabs` | `epsrel` | Workspace | QAG rule | Error policy |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `hilb`, built-in DOS, QAG (default) | - | `1e-14` | `1e-10` | `1000` per worker | `15` | `warn` |
| `hilb`, tabulated DOS, QAG (default) | `cspline` | `1e-14` | `1e-10` | `1000` per worker | `15` | `warn` |
| `hilb`, tabulated DOS, analytic | `cspline` | - | - | - | - | - |
| `kk`, QAG (default) | `akima` | `1e-12` | `1e-8` | `1000` per worker | `15` | `ignore` |
| `kk`, analytic | `akima` | - | - | - | - | - |
| `integ` | `akima` | `1e-12` | `1e-8` | `1000` | `15` | `warn` |
| `resample` | `akima` | - | - | - | - | - |
| `adapt`, `nrgchain` density | `linear` | - | - | - | - | - |
| `adapt` integral | `linear` | `0` | `allowed_error` (`1e-10`) | `1000` | - | `fail` |

`adapt` and `nrgchain` accept `density_interpolation=linear|steffen` in the
shared parameter file and have no interpolation command-line override. Both
tools must use the same setting because `adapt` constructs representative
energies while `nrgchain` independently constructs shell weights. The legacy
Mathematica discretization remains linear.

The interpolation methods are GSL `linear` (piecewise linear, at least two
points), `cspline` (natural cubic spline, at least three points), `akima`
(local piecewise cubic, at least five points), and `steffen`
(monotonicity-preserving piecewise cubic, at least three points). `cspline` and
`akima` can overshoot between samples and do not preserve positivity. `steffen`
stays within each interval's endpoint range and therefore preserves positivity
of nonnegative samples, but its second derivative can be discontinuous and its
shape can be conservative near extrema.

`hilb -d`, analytic `kk`, and `integ` materialize the selected GSL interpolant
as a polynomial on each input interval. `hilb --algorithm analytic -d` and
analytic `kk` evaluate the corresponding Cauchy transforms analytically. Hilb's
analytic algorithm requires `-d`. Its default QAG algorithm applies to both the
built-in and tabulated DOS; tabulated QAG evaluates the selected interpolant as
zero outside its support over the symmetric enclosing interval `[-B,B]` and
integrates each polynomial interval in a local normalized coordinate.
Hilb silently ignores explicit QAG controls with `--algorithm analytic`, while
analytic `kk` rejects them.

`integ` analytically evaluates its total, bounded, positive, negative,
absolute, negative-absolute, first-moment, and zero-temperature Fermi modes.
Absolute modes locate roots in each polynomial interval before integrating.
Only positive-temperature Fermi mode uses QAG, separately on every interval
and with additional temperature-scaled breakpoints around zero.

For built-in and tabulated `hilb` QAG, QAG `kk`, and positive-temperature
`integ` Fermi mode, QAG applies the selected Gauss-Kronrod rule on each adaptive
subinterval. The rule number is the number of Kronrod points. Higher-order rules
can be more efficient for smooth integrands; lower-order rules can leave more
work to adaptive subdivision near local difficulties. The workspace limit
bounds the number of stored subintervals, independently for each parallel QAG
worker where applicable. Tolerances must be finite and nonnegative, the
workspace limit must meet the algorithm-specific minimum, and the tolerance
pair must be accepted by GSL. Limits that could overflow GSL's internal
allocations are rejected.

The error policies apply when adaptive integration fails its effective error
target. `integ` rescales and sums its interval error estimates before applying
the policy to the final Fermi result:

- `ignore` suppresses the diagnostic and continues with GSL's returned result.
- `warn` reports the failure to standard error and continues with GSL's returned result.
- `fail` stops the calculation with a nonzero exit status.

Input errors and later validity checks remain fatal under every policy.

`adapt` uses CQUAD for the integral representative-energy method. CQUAD has no
selectable QAG quadrature rule, so `adapt` has no `--quadrature-rule` option.
Its `--epsrel` overrides CQUAD only: it does not change the parameter-file
`allowed_error`, which remains the ODE error control and the default CQUAD
relative tolerance.

The `Sum` reported by `hilb -d` and analytic `kk` is the analytic integral of
the interval polynomials. QAG `kk` obtains it from the selected GSL spline. See
the [`integ` documentation](integ/README.md) for its complete input, quantity,
and algorithm contract.

For `hilb` and `kk`, an explicit `--jobs` value takes precedence over
`OMP_NUM_THREADS`. Otherwise the entire comma-separated environment value must
be a valid list of positive integers, after which its first value is used; an
unset variable selects one worker. The actual count is
`min(requested workers, output points)`. Inputs are staged before parallel
calculation, output order is deterministic, and each actual QAG worker owns a
separate GSL workspace.

Successful regular-file `hilb -o` runs, Hilb DMFT runs with two regular-file
outputs, and successful file-mode `kk` runs report `Time elapsed: ... s` on
standard output after output writes. With `-v`, a second `Performance:` line
adds wall time, process CPU time, effective parallelism, throughput, actual
worker count, and algorithm. Unavailable process CPU time is shown as
`cpu=n/a effective_parallelism=n/a`. Hilb emits no timing when numerical output
goes to standard output, including `-o /dev/stdout` or a DMFT `/dev/stdout`
destination; `kk` stream mode likewise remains numeric-only. Hilb's worker pool
is CLI-internal: no public batch API is added, and its scalar APIs emit no CLI
timing or performance output.

## Installed Integration Front Ends

`integrate`, `integrateab`, `integratepos`, `integrateneg`, `integrateabs`,
`integratenegabs`, and `integrateeps` are installed compatibility front ends to
`integ`. They force linear interpolation. `integrateab filename [A B]` defaults
to `[-1,1]`; the other names select total, positive, negative, absolute,
negative-absolute, or first-moment integration as indicated by their names.

Compared with the historical scripts, the front ends use strict two-column
parsing and one sorted logical table, clip segments at subdomain boundaries,
apply absolute value after interpolation, evaluate the first moment
analytically, and print at `max_digits10` precision.
