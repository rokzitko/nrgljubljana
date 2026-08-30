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
| `-i METHOD`, `--interpolation METHOD` | `hilb`, `kk`, `integ`, `resample` | Select `linear`, `cspline`, `akima`, or `steffen` interpolation. |
| `--epsabs VALUE` | `hilb`, `integ`, `adapt` integral | Set the absolute integration tolerance. |
| `--epsrel VALUE` | `hilb`, `integ`, `adapt` integral | Set the relative integration tolerance. |
| `--workspace-limit N` | `hilb`, `integ`, `adapt` integral | Set the GSL workspace capacity; QAG requires `N >= 1` and CQUAD requires `N >= 3`. |
| `--quadrature-rule RULE` | `hilb`, `integ` | Select QAG rule `15`, `21`, `31`, `41`, `51`, or `61`. |
| `--gsl-error-policy POLICY` | `hilb`, `integ`, `adapt` integral | Select `ignore`, `warn`, or `fail`. |

Configuration reports use `auto -> VALUE` when a value is derived from input
data, another setting, or an environment-dependent default. Verbosity does not
change machine-readable standard output or numerical output files.

`hilb` retains `-a` and `-r` as aliases for `--epsabs` and `--epsrel`.
`integ` retains `-w` as an alias for `--gsl-error-policy ignore`; its `-a`
option selects the integral of the absolute value and is not an `epsabs`
alias.

## Compatibility defaults

These defaults preserve each tool's behavior when the new controls are not
specified. A dash means that the control is not available.

| Tool or mode | Interpolation | `epsabs` | `epsrel` | Workspace | QAG rule | Error policy |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `hilb`, built-in DOS | - | `1e-14` | `1e-10` | `1000` | `15` | `warn` |
| `hilb`, tabulated DOS | `cspline` | - | - | - | - | - |
| `kk` | `akima` | - | - | - | - | - |
| `integ` | `akima` | `1e-12` | `1e-8` | `1000` | `15` | `warn` |
| `resample` | `akima` | - | - | - | - | - |
| `adapt` integral | `linear` (fixed) | `0` | `allowed_error` (`1e-10`) | `1000` | - | `fail` |

`adapt` uses its existing internal linear interpolation and has no
`--interpolation` option.

The interpolation methods are GSL `linear` (piecewise linear, at least two
points), `cspline` (natural cubic spline, at least three points), `akima`
(local piecewise cubic, at least five points), and `steffen`
(monotonicity-preserving piecewise cubic, at least three points). `cspline` and
`akima` can overshoot between samples and do not preserve positivity. `steffen`
stays within each interval's endpoint range and therefore preserves positivity
of nonnegative samples, but its second derivative can be discontinuous and its
shape can be conservative near extrema.

`hilb -d` and `kk` materialize the selected GSL interpolant as a polynomial on
each input interval. Ordinary integrals and Cauchy transforms of those
polynomials are then evaluated analytically, without adaptive quadrature. The
QAG controls accepted by `hilb` affect only its built-in density; they are
validated but have no numerical effect with `-d`. `kk` has no QAG controls.

For the built-in `hilb` density and for `integ`, QAG applies the selected
Gauss-Kronrod rule on each adaptive subinterval. The rule number is the number
of Kronrod points. Higher-order rules can be more efficient for smooth
integrands; lower-order rules can leave more work to adaptive subdivision near
local difficulties. The workspace limit bounds the number of stored
subintervals. Tolerances must be finite and nonnegative, the workspace limit
must meet the algorithm-specific minimum, and the tolerance pair must be
accepted by GSL. Limits that could overflow GSL's internal allocations are
rejected.

The error policies apply when GSL returns a nonzero status or a nonfinite
result or error estimate:

- `ignore` suppresses the diagnostic and continues with GSL's returned result.
- `warn` reports the failure to standard error and continues with GSL's returned result.
- `fail` stops the calculation with a nonzero exit status.

Input errors and later validity checks remain fatal under every policy.

`adapt` uses CQUAD for the integral representative-energy method. CQUAD has no
selectable QAG quadrature rule, so `adapt` has no `--quadrature-rule` option.
Its `--epsrel` overrides CQUAD only: it does not change the parameter-file
`allowed_error`, which remains the ODE error control and the default CQUAD
relative tolerance.

The `Sum` reported by `hilb -d` and `kk` is the analytic integral of the
interval polynomials. In `integ`, the spline `Sum` check uses
`gsl_spline_eval_integ`; that operation is separate from QAG and is not
controlled by the QAG tolerances, rule, workspace, or error policy.
