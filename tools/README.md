# Numerical controls

The tabulated-data tools use the same option names where their numerical
controls overlap. Tool-specific options and input formats are documented in
each tool directory.

## Uniform options

| Option | Tools | Meaning |
| --- | --- | --- |
| `-i METHOD`, `--interpolation METHOD` | `hilb`, `kk`, `integ`, `resample` | Select `linear`, `cspline`, or `akima` interpolation. |
| `--epsabs VALUE` | `hilb`, `kk`, `integ`, `adapt` integral | Set the absolute integration tolerance. |
| `--epsrel VALUE` | `hilb`, `kk`, `integ`, `adapt` integral | Set the relative integration tolerance. |
| `--workspace-limit N` | `hilb`, `kk`, `integ`, `adapt` integral | Set the GSL workspace capacity; QAG requires `N >= 1` and CQUAD requires `N >= 3`. |
| `--quadrature-rule RULE` | `hilb`, `kk`, `integ` | Select QAG rule `15`, `21`, `31`, `41`, `51`, or `61`. |
| `--gsl-error-policy POLICY` | `hilb`, `kk`, `integ`, `adapt` integral | Select `ignore`, `warn`, or `fail`. |

`hilb` retains `-a` and `-r` as aliases for `--epsabs` and `--epsrel`.
`integ` retains `-w` as an alias for `--gsl-error-policy ignore`; its `-a`
option selects the integral of the absolute value and is not an `epsabs`
alias.

## Compatibility defaults

These defaults preserve each tool's behavior when the new controls are not
specified. A dash means that the control is not available.

| Tool or mode | Interpolation | `epsabs` | `epsrel` | Workspace | QAG rule | Error policy |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `hilb` | `cspline` | `1e-14` | `1e-10` | `1000` | `15` | `warn` |
| `kk` | `akima` | `1e-12` | `1e-8` | `1000` | `15` | `ignore` |
| `integ` | `akima` | `1e-12` | `1e-8` | `1000` | `15` | `warn` |
| `resample` | `akima` | - | - | - | - | - |
| `adapt` integral | `linear` (fixed) | `0` | `allowed_error` (`1e-10`) | `1000` | - | `fail` |

`adapt` uses its existing internal linear interpolation and has no
`--interpolation` option.

The interpolation methods are GSL `linear` (piecewise linear, at least two
points), `cspline` (natural cubic spline, at least three points), and `akima`
(local piecewise cubic, at least five points). Cubic interpolation is smooth
but can overshoot between samples and does not preserve positivity. Use
`linear` when preserving the range and positivity of nonnegative samples is
more important than smooth derivatives.

QAG applies the selected Gauss-Kronrod rule on each adaptive subinterval. The
rule number is the number of Kronrod points. Higher-order rules can be more
efficient for smooth integrands; lower-order rules can leave more work to
adaptive subdivision near local difficulties. The workspace limit bounds the
number of stored subintervals. Tolerances must be finite and nonnegative, the
workspace limit must meet the algorithm-specific minimum, and the tolerance
pair must be accepted by GSL. Limits that could overflow GSL's internal
allocations are rejected.

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

In `kk` and `integ`, the spline `Sum` check uses
`gsl_spline_eval_integ`. That operation is separate from QAG and is not
controlled by the QAG tolerances, rule, workspace, or error policy.
