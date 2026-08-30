# `integ`

`integ` integrates a two-column tabulation using a selected GSL interpolant and
adaptive QAG quadrature.

## Definitions

Let $s(x)$ be the selected interpolant through the input samples, with input
domain

$$
D = [x_{\min}, x_{\max}].
$$

The default output is the total integral

$$
I = \int_D s(x)\,dx
  = \int_{x_{\min}}^{x_{\max}} s(x)\,dx.
$$

For the positive- and negative-frequency modes, define

$$
D_+ = D \cap (0, \infty),
\qquad
D_- = D \cap (-\infty, 0).
$$

The corresponding integrals are

$$
I_+ = \int_{D_+} s(x)\,dx,
\qquad
I_- = \int_{D_-} s(x)\,dx.
$$

The absolute-value mode integrates the absolute value of the interpolant:

$$
I_{\mathrm{abs}} = \int_D \lvert s(x) \rvert\,dx.
$$

Thus, the absolute value is applied after interpolation, rather than to the
input samples before constructing the interpolant. The Fermi-weighted mode is

$$
I_{\mathrm{F}}(T)
= \int_D s(x) f_T(x)\,dx,
\qquad
f_T(x) = \frac{1}{1 + \exp(x/T)}.
$$

## Usage

```text
integ [options] input [-p|-n|-a|-f]
```

The program computes the total, positive-frequency, negative-frequency,
absolute-value, and Fermi-weighted integrals, then prints the selected value
with 16 significant digits.

## Options

| Option | Meaning |
| --- | --- |
| `-h`, `--help` | Print help. |
| `-v` | Print resolved configuration and verbose diagnostics to standard error. |
| `-vv` | Also print detailed integration diagnostics to standard error. |
| `-V`, `--version` | Print the project version. |
| `-t T`, `-T T` | Set the Fermi-Dirac temperature. Default: `1e-99`. |
| `-p` | Print the integral over positive `x`. |
| `-n` | Print the integral over negative `x`. |
| `-a` | Print the integral of the absolute value. |
| `-f` | Print the Fermi-Dirac-weighted integral. |
| `-i METHOD`, `--interpolation METHOD` | Select `linear`, `cspline`, `akima`, or `steffen`. Default: `akima`. |
| `--epsabs VALUE` | Set the absolute QAG tolerance. Default: `1e-12`. |
| `--epsrel VALUE` | Set the relative QAG tolerance. Default: `1e-8`. |
| `--workspace-limit N` | Set the positive QAG workspace subdivision limit. Default: `1000`. |
| `--quadrature-rule RULE` | Select `15`, `21`, `31`, `41`, `51`, or `61`. Default: `15`. |
| `--gsl-error-policy POLICY` | Select `ignore`, `warn`, or `fail`. Default: `warn`. |
| `-w` | Legacy alias for `--gsl-error-policy ignore`. |

`-a` selects an output quantity; unlike `hilb -a`, it is not an absolute
tolerance option. Use `--epsabs` for the tolerance.

`linear` is piecewise linear and requires at least two points. `cspline` is a
natural cubic spline and requires at least three points. `akima` is a local
piecewise-cubic interpolant and requires at least five points. `steffen` is a
monotonicity-preserving piecewise-cubic interpolant and requires at least three
points. `cspline` and `akima` can overshoot between samples. `steffen` stays
within each interval's endpoint range, but its second derivative can be
discontinuous and its shape can be conservative near extrema.

The QAG rule number is the number of Gauss-Kronrod points used on each adaptive
subinterval. Higher-order rules can help on smooth functions; lower-order rules
can leave more work to adaptive subdivision near local difficulties. The
workspace limit bounds the number of stored subintervals. Tolerances must be
finite and nonnegative and must form a combination accepted by GSL.

For a nonzero QAG status or a nonfinite result or error estimate, `ignore`
silently uses GSL's returned result, `warn` reports the failure to standard
error and continues, and `fail` exits unsuccessfully. Input and setup errors
remain fatal.

The verbose spline `Sum` is computed by `gsl_spline_eval_integ`, not QAG. The
QAG tolerances, rule, workspace, and error policy control the five requested
quadratures but do not control that spline integration operation.
