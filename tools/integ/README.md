# `integ`

`integ` integrates one selected quantity from a two-column tabulation. Let
$s(x)$ be the selected GSL interpolant and let its tabulated support be

$$
D=[x_{\min},x_{\max}].
$$

A successful calculation writes one unlabelled scalar to standard output using
`std::numeric_limits<double>::max_digits10` significant digits. Diagnostics and
errors are written to standard error.

## Usage

```text
integ [options] [input ...]
```

With no input arguments, or with an input argument of `-`, `integ` reads
standard input. `-` may appear at most once. If several files and/or standard
input are specified, all records form one logical table and are sorted by `x`
before interpolation.

Input is line-oriented. Blank lines and full-line comments whose first
non-whitespace character is `#` are skipped. Every other physical line must
contain exactly two whitespace-separated finite representable numbers. Inline
comments, partial or extra records, malformed fields, and values that overflow
or underflow to zero are rejected with the source and line number.

After sorting, `x` values must be strictly increasing and the total support
width must be finite. Duplicate `x` values are therefore rejected.

## Quantities

Exactly one quantity may be selected. Conflicting or repeated selectors fail;
the paired bounds together count as one selector. If none is given, the total
integral is selected.

| Selector | Quantity |
| --- | --- |
| default, `--total` | $\displaystyle \int_D s(x)\,dx$ |
| `--lower-bound A --upper-bound B` | $\displaystyle \int_{D\cap[A,B]}s(x)\,dx$; both finite bounds are required, $A\le B$, and an empty overlap gives zero |
| `-p`, `--positive` | $\displaystyle \int_{D\cap(0,\infty)}s(x)\,dx$ |
| `-n`, `--negative` | $\displaystyle \int_{D\cap(-\infty,0)}s(x)\,dx$ |
| `-a`, `--absolute` | $\displaystyle \int_D\lvert s(x)\rvert\,dx$ |
| `--negative-absolute` | $\displaystyle \int_{D\cap(-\infty,0)}\lvert s(x)\rvert\,dx$ |
| `-e`, `--energy-moment` | $\displaystyle \int_D x\,s(x)\,dx$ |
| `-f`, `--fermi` | $\displaystyle \int_D s(x)f_T(x)\,dx$ |

For $T>0$, the Fermi factor is

$$
f_T(x)=\frac{1}{1+\exp(x/T)}.
$$

Only the selected quantity is computed.

## Options

| Option | Meaning |
| --- | --- |
| `-h`, `--help` | Print help. |
| `-v` | Print the resolved configuration and diagnostics to standard error. |
| `-vv` | Also print interval-level QAG diagnostics when QAG is used. |
| `-V`, `--version` | Print the project version. |
| `-i METHOD`, `--interpolation METHOD` | Select `linear`, `cspline`, `akima`, or `steffen`. Default: `akima`. |
| `-t T`, `-T T` | Set the finite, nonnegative Fermi temperature. Default: `1e-99`. |
| `--epsabs VALUE` | Set the absolute QAG tolerance. Default: `1e-12`. |
| `--epsrel VALUE` | Set the relative QAG tolerance. Default: `1e-8`. |
| `--workspace-limit N` | Set the QAG workspace subdivision limit. Default: `1000`. |
| `--quadrature-rule RULE` | Select QAG rule `15`, `21`, `31`, `41`, `51`, or `61`. Default: `15`. |
| `--gsl-error-policy POLICY` | Select `ignore`, `warn`, or `fail`. Default: `warn`. |
| `-w` | Legacy alias for `--gsl-error-policy ignore`. |

`-a` is a quantity selector, not an absolute-tolerance option. Use `--epsabs`
for the QAG tolerance.

`linear` is piecewise linear and requires at least two points. `cspline` is a
natural cubic spline and requires at least three points. `akima` is a local
piecewise-cubic interpolant and requires at least five points. `steffen` is a
monotonicity-preserving piecewise-cubic interpolant and requires at least three
points. `cspline` and `akima` can overshoot between samples. `steffen` stays
within each interval's endpoint range, but its second derivative can be
discontinuous and its shape can be conservative near extrema.

## Numerical Algorithm

The selected interpolant is materialized once as a shared polynomial of degree
at most three on each input interval. Total, bounded, positive, negative,
absolute, negative-absolute, energy-moment, and zero-temperature Fermi
integration use analytic polynomial antiderivatives. Absolute integration
finds the real roots in every interval and splits the interval at those roots,
so the absolute value is applied after interpolation rather than to the input
samples.

Here, "analytic" or "exact" means that these modes have no quadrature error.
Interpolation error, root-finding error, and floating-point construction,
evaluation, and summation can still introduce rounding error.

For $T>0$, Fermi integration alone uses adaptive GSL QAG. It integrates every
polynomial interval separately in a dimensionless local coordinate and adds
temperature-scaled breakpoints around zero so narrow boundary layers and every
representable positive-energy tail remain resolved. Polynomial amplitude,
physical width, and the leading Fermi scale are applied outside the QAG
callback to avoid premature overflow or underflow. At $T=0$, Fermi mode maps
exactly to the analytic negative-domain integral and does not use QAG. The QAG
tolerances, workspace limit, rule, and error policy affect only
positive-temperature Fermi mode; they are still parsed and validated in other
modes. The default temperature, `1e-99`, is positive and therefore takes the
QAG path.

The QAG rule number is the number of Gauss-Kronrod points used on each adaptive
subinterval. Local error estimates are rescaled and summed, then checked against
the requested tolerance for the final result so cancellation between intervals
cannot hide a missed global target. `ignore` suppresses that failure diagnostic,
`warn` reports it and continues, and `fail` exits unsuccessfully. A nonfinite
final result, input error, or setup error remains fatal under every policy.

## Legacy Front Ends

The legacy commands are installed alongside `integ` and now invoke it with
`--interpolation linear`:

| Command | Selected quantity |
| --- | --- |
| `integrate` | total |
| `integrateab filename [A B]` | bounded; the default is `[-1,1]` |
| `integratepos` | positive |
| `integrateneg` | negative |
| `integrateabs` | absolute |
| `integratenegabs` | negative-absolute |
| `integrateeps` | energy moment $\int_D x\,s(x)\,dx$ |

These front ends preserve the command names but intentionally correct several
historical semantics:

- absolute values are taken after linear interpolation, including roots inside a segment;
- bounded and sign-restricted domains clip segments that cross a boundary instead of dropping the crossing part;
- `integrateeps` evaluates the first moment analytically from the interpolant, with no trapezoidal product approximation;
- input uses the strict parser described above, and all supplied records form one sorted table;
- results use `max_digits10` precision.
