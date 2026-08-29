# `hilb`

`hilb` evaluates energy-weighted Hilbert transforms of a density of states,

$$
H_n(z)=\int_{-B}^{B}\frac{E^n\rho(E)}{z-E}\,dE,
\qquad z=x+iy.
$$

The power $n$ is a nonnegative integer selected with `-n`. Its default is
`0`, which gives the original transform $H_0(z)$. The denominator convention
is `z-E`, so for a positive imaginary part of $z$ and a nonnegative effective
density $E^n\rho(E)$, the imaginary part of $H_n$ is negative. Odd powers on
a band containing negative energies do not satisfy this sign condition.

The command-line tool uses a real density of states. The low-level C++
functions in `hilb.hpp` can also transform a complex function
$\rho(E)=\rho_r(E)+i\rho_i(E)$.

All numeric output uses `std::numeric_limits<double>::max_digits10`
significant digits so that finite `double` values can be read back without
losing precision.

## Usage

`hilb` has three positional-argument modes.

### One complex argument

```text
hilb [options] x y
```

This evaluates the transform at

$$
z=(x+\Delta_x)+i(y+\Delta_y),
$$

where $\Delta_x$ and $\Delta_y$ are set by `-x` and `-y`. Without `-G`, the
program prints one value, `Im H_n(z)`. With `-G`, it prints the complete
complex result in C++ complex-number format:

```text
(Re(H_n),Im(H_n))
```

Example:

```sh
hilb -G -n 2 0.2 0.01
```

Use `--` before a negative positional `x` so that `getopt` does not interpret
it as an option:

```sh
hilb -n 1 -- -0.2 0.01
```

### File of complex arguments

```text
hilb [options] inputfile
```

Each input record consists of three whitespace-separated numbers:

```text
label x y
```

Input is line-oriented: each nonblank physical line must contain exactly one
three-field record. Blank lines are skipped. Partial records, extra fields,
comments, nonnumeric fields, and nonfinite values are errors reported with the
filename and line number.

The label is copied to the output. Without `-G`, each output row is

```text
label Im(H_n)
```

With `-G`, each row is

```text
label Re(H_n) Im(H_n)
```

The `-o FILE` option writes these data to a file instead of standard output.
The program's two `# hilb ...` identification lines remain on standard output.

### DMFT self-energy conversion

```text
hilb [options] resigma.dat imsigma.dat reaw.dat imaw.dat
```

The two input files contain rows of the form

```text
omega ReSigma
omega ImSigma
```

Each nonblank physical line must contain exactly one two-field record. Blank
lines are skipped, while comments, partial records, extra fields, nonnumeric
fields, and nonfinite values are errors. The two files must contain the same
number of data records; otherwise `hilb` reports the first unpaired record and
exits with failure.

The frequency labels in corresponding rows must agree within the tolerance set
by `-f`, whose default is $10^{-6}$. The imaginary self-energy is capped as

$$
\mathrm{Im}\,\Sigma\leftarrow
\min(\mathrm{Im}\,\Sigma,-c),
$$

where the positive clipping magnitude $c$ is set by `-c`. Define the smallest
positive normal `double` as

$$
m_{\mathrm{normal}} = 2.2250738585072014\times10^{-308}.
$$

The default clipping magnitude is

$$
c_{\min}=\sqrt{m_{\mathrm{normal}}}
\simeq 1.4916681462400413\times10^{-154}.
$$

Here $m_{\mathrm{normal}}$ is `std::numeric_limits<double>::min()`. This is the
smallest default whose square remains a normal nonzero `double`, which is
required by the current small-imaginary-part formulas. Values of `-c` below
$c_{\min}$ are rejected. The transform is evaluated at

$$
z=\omega-\mathrm{Re}\,\Sigma(\omega)
  +i[-\mathrm{Im}\,\Sigma(\omega)],
$$

followed by any shifts requested with `-x` and `-y`. If one or more input
values are changed by the cap, `hilb` writes one summary to standard error with
the number changed, the cap, and the first affected frequency and value. Since
the shifts are applied after clipping, a shift that makes the final
$|\mathrm{Im}\,z|$ smaller than the safe minimum is rejected.

Without `-G`, both real and imaginary results are divided by $-\pi$ before
being written to `reaw.dat` and `imaw.dat`. With `-G`, the raw real and
imaginary parts of $H_n$ are written. For `n=0` this is the usual local Green
function transform. For `n>0` it is the corresponding energy-weighted
transform. The `-o` option is invalid in this mode; the two output paths are
the final positional arguments.

## Options

| Option | Meaning |
| --- | --- |
| `-h`, `--help` | Print command-line help. |
| `-V` | Print `hilb 2026.09`. |
| `-d FILE` | Read a tabulated density of states from `FILE`. Otherwise use the built-in semicircular density. |
| `-i METHOD`, `--interpolation METHOD` | Interpolate a tabulated density with `linear`, `cspline`, `akima`, or `steffen`. Default: `cspline`. |
| `-n N` | Multiply the density by $E^N$. `N` must be a nonnegative integer. Default: `0`. |
| `-G` | Print or write both parts of the raw complex transform instead of the mode-specific default. |
| `-v` | Print additional information to standard output. |
| `-s SCALE` | Set the positive scale $s$ of the built-in semicircular density, with $B=1/s$. |
| `-B B` | Set the positive half-bandwidth $B$ of the built-in semicircular density, with $s=1/B$. |
| `-o FILE` | Write results to `FILE` in the single-point and argument-file modes. |
| `-x DX` | Add `DX` to the real part of every transform argument. |
| `-y DY` | Add `DY` to the imaginary part of every transform argument. |
| `-c C` | Set the positive clipping magnitude for $\mathrm{Im}\,\Sigma$ in DMFT mode. Default: $c_{\min}\simeq1.4917\times10^{-154}$. |
| `-t T` | Set the direct/singularity-subtracted threshold for the built-in density. Default: $10^{-3}$. |
| `-a A`, `--epsabs A` | Set the absolute QAG tolerance for the built-in density. Default: $10^{-14}$. |
| `-r R`, `--epsrel R` | Set the relative QAG tolerance for the built-in density. Default: $10^{-10}$. |
| `--workspace-limit N` | Set the positive QAG workspace subdivision limit for the built-in density. Default: `1000`. |
| `--quadrature-rule RULE` | Select built-in-density Gauss-Kronrod rule `15`, `21`, `31`, `41`, `51`, or `61`. Default: `15`. |
| `--gsl-error-policy POLICY` | Handle built-in-density QAG failures with `ignore`, `warn`, or `fail`. Default: `warn`. |
| `-f F` | Set the absolute frequency-label tolerance in DMFT mode. Default: $10^{-6}$. |

Use either `-s` or `-B`; if both are present, the last one takes effect. They
configure only the built-in density and do not rescale data loaded with `-d`.
The shifts affect $z$, not the integration energy $E$ or the factor $E^n$.
Verbose diagnostics are written to standard output and can therefore be mixed
with data unless `-o` is used.

The short options `-a` and `-r` are retained as legacy aliases for the uniform
long options `--epsabs` and `--epsrel`.

Malformed options, invalid numeric arguments, and positional argument counts
other than 1, 2, or 4 produce a nonzero exit status. Help and version requests
exit successfully without requiring positional arguments.

## Density Of States

### Built-in semicircular density

The default is the normalized semicircular, or Bethe-lattice, density

$$
\rho_B(E)=
\begin{cases}
\displaystyle\frac{2}{\pi B^2}\sqrt{B^2-E^2}, & |E|<B,\\
0, & |E|\ge B.
\end{cases}
$$

The default half-bandwidth is $B=1$. Internally, the scale is $s=1/B$ and the
same density is evaluated as

$$
\rho_B(E)=\frac{2s}{\pi}\sqrt{1-(Es)^2}.
$$

### Tabulated density

A file selected with `-d` contains two whitespace-separated numeric columns:

```text
E rho(E)
```

Each nonblank physical line must contain exactly one two-field record. Blank
lines are skipped. Headers, comments, partial records, extra fields,
nonnumeric values, and nonfinite values are rejected with a filename and line
number.

Records are sorted by energy after loading. The implementation constructs one
GSL interpolant through the supplied density values, materializes its
intervalwise polynomial coefficients, and reuses them for normalization and
every transform. `linear` is piecewise linear and requires at least two points;
`cspline` is a natural cubic spline and requires at least three points; `akima`
is a local piecewise-cubic interpolant and requires at least five points;
`steffen` is a monotonicity-preserving piecewise-cubic interpolant and requires
at least three points. The default is `cspline`. The interpolant is zero outside
the tabulated interval `[Emin,Emax]`, which lies within the symmetric enclosing
interval

$$
[-B,B],\qquad B=\max(|E_{\min}|,|E_{\max}|).
$$

The data are not normalized automatically. In verbose mode, `hilb` reports the
analytic integral of the interval polynomials and checks that it is finite.
Energies must be finite and strictly increasing after sorting; duplicate
energies are rejected. `cspline` and `akima` can overshoot between samples and
do not preserve positivity. `steffen` stays within each interval's endpoint
range and therefore preserves positivity of nonnegative samples, at the cost
of a possibly discontinuous second derivative and a more conservative shape
near extrema.

For an energy-weighted transform, interpolation is performed first and each
interval polynomial is multiplied by the power exactly:

$$
g(E)=E^n S_\rho(E).
$$

This is intentionally different from constructing a spline through the
weighted samples $E_i^n\rho(E_i)$.

## Numerical Algorithm

Write the effective, possibly complex numerator as

$$
g(E)=E^n\rho(E)=g_r(E)+ig_i(E).
$$

### Tabulated densities

For a density selected with `-d`, and for the tabulated C++ interfaces, the
selected GSL interpolant is represented as a polynomial on every input
interval. The transform of each weighted polynomial is evaluated analytically
using polynomial division and complex logarithms. A global subtraction of the
interpolant at $\mathrm{Re}\,z$ stabilizes arguments close to the real axis;
compensated summation combines the intervals, and a moment expansion avoids
loss of precision far from the tabulated support.

This path performs no adaptive quadrature. The `-t`, QAG tolerance, workspace,
quadrature-rule, and GSL-error-policy options are validated for compatibility
but have no numerical effect with `-d`. No GSL integration workspace is
allocated. The transform still requires a nonzero imaginary part as described
under [Identities And Limitations](#identities-and-limitations).

### Built-in and callable densities

For the built-in semicircular density and arbitrary callable C++ interfaces,
the real and imaginary parts of the transform are integrated separately. The
algorithm switches between direct adaptive quadrature and a
singularity-subtracted form according to the magnitude of $y=\mathrm{Im}\,z$.

#### Direct quadrature

When $|y|\ge T$, the defining integral is evaluated directly, where $T$ is set
by `-t` and defaults to $10^{-3}$. With $z=x+iy$, the two real integrands are

$$
\mathrm{Re}\,\frac{g(E)}{z-E}
=\frac{g_r(E)(x-E)+g_i(E)y}{(x-E)^2+y^2},
$$

$$
\mathrm{Im}\,\frac{g(E)}{z-E}
=\frac{-g_r(E)y+g_i(E)(x-E)}{(x-E)^2+y^2}.
$$

Each is integrated over `[-B,B]` with GSL `gsl_integration_qag` using the
selected Gauss-Kronrod rule. The default is the 15-point rule.

#### Small imaginary part

When $|y|<T$, direct quadrature would need to resolve a narrow structure near
$E=x$. The implementation instead subtracts the complete weighted numerator
at that point:

$$
H_n(z)=
\int_{-B}^{B}\frac{g(E)-g(x)}{z-E}\,dE+g(x)Q(z),
$$

where the constant part is integrated analytically,

$$
Q(z)=\int_{-B}^{B}\frac{dE}{z-E}=Q_R+iQ_I,
$$

$$
Q_R=\frac{1}{2}\log\frac{(B+x)^2+y^2}{(B-x)^2+y^2},
$$

$$
Q_I=\arctan\frac{x-B}{y}-\arctan\frac{x+B}{y}.
$$

For complex $g(x)=a+ib$, the analytic contribution is assembled as

$$
\mathrm{Re}[g(x)Q]=aQ_R-bQ_I,
\qquad
\mathrm{Im}[g(x)Q]=aQ_I+bQ_R.
$$

The remaining integral is rescaled with

$$
W=\frac{E-x}{|y|},
$$

split into its positive and negative parts, and transformed again with

$$
W=\pm e^r,\qquad dE=|y|e^r\,dr.
$$

This logarithmic coordinate supplies more integration points near the narrow
feature. When $x$ is inside the band, the lower limit $r=-\infty$ is truncated
to `-36.8`, for which $e^r$ is approximately $10^{-16}$.

### QAG settings

The executable uses these defaults:

| Setting | Value |
| --- | ---: |
| GSL routine | `gsl_integration_qag` |
| Local rule (`--quadrature-rule`) | `15` |
| Workspace subdivision limit (`--workspace-limit`) | `1000` |
| Absolute tolerance (`-a`, `--epsabs`) | $10^{-14}$ |
| Relative tolerance (`-r`, `--epsrel`) | $10^{-10}$ |
| GSL error policy (`--gsl-error-policy`) | `warn` |
| Direct/subtracted threshold (`-t`) | $\lvert\mathrm{Im}\,z\rvert=10^{-3}$ |

The rule can be `15`, `21`, `31`, `41`, `51`, or `61`; its number is the number
of Kronrod points used on each adaptive subinterval. Higher-order rules can be
more efficient for smooth integrands, while lower-order rules can better focus
adaptive subdivision around local difficulties. The workspace limit bounds
the number of stored subintervals. Tolerances must be finite and nonnegative
and must form a combination accepted by GSL.

The callable C++ `hilbert_transform` interface exposes `lim_direct`, `epsabs`,
and `epsrel`. An overload accepting an `integrator&` lets callers select the
workspace, rule, and policy and reuse the workspace. The executable uses one
workspace for all requested built-in-density transforms.

The GSL global error handler is disabled. A nonzero QAG status or a nonfinite
result or error estimate is handled by `--gsl-error-policy`: `ignore` silently
retains GSL's returned result, `warn` retains it and writes at most one summary
line to standard error with counts and the first failure, and `fail` aborts
with a nonzero exit status. Input errors and unrelated runtime errors remain
fatal under every policy.

## Identities And Limitations

If

$$
\mu_k=\int_{-B}^{B}E^k\rho(E)\,dE,
$$

then the weighted transforms satisfy the useful recurrence

$$
H_{n+1}(z)=zH_n(z)-\mu_n.
$$

This also shows why $H_n$ is not simply $z^nH_0$.

Numerical limitations include:

- Values with $|\mathrm{Im}\,z|<c_{\min}$ are rejected. Squaring a smaller
  value can underflow in the callable small-imaginary-part formulas. `hilb`
  does not expose a real-axis principal-value mode.
- Large powers can overflow or underflow in $E^n$ and can amplify cancellation
  in the integral.
- The clipping magnitude and callable direct/subtracted threshold are absolute
  energy scales and are not scaled with the bandwidth. The absolute QAG
  tolerance has the units of the transformed integral; the relative tolerance
  is dimensionless.
- Callable-path singularity subtraction assumes that $g(E)$ is sufficiently
  smooth near $E=x$. Discontinuous callable densities are difficult cases.
- Descriptive headers and comments are not part of the numeric input grammar.
