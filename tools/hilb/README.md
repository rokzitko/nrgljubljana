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

## Usage

`hilb` has three positional-argument modes.

### One complex argument

```text
hilb [options] x y
```

This evaluates the transform at

$$
z=(x+\mathtt{shiftx})+i(y+\mathtt{shifty}),
$$

where `shiftx` and `shifty` are set by `-x` and `-y`. Without `-G`, the
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

The frequency labels in corresponding rows must agree within $10^{-6}$. The
imaginary self-energy is forced to be no greater than $-10^{-8}$, and the
transform is evaluated at

$$
z=\omega-\operatorname{Re}\Sigma(\omega)
  +i[-\operatorname{Im}\Sigma(\omega)],
$$

followed by any shifts requested with `-x` and `-y`.

Without `-G`, both real and imaginary results are divided by $-\pi$ before
being written to `reaw.dat` and `imaw.dat`. With `-G`, the raw real and
imaginary parts of $H_n$ are written. For `n=0` this is the usual local Green
function transform. For `n>0` it is the corresponding energy-weighted
transform. The `-o` option is not used in this mode; the two output paths are
the final positional arguments.

## Options

| Option | Meaning |
| --- | --- |
| `-h` | Print command-line help. |
| `-d FILE` | Read a tabulated density of states from `FILE`. Otherwise use the built-in semicircular density. |
| `-n N` | Multiply the density by $E^N$. `N` must be a nonnegative integer. Default: `0`. |
| `-G` | Print or write both parts of the raw complex transform instead of the mode-specific default. |
| `-v` | Print additional information to standard output. |
| `-s SCALE` | Set the scale of the built-in semicircular density, with $B=1/\mathtt{SCALE}$. |
| `-B B` | Set the half-bandwidth of the built-in semicircular density, with $\mathtt{SCALE}=1/B$. |
| `-o FILE` | Write results to `FILE` in the single-point and argument-file modes. |
| `-x DX` | Add `DX` to the real part of every transform argument. |
| `-y DY` | Add `DY` to the imaginary part of every transform argument. |

Use either `-s` or `-B`; if both are present, the last one takes effect. They
configure only the built-in density and do not rescale data loaded with `-d`.
The shifts affect $z$, not the integration energy $E$ or the factor $E^n$.
Verbose diagnostics are written to standard output and can therefore be mixed
with data unless `-o` is used.

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

The default half-bandwidth is $B=1$. Internally, `scale=1/B` and the same
density is evaluated as
$2\,\mathtt{scale}\sqrt{1-(E\,\mathtt{scale})^2}/\pi$.

### Tabulated density

A file selected with `-d` contains two whitespace-separated numeric columns:

```text
E rho(E)
```

Records are sorted by energy after loading. The implementation constructs a
GSL natural cubic spline (`gsl_interp_cspline`) through the supplied density
values and returns zero outside the tabulated interval
`[Emin,Emax]`. Quadrature is performed over the symmetric enclosing interval

$$
[-B,B],\qquad B=\max(|E_{\min}|,|E_{\max}|).
$$

The data are not normalized automatically. In verbose mode, `hilb` reports
the integral of the interpolated density and only checks that it is finite.
At least three points with distinct energies are required by the cubic spline.
Duplicate energies, headers, comments, and extra columns are not supported by
the token-based reader. Cubic interpolation can overshoot between data points
and does not preserve positivity.

For an energy-weighted transform, interpolation is performed first and the
power is applied to each interpolated value during quadrature:

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

The real and imaginary parts of the transform are integrated separately. The
algorithm switches between direct adaptive quadrature and a
singularity-subtracted form according to the magnitude of $y=\operatorname{Im}z$.

### Direct quadrature

When $|y|\ge 10^{-3}$, the defining integral is evaluated directly. With
$z=x+iy$, the two real integrands are

$$
\operatorname{Re}\frac{g(E)}{z-E}
=\frac{g_r(E)(x-E)+g_i(E)y}{(x-E)^2+y^2},
$$

$$
\operatorname{Im}\frac{g(E)}{z-E}
=\frac{-g_r(E)y+g_i(E)(x-E)}{(x-E)^2+y^2}.
$$

Each is integrated over `[-B,B]` with GSL `gsl_integration_qag` using the
15-point Gauss-Kronrod rule.

### Small imaginary part

When $|y|<10^{-3}$, direct quadrature would need to resolve a narrow structure
near $E=x$. The implementation instead subtracts the complete weighted
numerator at that point:

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
\operatorname{Re}[g(x)Q]=aQ_R-bQ_I,
\qquad
\operatorname{Im}[g(x)Q]=aQ_I+bQ_R.
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

### Quadrature settings

The settings are fixed in the executable:

| Setting | Value |
| --- | ---: |
| GSL routine | `gsl_integration_qag` |
| Local rule | 15-point Gauss-Kronrod |
| Workspace subdivision limit | 1000 |
| Absolute tolerance | $10^{-14}$ |
| Relative tolerance | $10^{-10}$ |
| Direct/subtracted threshold | $|\operatorname{Im}z|=10^{-3}$ |

There are no command-line controls for these settings. The C++
`hilbert_transform` interface exposes the direct/subtracted threshold as
`lim_direct`; the GSL tolerances remain fixed. The GSL global error handler is
disabled. The default integration wrapper returns GSL's result without
reporting its error estimate or a nonzero status code.

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

- `Im(z)=0` is unsupported. The small-imaginary-part path divides by `|y|`
  and is not a direct principal-value integrator.
- Large powers can overflow or underflow in $E^n$ and can amplify cancellation
  in the integral.
- The fixed threshold $10^{-3}$ is expressed in the same energy units as the
  input and is not scaled with the bandwidth.
- Singularity subtraction assumes that $g(E)$ is sufficiently smooth near
  $E=x$. Discontinuous densities are difficult cases.
- At exactly $x=\pm B$ in the small-imaginary-part branch, the current boundary
  special case retains only the analytic constant contribution.
- The argument and DOS readers stop when numeric extraction fails; they do not
  support descriptive headers or comments.
