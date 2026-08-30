# `broaden`

`broaden` converts raw, discrete spectral weights produced by NRG calculations
into a continuous spectrum. It can average several $z$-shifted spectra, apply
one of several primary broadening kernels, and optionally perform a final
Gaussian or Fermi-derivative convolution.

The finite-temperature hybrid kernels used here are implemented in
`c++/broadening.hpp` and shared with the main NRG executable. The two programs
select the same legacy kernel by default. The standalone tool additionally
provides the normalization-conserving `-n` mode.

## Usage

```text
broaden [options] <name> <Nz> <alpha> <T> [omega0_ratio]
```

The positional arguments are:

| Argument | Meaning |
| --- | --- |
| `name` | Name of the binary file containing raw delta-peak energies and weights. |
| `Nz` | Number of spectra to average. Must be at least one. |
| `alpha` | Positive width of the logarithmic broadening kernel. It is also the width of the primary Gaussian selected by `-g`. |
| `T` | Positive temperature. It sets the crossover and optional final-convolution scales. |
| `omega0_ratio` | Optional positive ratio defining $\omega_0=\mathtt{omega0\_ratio}\,T$. See [Crossover scale](#crossover-scale). |

Without `-o`, the input files are read as

```text
1/<name>, 2/<name>, ..., <Nz>/<name>.
```

With `-o` and `Nz=1`, `name` is read directly from the current working
directory. Peaks at identical frequencies are combined and their weights are
divided by `Nz`. This is $z$-averaging, not normalization to unit spectral
weight.

The primary output is always `spec.dat`, with two text columns containing the
output frequency and broadened value. Option `-c` additionally writes
`cumulative.dat`.

## Input format

The input is a native binary sequence of `double` values. By default, each
record contains

```text
frequency weight
```

This is the format written by the main executable when `savebins=true` and
`reim=false`, including the default output of a complex calculation. When the
main executable writes a spectrum with `reim=true`, each record contains

```text
frequency real_weight imaginary_weight
```

Use `-2` to broaden `real_weight` or `-3` to broaden `imaginary_weight` from
this three-double format. Binary files are architecture-native and contain no
header or byte-order marker.

## Broadening definition

Write the raw spectrum as a sum of delta peaks,

$$
A_{\mathrm{raw}}(\epsilon)
=\sum_i w_i\,\delta(\epsilon-\epsilon_i).
$$

The tool evaluates the broadened spectrum on its output mesh as

$$
A(\omega)=\sum_i w_i K(\omega,\epsilon_i),
$$

where $\omega$ is an output frequency and $\epsilon_i$ is the position of an
input delta peak.

### Modified log-Gaussian

For nonzero frequencies of the same sign, the logarithmic kernel is

$$
L_\alpha(\omega,\epsilon)=
\frac{1}{\alpha |\omega|\sqrt{\pi}}
\exp\left[
-\left(
\frac{\ln|\omega/\epsilon|}{\alpha}-\frac{\alpha}{4}
\right)^2
\right].
$$

It is zero for opposite signs and when $\epsilon=0$. For a nonzero peak
frequency, the shift
$\gamma=\alpha/4$ makes the kernel symmetric under
$\omega\leftrightarrow\epsilon$. On the full frequency axis,

$$
\int_{-\infty}^{\infty}L_\alpha(\omega,\epsilon)\,d\omega=1.
$$

### Low-energy Gaussian

The hybrid kernel uses

$$
G_{\omega_0}(\omega,\epsilon)=
\frac{1}{\omega_0\sqrt{\pi}}
\exp\left[-\left(\frac{\omega-\epsilon}{\omega_0}\right)^2\right].
$$

This Gaussian has unit area and standard deviation $\omega_0/\sqrt{2}$. Its
width convention differs from that of the optional final Gaussian convolution.

### Crossover function

The interpolation between the two kernels is controlled by

$$
h(x)=
\begin{cases}
1, & |x|\geq\omega_0,\\[2mm]
\displaystyle
\exp\left[-\left(
\frac{\ln(|x|/\omega_0)}{\alpha}
\right)^2\right], & 0<|x|<\omega_0,\\[2mm]
0, & x=0.
\end{cases}
$$

The two hybrid modes differ only in whether $x$ is the output frequency
$\omega$ or the delta-peak frequency $\epsilon$.

## Hybrid kernel modes

| Selection | Crossover argument | Kernel | Full-axis normalization |
| --- | --- | --- | --- |
| Default | Output frequency $\omega$ | $K_{\mathrm{out}}=h(\omega)L_\alpha+[1-h(\omega)]G_{\omega_0}$ | Not generally conserved |
| `-n` | Peak frequency $\epsilon$ | $K_{\mathrm{peak}}=h(\epsilon)L_\alpha+[1-h(\epsilon)]G_{\omega_0}$ | Conserved analytically for the unshifted full-axis kernel |

### Default output-frequency mode

The default kernel is

$$
K_{\mathrm{out}}(\omega,\epsilon)=
h(\omega)L_\alpha(\omega,\epsilon)
+[1-h(\omega)]G_{\omega_0}(\omega,\epsilon).
$$

Because the mixing coefficient depends on the integration variable,

$$
\int_{-\infty}^{\infty}K_{\mathrm{out}}(\omega,\epsilon)\,d\omega
=1+
\int_{-\infty}^{\infty}h(\omega)
\left[L_\alpha(\omega,\epsilon)-G_{\omega_0}(\omega,\epsilon)\right]d\omega,
$$

which is generally not one. This is the historical default because it has
produced smoother spectra in conjunction with the self-energy trick. Keeping
this mode as the default preserves existing `broaden` and main-NRG results.

### Normalization-conserving mode

Option `-n` selects the kernel proposed by Weichselbaum and von Delft,

$$
K_{\mathrm{peak}}(\omega,\epsilon)=
h(\epsilon)L_\alpha(\omega,\epsilon)
+[1-h(\epsilon)]G_{\omega_0}(\omega,\epsilon).
$$

Here the mixing coefficient is constant with respect to $\omega$, so

$$
\begin{aligned}
\int_{-\infty}^{\infty}K_{\mathrm{peak}}(\omega,\epsilon)\,d\omega
&=h(\epsilon)\int L_\alpha\,d\omega
+[1-h(\epsilon)]\int G_{\omega_0}\,d\omega\\
&=h(\epsilon)+[1-h(\epsilon)]=1.
\end{aligned}
$$

For a peak exactly at $\epsilon=0$, the implementation has $h(0)=0$ and uses
the unit-area Gaussian directly.

The `-n` option selects this formula. It does not measure `spec.dat` and divide
the result by its numerical integral. Finite meshes and subsequent operations
can therefore still change the estimated output weight.

## Crossover scale

When `omega0_ratio` is present,

$$
\omega_0=\mathtt{omega0\_ratio}\,T.
$$

When it is omitted, `broaden` uses

$$
\omega_0=10^{-9}T,
$$

which effectively disables the Gaussian crossover at ordinary nonzero mesh
frequencies. In that case both hybrid modes reduce almost everywhere to the
modified log-Gaussian.

This default differs from the main NRG executable. In the main executable,

$$
\omega_0=
\begin{cases}
\mathtt{omega0}, & \mathtt{omega0}>0,\\
\mathtt{omega0\_ratio}\,T, & \mathtt{omega0}<0,
\end{cases}
$$

An explicit `omega0=0` is rejected. The default parameters give
$\omega_0=T$. To reproduce that scale with the standalone tool, pass `1`
explicitly as `omega0_ratio`.

## Other broadening modes

### Primary Gaussian

Option `-g` bypasses the hybrid kernel and uses

$$
K_g(\omega,\epsilon)=
\frac{1}{\alpha\sqrt{\pi}}
\exp\left[-\left(\frac{\omega-\epsilon}{\alpha}\right)^2\right].
$$

In this mode `-n` has no effect because no crossover function is evaluated.

### Final Gaussian convolution

Option `-f gamma` applies a second pass after primary broadening. Its standard
deviation is $\sigma=\mathtt{gamma}\,T$, and its kernel is

$$
C_G(x,y;\sigma)=
\frac{1}{\sqrt{2\pi}\sigma}
\exp\left[-\frac{(x-y)^2}{2\sigma^2}\right].
$$

Unlike $G_{\omega_0}$ and $K_g$, this kernel uses the conventional standard
deviation parameterization.

### Final Fermi-derivative convolution

Option `-x gamma` applies a convolution of width
$\sigma=\mathtt{gamma}\,T$ with

$$
C_{\mathrm{FD}}(x,y;\sigma)=
\frac{1}{2\sigma[1+\cosh((x-y)/\sigma)]}
=\frac{1}{4\sigma}\operatorname{sech}^2
\left(\frac{x-y}{2\sigma}\right).
$$

If both `-f` and `-x` are given, the final Gaussian pass is applied first and
the Fermi-derivative pass second. Each pass updates only interior mesh points
with $|x|<100\sigma$; endpoints and points outside that range retain the result
of the preceding stage.

## Accumulation point

Option `-a a` shifts the logarithmic accumulation points to $\pm a$. Its
intended range is $0\leq a<\mathtt{broaden\_max}$. In that range, for the
log-Gaussian term only, the tool replaces $(\omega,\epsilon)$ by

$$
(\widetilde\omega,\widetilde\epsilon)=
\begin{cases}
(\omega-a,\epsilon-a), & \omega>a\ \text{and}\ \epsilon>a,\\
(\omega+a,\epsilon+a), & \omega<-a\ \text{and}\ \epsilon<-a,\\
(\omega,\epsilon), & \text{otherwise}.
\end{cases}
$$

The Gaussian and crossover function continue to use the original frequencies.
For a positive accumulation point, the generated logarithmic mesh approaches
$+a$ from above and $-a$ from below rather than approaching zero.

The current command-line parser does not reject accumulation values outside
the intended range. Their mesh and coordinate-shift behavior should not be
relied upon.

The full-axis normalization proof for `-n` assumes the unshifted kernel
($a=0$). A nonzero accumulation point changes the log-Gaussian coordinates and
the represented output domain, so `-n` should not be interpreted as an exact
normalization guarantee in that case.

## Options

| Option | Meaning |
| --- | --- |
| `-h` as the only argument | Print command-line help and exit. |
| `-v` | Print the resolved configuration and input, mesh, total-weight, and kernel-area diagnostics to standard error. |
| `-vv` | Enable very verbose diagnostics. |
| `-V` or `--version` | Print the project version and exit immediately. |
| `-m min` | Set the smallest geometric mesh scale. Default: `1e-7`. |
| `-M max` | Set the largest mesh frequency. Default: `2.0`. |
| `-r ratio` | Set the ratio between geometric mesh scales. Must exceed one. Default: `1.01`. |
| `-o` | For `Nz=1`, read `name` directly instead of `1/name`. |
| `-2` | Read three-double records and broaden the second column. |
| `-3` | Read three-double records and broaden the third column. |
| `-n` | Use the peak-frequency, normalization-conserving hybrid kernel. |
| `-s` | Print raw spectral sums weighted by fermionic and bosonic thermal factors. |
| `-c` | Write the cumulative raw spectral weight to `cumulative.dat`. |
| `-g` | Use the primary Gaussian kernel of width `alpha`. |
| `-f gamma` | Apply a final Gaussian convolution with $\sigma=\mathtt{gamma}\,T$. |
| `-x gamma` | Apply a final Fermi-derivative convolution with $\sigma=\mathtt{gamma}\,T$. |
| `-a a` | Shift logarithmic accumulation points to $\pm a$. Default: `0`. |
| `-l cutoff` | Zero input weights with $\lvert\epsilon\rvert<\mathtt{cutoff}$. |
| `-h cutoff` with other arguments | Zero input weights with $\lvert\epsilon\rvert>\mathtt{cutoff}$. |
| `-P` | Keep only nonnegative input frequencies. |
| `-N` | Keep only negative input frequencies. |
| `-A` | Generate only the positive output mesh. |
| `-B` | Generate only the negative output mesh. |
| `-L filename` | Load the output mesh from the first column of a two-column text file. |

## Examples

Broaden one file directly with the legacy kernel and the same default crossover
scale as the main executable:

```sh
broaden -o spec.bin 1 0.3 0.01 1
```

Use the normalization-conserving kernel:

```sh
broaden -n -o spec.bin 1 0.3 0.01 1
```

Average eight files named `1/spec.bin` through `8/spec.bin`:

```sh
broaden -n spec.bin 8 0.2 0.01 0.5
```

Evaluate on a custom mesh:

```sh
broaden -n -L mesh.dat -o spec.bin 1 0.3 0.01 1
```

## Numerical normalization

The analytical statements above concern integration over the complete real
axis. The numerical integral reported by

```text
Estimated weight (trapezoidal rule)=...
```

can differ from the raw weight for several independent reasons:

- `broaden_min` and `broaden_max` truncate the output domain.
- `-A` and `-B` retain only one side of the spectrum.
- Gaussian weight can cross zero or extend beyond the mesh boundaries.
- A nonzero accumulation point changes the logarithmic coordinates and mesh domain.
- A loaded mesh may be too coarse or too narrow to resolve a kernel.
- The final convolution passes use discrete quadrature, retain the endpoints, and do not divide by the sampled kernel area.
- Filtering options intentionally remove raw spectral weight.

Custom meshes must be ordered by increasing frequency for integration and final
convolution. They should not contain exactly zero because the implemented
log-Gaussian expression contains $1/|\omega|$; the generated geometric mesh
never includes zero.

Verbose mode evaluates several individual kernels on the selected mesh and
prints their trapezoidal areas. These are diagnostics only; the output is never
renormalized automatically.

## Relationship to the main NRG executable

The mathematical primitives are shared in `c++/broadening.hpp`, but the two
programs retain separate spectrum accumulation, mesh, and file-I/O code.

| Program and mode | Hybrid crossover | Accumulation passed to kernel | Default $\omega_0$ |
| --- | --- | --- | --- |
| Main NRG executable | Output frequency | Zero | $T$ |
| `broaden` default | Output frequency | Value from `-a` | $10^{-9}T$ unless the ratio is supplied |
| `broaden -n` | Peak frequency | Value from `-a` | $10^{-9}T$ unless the ratio is supplied |

Thus, for zero accumulation and equal $\alpha$ and $\omega_0$, the default
standalone kernel is the same kernel used by the main executable. The main
executable remains fixed to the historical output-frequency mode.

## Reference

- Andreas Weichselbaum and Jan von Delft, "Sum-rule Conserving Spectral
  Functions from the Numerical Renormalization Group," *Physical Review
  Letters* **99**, 076402 (2007),
  [arXiv:cond-mat/0607497](https://arxiv.org/abs/cond-mat/0607497).
