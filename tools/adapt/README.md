# `adapt`

`adapt` calculates the logarithmic discretization functions used by the numerical renormalization group (NRG). It supports fixed and adaptive meshes and treats the positive and negative frequency branches separately.

## Discretization functions

In units of the half-bandwidth, the mesh points and representative energies are

$$
\epsilon(x) = g(x)\Lambda^{2-x}, \qquad
E(x) = f(x)\Lambda^{2-x}.
$$

Here, $x=j+z$, where $j=1,2,3,\ldots$ is the interval index and $z\in[0,1]$ is the twist parameter. The lower-case $\epsilon(x)$ is the guiding function that defines interval boundaries. The upper-case representative energy $E(x)$, called `Epsilon` in program diagnostics, is the energy assigned to that interval.

The output files contain the coefficient functions, not the physical energies:

- `GSOL.dat` and `GSOLNEG.dat` contain $g(x)$ for adaptive positive and negative meshes.
- `FSOL.dat` and `FSOLNEG.dat` contain $f(x)$ for positive and negative representative energies.

Downstream programs reconstruct $\epsilon(x)$ and $E(x)$ using the powers of $\Lambda$ shown above. This file format is the same for both representative-energy algorithms.

## Representative-energy algorithms

Let

$$
W(\omega)=\int_0^\omega \rho(u)\,du
$$

be the cumulative hybridization weight on one frequency branch, and let

$$
w(x)=W(\epsilon(x))-W(\epsilon(x+1))
$$

be the weight of one discretization interval.

### Differential-equation method

The default and historically used algorithm solves

$$
\frac{dE(x)}{dx}=-\frac{w(x)}{\rho(E(x))}, \qquad E(1)=1.
$$

Internally, `adapt` evolves $f(x)$ with a fourth-order Runge-Kutta method, uses a second-order Heun step for local error estimation, and reduces the step size when necessary. The relevant controls are `dx_fine`, `dx_fast`, `xfine`, `allowed_error`, and `max_subdiv`.

Select this method by omitting `f_method` or by setting

```ini
f_method=ode
```

The default path invokes the original ODE calculation and output-grid loop without changing their floating-point operations. Existing parameter files therefore retain their previous output.

### Integral method

The representative energy is also given directly:

$$
E(x)=W^{-1}\left[\int_x^{x+1}W(\epsilon(x'))\,dx'\right].
$$

This follows by writing $\rho(E)dE/dx=dW(E)/dx=-w(x)$ and integrating from $x$ to infinity. It avoids propagating the representative energy as an initial-value problem.

The implementation normalizes $W$ by $W(1)$, evaluates the integral with GSL
CQUAD, and obtains $W^{-1}$ by monotonic bisection. Normalization does not
change the result. If $W$ has a plateau caused by an exactly zero-density
interval, the upper edge of that plateau is used as the generalized inverse.
By default CQUAD uses `epsabs=0` and `epsrel=allowed_error`; the parameter
default is `allowed_error=1e-10`. The reported `max_error` is the largest
estimated absolute CQUAD error for the normalized integral. `allowed_error`
must be positive and finite even when `--epsrel` overrides the CQUAD value,
because it can still control the adaptive ODE stage used to construct the mesh.

Select this method in the parameter file with

```ini
f_method=integral
```

or override the parameter file on the command line with `--integral`.

## Usage

```text
adapt [-h|--help] [-v|-vv] [--flat GG] [--integral] [--epsabs VALUE]
      [--epsrel VALUE] [--workspace-limit N]
      [--gsl-error-policy ignore|warn|fail] [P|N] [param_filename]
```

- `P`, `POS`, or `POSITIVE` selects positive frequencies. This is the default.
- `N`, `NEG`, or `NEGATIVE` selects negative frequencies.
- `param_filename` selects the parameter file. The default is `param`.
- `--flat GG` uses a constant positive hybridization $\rho(\omega)=GG$ without reading `Delta.dat`.
- `--integral` selects the integral representative-energy method, overriding `f_method=ode` if present.
- `--epsabs VALUE` sets the absolute CQUAD tolerance. Default: `0`.
- `--epsrel VALUE` overrides the relative CQUAD tolerance. Default: `allowed_error`.
- `--workspace-limit N` sets the CQUAD workspace capacity; `N` must be at least `3`. Default: `1000`.
- `--gsl-error-policy POLICY` selects `ignore`, `warn`, or `fail` for CQUAD failures. Default: `fail`.
- `-v` writes the resolved configuration to standard error.
- `-vv` additionally enables very verbose diagnostics.
- `-V` or `--version` prints the project version and exits immediately.
- `-h` or `--help` prints the command synopsis.

Options and positional arguments may be given in either order. Examples:

```sh
adapt P
adapt N custom.param
adapt --integral P
adapt --flat 0.01 --integral N custom.param
```

Input and output paths are relative to the working directory. Sample parameter files are under `test/tools/adapt`.

## Integral numerical controls

The four numerical command-line controls apply only to the integral method;
they do not alter the ODE calculation. Values for `--epsabs` and `--epsrel`
must be finite and nonnegative and must form a combination accepted by GSL.
The workspace limit must be at least three and small enough for GSL's internal
allocation sizes.

`--epsrel` overrides CQUAD only. It does not change the parameter-file
`allowed_error`: that value remains the ODE step-error threshold and, when
`--epsrel` is omitted, the default CQUAD relative tolerance. Similarly,
`--epsabs` defaults to zero independently of `allowed_error`.

CQUAD chooses its quadrature rules internally and has no selectable QAG rule.
Consequently, `adapt` does not accept `--quadrature-rule`.

The error policy applies when CQUAD returns a nonzero status or a nonfinite
result or error estimate. `ignore` uses the returned result without a
diagnostic, `warn` reports the failure to standard error and continues, and
`fail` aborts with a nonzero exit status. Input errors and later checks for
invalid cumulative weights or representative energies remain fatal under
every policy.

## Numerical comparison

Only the ODE and integral representative-energy algorithms are compared here. Both calculations use the same fixed, non-adaptive mesh (`adapt=false`, so $g(x)=1$). All results use $\Lambda=2$.

### Analytic flat-band benchmark

For a flat band, $\rho(\omega)=\rho_0$, the cumulative weight is $W(\omega)=\rho_0\omega$. On the fixed mesh, $\epsilon(x)=\Lambda^{2-x}$ for $x\ge2$, so the representative energy is known exactly and is independent of $\rho_0$:

$$
E_{\rm exact}(x)
= W^{-1}\left[\int_x^{x+1}W(\epsilon(x'))\,dx'\right]
= \frac{1-\Lambda^{-1}}{\ln\Lambda}\Lambda^{2-x}.
$$

For $\Lambda=2$ this gives

$$
E_{\rm exact}(x)=\frac{2^{1-x}}{\ln 2}, \qquad
f_{\rm exact}(x)=\frac{1}{2\ln 2}=0.7213475204444817.
$$

The benchmark used `adapt --flat 1 P` for ODE and `adapt --flat 1 --integral P` for the integral method, with `outputstep=1/64`, `xfine=5`, `dx_fine=1e-5`, `dx_fast=1e-4`, and `allowed_error=1e-10`. Energies are in units of the half-bandwidth $D$. The final column is normalized to the exact result.

| $x$ | Mesh scale $\epsilon(x)/D=2^{2-x}$ | Exact $E(x)/D$ | ODE $E(x)/D$ | Integral $E(x)/D$ | $|E_{\rm ODE}-E_{\rm integral}|/E_{\rm exact}$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $10$ | $3.906250000000\times10^{-3}$ | $2.817763751736\times10^{-3}$ | $2.817763751453\times10^{-3}$ | $2.817763751736\times10^{-3}$ | $1.00567\times10^{-10}$ |
| $20$ | $3.814697265625\times10^{-6}$ | $2.751722413805\times10^{-6}$ | $2.751722130434\times10^{-6}$ | $2.751722413805\times10^{-6}$ | $1.02979\times10^{-7}$ |
| $30$ | $3.725290298462\times10^{-9}$ | $2.687228919731\times10^{-9}$ | $2.686945548992\times10^{-9}$ | $2.687228919731\times10^{-9}$ | $1.05451\times10^{-4}$ |
| $40$ | $3.637978807092\times10^{-12}$ | $2.624246991925\times10^{-12}$ | $2.340876252217\times10^{-12}$ | $2.624246991925\times10^{-12}$ | $1.07982\times10^{-1}$ |

At these points no difference between the integral and analytic values was detected at 16-digit output precision (relative error below $10^{-15}$). The ODE result carries an approximately constant absolute offset of $2.83\times10^{-13}D$. Its relative effect is negligible at $x=10$, but grows as the physical energy scale decreases and reaches about 10.8% at $x=40$.

## References

- Rok Zitko, "Adaptive logarithmic discretization for numerical renormalization group methods", *Computer Physics Communications* **180**, 1271-1276 (2009).
- Rok Zitko and Thomas Pruschke, "Energy resolution and discretization artefacts in the numerical renormalization group", *Physical Review B* **79**, 085106 (2009).

## Acknowledgements

The `adapt` solver was developed during the author's stay at the Institute for Theoretical Physics, University of Goettingen, Germany. Fruitful discussions with Prof. Thomas Pruschke, computer support by GWDG, and support by the German Science Foundation through SFB 602 are acknowledged.
