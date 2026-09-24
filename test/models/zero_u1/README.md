# Zero-U1 Matsubara Contract

`ZERO` consists of two uncoupled, noninteracting Wilson chains. This fixture
uses `U1`, flat-band Z discretization, `Lambda=2`, `z=1`, `Nmax=2`, and
`T=0.01`. With `keep=5000`, all 4096 final many-body states are retained.
Each channel/spin therefore has the endpoint Green's function of an exact
three-site, zero-onsite chain.

Both backend variants run the usual strict golden comparisons for every
output except these two files:

- `spec_FDMmats_dens_A_f_u-A_f_u.dat`
- `spec_FDMmats_dens_A_f_d-A_f_d.dat`

`validate_matsubara.pl` validates those files instead. It requires five matched
three-column finite numeric records, retaining the standard strict spectral
comparison against the existing references for frequency and imaginary part:
absolute tolerance `1e-12`, relative frequency tolerance `1e-5`, and relative
spectral tolerance `2e-2`. Records are not sorted or interpolated.

Particle-hole symmetry makes every real component exactly zero. Each actual
real component must satisfy `abs(Re G) <= 1e-12`, the unchanged standard
absolute bound, rather than reproduce a reference's nonzero roundoff. For
example, the historical first down-spin real value is `-9.594024077e-13`,
while an RKPW run gives `+2.161291979e-13`. Both obey the zero bound, but their
difference exceeds it. Current legacy runs also differ from the historical
reference. Repeated RKPW runs, including initializer precision 1000, give the
same result. No golden files or tolerance overrides are changed.

An additional independent check uses the analytic flat-band hoppings, not
coefficients read from production output:

```text
t0^2 = 1 / (7 log(2)^2)
t1^2 = 18 / (217 log(2)^2)
omega_n = (2 n + 1) pi T, n = 0,...,4
G(i omega_n) = -i (omega_n^2 + t1^2)
                 / (omega_n (omega_n^2 + t0^2 + t1^2))
```

Frequency and imaginary part must also match this expression within `5e-10`
relative error, allowing the ten-significant-digit output rounding. This
checks the model and normalization independently of both chain backends and
the historical spectral references. It is intentionally specific to this
fixture, not a general relaxation of spectral comparisons.

The neutral `zero_u1_matsubara_contracts` test uses static analytic data and
rejects nonzero real values beyond the bound, wrong grids and imaginary
parts, mismatched/reordered rows, wrong column counts, nonfinite numbers,
and unrepresentable nonzero inputs. It runs no initializer or solver.

```sh
perl test/models/zero_u1/test_validate_matsubara.pl
perl test/models/zero_u1/validate_matsubara.pl test/models/zero_u1/ref /absolute/run/directory
```
