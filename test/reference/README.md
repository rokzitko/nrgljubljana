# Multiprecision transform references

This directory contains offline reference data for the canonical Cauchy
transform used by `piecewise_polynomial.hpp` and `hilb`:

\[
p(E)=\sum_k c_{i,k}u^k,\qquad
u=\frac{E-a_i}{b_i-a_i},\qquad
H_n(z)=\int_L^R\frac{E^n p(E)}{z-E}\,dE.
\]

The coefficients are local normalized-coordinate coefficients, and `p` is
zero outside `[L,R]`. Energy weighting means interpolate first and then
multiply the resulting polynomial by `E^n`. Complex references use the
principal branch of `log`; for real `p`, upper/lower conjugate arguments give
conjugate results. Rows marked `real_pv` use the Cauchy principal value inside
the support and the ordinary real integral outside it. Internal-knot PV rows
are finite because both generated polynomials are continuous. Endpoints are
intentionally excluded.

The KK references use
$-\operatorname{PV}H_0(x)/\pi$, matching `kk`'s denominator convention
$1/(E-x)$.

`quadratic` is an explicit continuous piecewise quadratic. `cubic_cspline` is
the natural cubic spline through the headerless uneven-grid data in
`generated/cubic_cspline_dos.dat`; its endpoint second derivatives vanish, so
the file is also a direct `hilb -d ... -i cspline` fixture. `legendre3` has
three exactly vanishing leading far-field moments, while `ill_conditioned`
has a small nonzero leading moment formed by cancellation of coefficients of
order `1e100`. Together they distinguish true cancellation from arithmetic
roundoff in the far-field series.

The DMFT fixture uses

\[
z=\omega-\operatorname{Re}\Sigma(\omega)
  +i[-\operatorname{Im}\Sigma(\omega)],\qquad
(\mathrm{reaw},\mathrm{imaw})=-H_0(z)/\pi.
\]

`dmft_resigma.dat`, `dmft_imsigma.dat`, `dmft_reaw.ref`, and
`dmft_imaw.ref` are direct two-column CLI input/reference files. The same
case, including `z` and the raw transform, is recorded in `dmft_mapping.tsv`.

## Generated files

All `.tsv` files are UTF-8/ASCII, tab-separated, and begin with comment lines
that name their columns. Decimal fields contain 60 significant digits (except
exact zero); calculations use 400 decimal digits.

- `polynomials.tsv`: interval bounds and `c0` through `c3`; quadratic `c3` is zero.
- `moments.tsv`: exact polynomial moments `integral E^n p(E) dE` for `n=0..3`.
- `cauchy_transforms.tsv`: complex and real-PV transform cases.
- `dmft_mapping.tsv`: the self-energy-to-transform mapping case.
- `cubic_cspline_dos.dat`: headerless tabulated DOS for the CLI.
- `dmft_*.dat` and `dmft_*.ref`: headerless DMFT CLI inputs and references.
- `kk_nonlinear_input.dat` and `kk_nonlinear.tsv`: symmetric nonlinear input and independent Akima/cspline KK values.

Normal CTest runs consume only the committed generated files. They do not run
Python and do not require mpmath.

## Regeneration

From the repository root:

```sh
python3 -m venv test/reference/.venv
test/reference/.venv/bin/python -m pip install -r test/reference/requirements.txt
test/reference/.venv/bin/python test/reference/generate.py
test/reference/.venv/bin/python test/reference/generate.py --check
```

The generator requires exactly mpmath 1.3.0, contains no random or
platform-dependent inputs, writes no timestamps, and validates analytic
results against independent 400-digit quadrature. `--check` performs a
byte-for-byte comparison without rewriting the committed data.
