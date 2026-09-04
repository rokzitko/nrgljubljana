# Floquet Model Construction With nrginit

The Mathematica initializer builds a finite one-frequency Floquet factor and
appends it to the model basis. The model owns the drive Hamiltonian and the mode
operators exported to the runtime.

## Configuration

A Floquet initialization supplies the mode cutoff, angular frequency, runtime
switch, and mode operator:

```ini
[extra]
Omega=0.234

[param]
model=model.m
symtype=SPU1
floquet=true
options=Ncut=1
ops=m
diagratio=1
```

`Ncut` is a non-negative decimal integer. `Omega` must round to a finite,
positive C++ `double`; a positive input that rounds to zero is rejected. The
same `param` file is read by the initializer and the C++ runtime.
`Omega` uses decimal-literal syntax and does not accept the initializer's
`!expression` extension.

## Floquet Basis

`floquetbasis[Ncut]` returns

\[
\left\{|-N_{\mathrm{cut}}\rangle,\ldots,|0\rangle,\ldots,
  |N_{\mathrm{cut}}\rangle\right\}.
\]

`transformtoFL[basis,{Ncut}]` forms the tensor-product basis. Existing invariant
labels are unchanged, and the Floquet index is appended as the final ket and bra
coordinate. Floquet extension is the last built-in basis stage, so spin,
orbital, and phonon factors precede it. A `hook_basis` implementation receives
the completed non-Floquet basis; the initializer applies the Floquet factor
after the hook, so hook-provided auxiliary coordinates also precede it. The
hook must leave `bvc` as a canonical occupation-vector basis: occupation
entries must remain valid, states must be nonzero, and auxiliary factors must
form one consolidated ket without unresolved `Null` coordinates.

The initializer sets the global lowercase symbol `ncut` from the exact option
token `Ncut=...`. A model that participates in this basis extension sets
`MAKEFLOQUET = 1`.

## Mode-Space Operators

The built-in constructors are:

| Mathematica expression | Operator |
| --- | --- |
| `floquetid[Ncut]` | \(I_F=\sum_m|m\rangle\langle m|\) |
| `floquetnumber[Ncut]` | \(\hat M=\sum_m m|m\rangle\langle m|\) |
| `floquetshift[k,Ncut]` | \(S_k=\sum_m|m+k\rangle\langle m|\), within the finite window |
| `floquetplus[Ncut]` | \(S_1\) |
| `floquetminus[Ncut]` | \(S_{-1}=S_1^\dagger\) |
| `floquetcos[Ncut]` | \(S_1+S_{-1}\) |
| `floquetlift[op,k,Ncut]` | \(op\otimes S_k\) |

These helpers resolve against the state on which they act and preserve all ket
coordinates preceding the final Floquet coordinate. This makes the mode number,
identity, and shifts independent of the number of earlier tensor factors.

For Fourier components \(H_k\), the finite Sambe Hamiltonian has the form

\[
H_F=\sum_k H_k\otimes S_k+\Omega I\otimes\hat M.
\]

`floquetcos` contains the sum of the two shifts and no factor of \(1/2\). Thus a
term \(A O\cos(\Omega t)\) contributes

\[
\frac{A}{2}O\otimes\texttt{floquetcos[Ncut]}.
\]

`floquetlift[op,k,Ncut]` is the direct construction for a Fourier component.
For a static operator that itself contains auxiliary ket/bra factors, use
`floquetlift[op,0,Ncut]` to extend it over the Floquet identity. Fermionic-only
static terms already preserve a trailing auxiliary ket during application.

## Three-Mode Example

For `Ncut=1`, in the ordered basis \((|-1\rangle,|0\rangle,|1\rangle)\),

\[
\Omega\hat M+g(S_1+S_{-1})=
\begin{pmatrix}
-\Omega & g & 0\\
g & 0 & g\\
0 & g & \Omega
\end{pmatrix}.
\]

Its eigenvalues are

\[
0,\qquad \pm\sqrt{\Omega^2+2g^2}.
\]

This matrix is also the analytic case used by the Floquet primitive regression
test.

## Custom Model File

A model can define the finite-space terms directly:

```mathematica
MAKEFLOQUET = 1;
snegrealconstants[g, Omega];

Hfl = g floquetlift[number[d[]], 1, ncut] +
      g floquetlift[number[d[]], -1, ncut] +
      Omega floquetnumber[ncut];
opm = floquetnumber[ncut];
opm2 = pow[opm, 2];

H = H0 + Hc + H1 + Hfl;
```

Here `g` is the coefficient of each of the two nearest-mode shifts. Static model
terms act diagonally in the appended mode coordinate.

## Mode Operator Export

The runtime truncation path reads exactly one even-parity singlet declaration
named `m`. A `modeloperators.m` file can export the model-owned definitions as:

```mathematica
Module[{t = {}},
  t = Join[t, mtSingletOp["m", opm]];
  t = Join[t, mtSingletOp["m^2", opm2]];
  t
]
```

`m^2` is an optional observable and is not used to form the truncation
criterion. The generated `m` declaration contains one Hermitian diagonal block
for every seed invariant sector. The runtime obtains the allowed mode interval
from the spectrum of those seed blocks.

The repository example is in
`test/nrginit+nrgrun/floquet_spu1/`, including its `param`, `model.m`, and
`modeloperators.m` files.

## Related Pages

- [Floquet formalism](floquet-formalism.md)
- [Floquet truncation](floquet-truncation.md)
- [nrginit workflow](nrginit-workflow.md)
- [Parameter reference](parameter-reference.md)
