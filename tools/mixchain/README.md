# mixchain

Discretization of a **matrix** hybridisation function Γ(ω) and its mapping onto a Wilson chain, following
J.-G. Liu, D. Wang and Q.-H. Wang, PRB 93, 035102 (2016).

`adapt` and `nrgchain` handle a scalar hybridisation function. `mixchain` handles the N×N Hermitian positive
semidefinite case, where the channels mix: each eigenvalue branch of Γ is discretized logarithmically into a star
Hamiltonian, and the star is then mapped to a chain with matrix coefficients.

**Status: under development.** The star stage is implemented; the chain stage is not. See `PLAN.md`.

## Input

One two-column file (ω, value) per component of Γ, named `<dos_prefix>_ij-re.dat` and `<dos_prefix>_ij-im.dat` for
i, j = 1…`channels`, with `dos_prefix` defaulting to `Gamma`. Each file covers both negative and positive
frequencies, as `Delta.dat` does for the scalar tools.

- All `-re` files are required, including the diagonal.
- The off-diagonal `-im` files are required either all or none. None means that Γ is real and the whole calculation
  runs in real arithmetic. A diagonal `-im` file is optional and must vanish.
- Every file must be tabulated on one and the same frequency grid: the input is never reinterpolated.
- Γ_ji is checked against conj(Γ_ij) relative to the largest element of Γ, and Γ is then symmetrized to
  (Γ + Γ†)/2.

**Normalization.** As for the `dos` file of `adapt` and `nrgchain`, the input is π times the spectral function of
the hybridisation function:

    Gamma_in(omega) = pi Gamma(omega) = -Im Delta(omega + i0+)   for N = 1,

with Δ(z) = Σ_k V_k V_k†/(z − ε_k) and Γ(ω) = Σ_k V_k V_k† δ(ω − ε_k). An overall factor in Γ leaves every chain
coefficient unchanged and scales only the hybridisation weight Θ = ∫Γ_in dω, so the convention matters only for the
impurity coupling V.

**Band.** `bandrescale` maps the band edge to 1: ω → ω/bandrescale and Γ → Γ·bandrescale, which leaves ∫Γdω
unchanged. The discretization mesh only ever reaches |ω| = 1, so any weight tabulated beyond the band edge is
discarded, and where the input stops short of the edge the density is continued at its last tabulated value. Both
are reported per diagonal element when they occur.

## Parameters

Read from the `[param]` block of the parameter file, with the names of `adapt` and `nrgchain` wherever the quantity
is the same.

| parameter | default | meaning |
| --- | --- | --- |
| `Lambda` | — | discretization parameter |
| `z` | 1 | twist parameter |
| `mMAX` | `2*Nmax` | largest interval index; the star has 2·`channels`·(`mMAX`+1) levels |
| `channels` | 1 | the dimension N of Γ |
| `dos_prefix` | `Gamma` | prefix of the input files |
| `bandrescale` | 1 | band rescaling |
| `adapt` | false | adaptive mesh (see below) |
| `hardgap`, `boundary` | false, 0 | accumulation point of the mesh, as a fraction of the rescaled band edge |
| `density_interpolation` | `linear` | `linear` or `steffen` |
| `mesh_weight` | `frobenius` | weight function of the adaptive mesh: `frobenius` or `trace` |
| `branch_ordering` | `tracked` | `tracked` or `sorted` |
| `allowed_error` | 1e-10 | error control of the integral method |

**Mesh.** The interval boundaries are ε(x), the interval of index m being [ε(z+m+2), ε(z+m+1)]. With `adapt=false`
this is the usual ε(x) = Λ^(2−x). With `adapt=true` it is ε(x) = W⁻¹(Λ^(2−x)), where W is the normalized cumulative
of a scalar weight function — ‖Γ(ω)‖_F or tr Γ(ω), see `mesh_weight` — so that the intervals follow the structure of
Γ. This is *not* the adaptive mesh of `adapt`, which solves an ODE for a guiding function g(x) and writes it to
`GSOL.dat`; `mixchain` reads no such file and solves no differential equation. Both forms are followed by the
`hardgap` rescaling, though with cleanly gapped data the adaptive mesh already accumulates at the gap edge by
itself.

**Branches.** Γ is diagonalized at the input nodes only. Going outward from ω = 0, the eigenvectors at neighbouring
nodes are matched by overlap, so that a branch keeps its identity through a crossing of eigenvalues; degenerate
subspaces are rotated to match the previous node. `branch_ordering=sorted` labels by descending eigenvalue instead,
which is what exposes the crossing artifact. Both modes report where the two labellings diverge.

## The star file

`star.dat`, written by the `s` mode and read by the `l` mode, describes

    H = sum_k E_k c_k^dag c_k + sum_k sum_i ( v_{k,i} d_i^dag c_k + h.c. )

with one row per bath level k.

Lines beginning with `#` are comments. Exactly one of them is read back: the line carrying `channels=`, whose
whitespace-separated `key=value` pairs form the header. Unknown keys are ignored, which is why the diagnostics and
the column legend can be comments.

Header keys:

| key | meaning |
| --- | --- |
| `channels` | the dimension N of Γ, and the number of components of every coupling vector |
| `mMAX` | the largest interval index; the file holds 2·`channels`·(`mMAX`+1) rows |
| `z` | the twist parameter of the mesh |
| `Lambda` | the discretization parameter |
| `bandrescale` | the rescaling applied to Γ; the energies are in the rescaled band, whose edge is 1 |
| `complex` | 1 if the coupling vectors are complex, 0 if real. It fixes the number of columns. |

Columns of a data row:

| column | meaning |
| --- | --- |
| `m` | interval index, 0…`mMAX`; a larger m lies closer to the Fermi level |
| `sign` | the frequency branch the level came from, `+` for ω > 0 and `-` for ω < 0 |
| `a` | the eigenvalue branch of Γ it came from, 0…`channels`−1 |
| `E` | the representative energy of that interval and branch, carrying the sign of its frequency branch |
| `v_i` | the coupling vector in the **channel** basis; for `complex=1` each component is a pair `Re_v_i Im_v_i` |

`a` is a label from the branch tracking, not a channel index: a branch is an eigenvector direction of Γ, which in
general points across several channels and rotates with ω, so a single row usually has all N components nonzero.
The two coincide only for a Γ that is diagonal and stays so.

The coupling vector is v = √w·u, with w the weight of that branch over that interval and u its eigenvector at E, in
the normalization of the input. Neither is renormalized, so Σ_k v_k v_k† is the hybridisation weight Θ of the star
as built.

`m`, `sign` and `a` are not used to construct the chain, which sees only the pairs (E, v). They are written so that
the file can be read by a person and checked when loaded: the row count, the ranges of `m` and `a`, and the
agreement of the sign column with the sign of E. The diagnostics are comments and are not read back, being
properties of Γ rather than of the star.

Example, for `channels=2` with real data:

```
# mixchain star
# channels=2 mMAX=2 z=1 Lambda=2 bandrescale=1 complex=0
# max_interval_deviation=0 at_omega=0
# m sign a E v1 v2
0 + 0 0.721347520444481703 0.547722557505166113 0
0 + 1 0.721347520444481703 0 0.316227766016837933
1 + 0 0.360673760222240851 0.387298334620741689 0
```

## Diagnostics

- Per interval, the matrix sum rule ‖Σ_a w_a u_a u_a† − ∫Γ‖ / ‖∫Γ‖, with the largest value and where it occurred.
  It detects both a mislabelled branch and eigenvectors that rotate too fast for the interval to resolve.
- Θ = Σ_k v_k v_k† against ∫Γ over the range the mesh covers. The trace agrees by construction; the off-diagonal
  elements only approximately.
- The frequencies at which the tracked and the sorted branch orderings diverge.
- Whether the mesh reaches below the innermost tabulated frequency, in which case the density there is the constant
  continuation of the input — exact for a flat band, an approximation for anything with structure at low frequency.
- How many intervals hold no node of the input tabulation, and from which frequency downwards. There the star
  follows the interpolant between two tabulated points rather than the data. It happens wherever the mesh resolves
  more finely than the input: at the bottom of the band, and around an accumulation point set by `hardgap` or found
  by the adaptive mesh — where, for a density that diverges at a gap edge, it changes the fall-off of the
  coefficients.
- The largest CQUAD error estimate of the integral method.
