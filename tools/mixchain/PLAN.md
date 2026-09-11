# Channel-mixing discretization: port plan

Working document for porting `ChannelMixingDiscretization` (Julia,
`~/.julia/dev/ChannelMixingDiscretization`) to a C++ tool, `mixchain`, in
`tools/mixchain/`.

The method follows J.-G. Liu et al., PRB 93, 035102 (2016). The input is a
matrix hybridization Γ(ω), N×N Hermitian positive semidefinite. Each eigenvalue
branch of Γ is logarithmically discretized into a star Hamiltonian, and block
Lanczos maps the star to a Wilson chain with matrix coefficients.

## Background findings

### Most of the algorithm already exists here

| Julia | C++ in this repository |
| --- | --- |
| Representative energy ε(x) = W⁻¹[∫ₓ^{x+1} W(ν(x′))dx′] | `adapt` with `f_method=integral` |
| Fixed grid with accumulation point Δ: ν(x) = Δ + (D−Δ)Λ^{2−x} | `adapt` with `hardgap=true`, `boundary=Δ/D` |
| Interval weights w(x) = W(ν(x)) − W(ν(x+1)), x = j+z | `nrgchain` `df_pos[m]`, `df_neg[m]` with j = m+1 |
| Scalar Lanczos (star → chain) | `nrgchain` `tridiag()`, GMP `mpf` at `preccpp` bits |

These parts are new: splitting Γ into eigenvalue branches, taking eigenvectors
at the representative energies, block Lanczos, and matrix output.

### The `nrg` solver has no general matrix-chain input

`c++/read-input.hpp` reads only these coefficient blocks:

| Block | Content | Used by |
| --- | --- | --- |
| `z` | `xi`, `zeta` per coefficient channel | all symmetries |
| `X` | `xiR` and `zetaR` per channel (the Hamiltonian uses `zetaR` of channel 0 only) | QS, QSZ with `rungs=true` |
| `Z` | `delta`, `kappa` per channel | SPSU2 family, P, PP, NONE |
| `z` with 4 sets | up-up, down-down, up-down, down-up | U1 with `pol2x2=true` |
| `W` | general matrices `eps`, `t` per shell | compiled out behind `CHAIN_COEF` |

The `W` layout (`read_matrix_table` in `c++/coef.hpp`) is: max index, size1,
size2, then the matrices. `nrginit` can already write it (`makewilsonchain` in
`nrginit/initial.m`, QS only, diagonal).

Status of the `W` block (`CHAIN_COEF`):

- Added on 2023-05-10 (`a022a840` in `c++/`, `33207008` in `nrginit`), then
  wrapped in `#ifdef CHAIN_COEF` on 2023-11-16 (`f27b1bd0`). Nothing defines
  `CHAIN_COEF`.
- Only the reading side exists. No symmetry type uses `coef.eps` or `coef.t`:
  every `make_matrix` builds the shell-coupling terms from the scalar tables
  through its macros (`OFFDIAG`, `DIAG`, `OFFDIAG_MIX`, `RUNGHOP`,
  `ANOMALOUS`, `ISOSPINX`).
- The design note in `coef.hpp` says the meaning of the matrix indices depends
  on the symmetry type, and that each symmetry type should provide a view
  adapter onto the matrices. Enabling the block therefore means writing one
  adapter per symmetry type, which is an upstream change to every
  `sym-*-impl.hpp` in use.
- The comments on the members look swapped: `eps` is labelled as the
  f†_N f_{N+1} term and `t` as f†_N f_N, while `nrginit` writes `eps` from
  `zeta` (on-site) and `t` from `xi` (hopping).
- Enabling it would let `mixchain` hand `nrg` one general block instead of
  per-symmetry scalar files. It would not remove the gauge requirement: an
  adapter reads particular matrix entries, so the chain must still be in the
  form that symmetry type assumes. V would still enter through the initial
  Hamiltonian, not through `W`.

Other places that read chain files:

- `nrginit/wilson.m`: `band=manual` (theta files; with `pol2x2` it takes the
  matrix square root of the theta matrix), `band=manual_V` (`V{i}{j}.dat`),
  and `tri=manual`, `manual_nambu`, `manual_nambu_new` (`xi`, `zeta`, `xiR`,
  `zetaR`, `scdelta`, `sckappa`).
- `tools/matrix/parser.cc`: `load_discretization_sc()` reads
  `V{i}{j}{ch}.dat`, `xi`, `zeta`, `scdelta`, `sckappa`.
- Templates in `~/sources/nrgljubljana_templates` (SIAM/U1 `pol2x2`,
  SIAM/SPSU2).

The impurity coupling V enters the initial Hamiltonian, not the `z` block.

### Parts of the Julia code not to port as they are

- `makeguiding(::AdaptiveGrid)` uses `2^(2-x)` instead of `Λ^(2-x)`, so the
  adaptive grid is wrong for Λ ≠ 2. The examples use `:fixed`.
- ε(x) is computed as `Rint(x+1) - Rint(x)` from a running integral that
  starts at x = 1. Rounding error is about 1e-16 in absolute terms while the
  difference is about Λ^(−x), so the relative error grows like 1e-16·Λ^x:
  roughly 1e-7 at x = 30 and 1e-4 at x = 40 for Λ = 2. `nrgchain`'s default
  `mMAX = 2·Nmax` needs x up to about 2·Nmax.
- The SPSU2 exporter (`nrgfilesSPSU2`) keeps only the real parts of T₁₁, T₁₂,
  E₁₁, E₁₂ and assumes T₂₂ = −T₁₁, E₂₂ = −E₁₁. The QR gauge makes T₂₁ = 0, so
  the export is exact only when κ = 0. Nothing checks this.
- The input is resampled onto a log mesh before discretization, so the data is
  interpolated twice.
- Degeneracy handling is disabled (`if true`), and `@show` debug output
  remains in `evaluatebranch`.
- `weightsvalues` raises eigenvalues below 1e-5 up to 1e-5, adding weight.

Prior art: `~/sources/general_wilson/general_wilson_old/nrgchainNxN.cpp`
(Ž. Osolin, 2014–2015), an N×N chain code based on a different, fitting-based
method. Not used as a basis for this port.

## Proposed design

### Structure

- Two libraries with a star file on disk between them: star discretization
  (double precision) and block Lanczos (multiprecision).
- One executable with `nrgchain`-style modes: `s` computes and saves the star,
  `l` loads a star and tridiagonalizes, default does both.
- `adapt`'s integral method and mesh code (fixed grid, `GSOL` adaptive grid,
  `hardgap`) move to `tools/common` and are shared by `adapt` and the new tool.
- One z per run, as in `nrgchain`, `nrginit` and `instantiate`.
- An in-process API (like `NRG::Tools::NrgChain::calculate_from_params`) so
  that `instantiate` can use it later.

### Interpolation

Repository convention: one interpolant per tabulated function; point values and
cumulative weights both come from it via its exact primitive;
`density_interpolation=linear` (default) or `steffen`; cubic spline and Akima
are rejected for densities; values are held constant beyond the table ends and
a point at zero is added.

- Diagonalize Γ only at the input nodes. Each eigenvalue branch becomes a
  scalar density handled by `TabulatedDensity`, so its weights are exactly what
  `nrgchain` would compute for that branch.
- Representative energies come from the shared integral-method code.
- Eigenvectors at a representative energy come from Γ interpolated element by
  element (real and imaginary parts separately) with the same method key.
  Linear interpolation keeps Γ positive semidefinite; Steffen does not
  guarantee it for off-diagonal elements, which is harmless because only
  eigenvectors are taken from this interpolant.
- Between nodes, the eigenvalues of the interpolated matrix differ from the
  interpolated eigenvalues. This is inherent to the method and is documented.

Results will not match Julia exactly and will differ most at late sites.
Validation is against analytic limits and the scalar tools.

### Matrix form

Definitions. The impurity has N orbitals dᵢ coupled to bath states cₖ with
energies εₖ, H_hyb = Σₖ Σᵢ (Vᵢₖ dᵢ† cₖ + h.c.). Let Vₖ = (V₁ₖ, …, V_Nₖ)ᵀ be
the coupling vector of bath state k. The hybridization function and its
spectral function are N×N matrices:

$$
\Delta(z) = \sum_k \frac{V_k V_k^\dagger}{z-\epsilon_k},
\qquad
G_0(z)^{-1} = z - E_d - \Delta(z),
$$

$$
\Gamma(\omega) = -\frac{1}{2\pi i}\Big[\Delta(\omega+i0^+) - \Delta(\omega+i0^+)^\dagger\Big]
= \sum_k V_k V_k^\dagger\,\delta(\omega-\epsilon_k).
$$

In components, Γᵢⱼ(ω) = −[Δᵢⱼ(ω+i0⁺) − Δⱼᵢ(ω+i0⁺)*]/(2πi); for N = 1,
Γ = −Im Δ/π. Γ(ω) is Hermitian positive semidefinite, and
∫Γ(ω) dω = Σₖ Vₖ Vₖ†.

- Input: Γ(ω). The basis (channels, spin, Nambu) is the user's choice; the
  core does not need to know it.
- Normalization: the same as the `dos` file of `adapt` and `nrgchain`, which
  is named `Delta.dat` but holds πΓ, that is Γ without the 1/π:

  $$
  \Gamma_{\rm in}(\omega) = \pi\,\Gamma(\omega)
  = -\frac{1}{2i}\Big[\Delta(\omega+i0^+) - \Delta(\omega+i0^+)^\dagger\Big],
  $$

  which is −Im Δ for N = 1. The coupling is recovered downstream as √(θ/π)
  with θ = ∫`Delta.dat` (`gammaPolCh` in `nrginit/initial.m:296`,
  `gammapolch` in `tools/matrix/parser.cc:2158`). An overall factor in Γ
  leaves every Eₙ and Tₙ unchanged and scales only
  Θ = ∫Γ_in dω = π Σₖ Vₖ Vₖ†, so the convention matters only for V.
- Input files (D5): one two-column file (ω, value) per component,
  `Gamma_ij-re.dat` and `Gamma_ij-im.dat` for i, j = 1…N:
  - all `-re` files are required, including the diagonal;
  - off-diagonal `-im` files: all or none. None means Γ is real and the tool
    uses real arithmetic; a partial set is an error. Diagonal `-im` files are
    optional and must be zero within tolerance;
  - Hermiticity is checked by comparing Γⱼᵢ with Γᵢⱼ* relative to ‖Γ‖: an
    error above tolerance, otherwise Γ is symmetrized to (Γ + Γ†)/2;
  - all files must have identical ω nodes; otherwise it is an error, with no
    reinterpolation;
  - for N = 1 the input is `Gamma_11-re.dat`, an ordinary `Delta.dat` under
    another name.
- Star: per interval and sign, N orbitals with energies Eᵢ and coupling
  vectors √wᵢ·uᵢ, the discrete counterparts of Vₖ above.
- Chain: Eₙ (Hermitian), Tₙ, and the impurity coupling V to the first site,
  H_hyb = Σᵢⱼ (Vᵢⱼ dᵢ† f₀ⱼ + h.c.), with V V† = Θ/π. Real arithmetic when Γ
  is real, complex otherwise.

### Gauge and output

The core is symmetry-agnostic. Lanczos is defined only up to fₙ → Uₙfₙ at each
site, and every choice gives a unitarily equivalent chain with the same
physics, provided V transforms with the first site (V → V·U₀), because the
impurity operators are fixed.

Symmetry matters only when writing output. `nrg` cannot read a general matrix
chain (the `W` block is compiled out). Each symmetry type reads a few scalar
tables and rebuilds the chain matrices using a hard-coded form, so the chain
must be written in the gauge where it has that form. Otherwise the extra
entries are dropped and the model is silently wrong.

Example: a normal metal in Nambu form, ψ = (f↑, f↓†),
Γ(ω) = diag(ρ(ω), ρ(−ω)). The polar gauge gives Tₙ = diag(ξ, +ξ); SPSU2
expects diag(ξ, −ξ). The two differ by f↓,n → (−1)ⁿ f↓,n, which is physically
harmless but makes spin up and spin down hop with different signs, breaking the
spin-symmetric form SPSU2's matrix elements assume.

| Target | Form `nrg` imposes | Gauge needed |
| --- | --- | --- |
| Uncoupled channels (`z`) | diagonal Eₙ, Tₙ | none; Γ must be diagonal, and the canonical gauge keeps it so |
| U1 `pol2x2` | general real 2×2 Tₙ, real symmetric Eₙ | any real gauge |
| QS/QSZ rungs (`X`) | general real Tₙ, real symmetric Eₙ | any real gauge |
| SPSU2 family (`Z`) | Tₙ = [[ξ, κ], [κ*, −ξ]], Eₙ = [[ζ, δ], [δ*, −ζ]] | a specific rule; exists only if Γ has the Nambu particle-hole structure |

Not gauge freedom, and therefore important:

- Realness: real Γ must give real coefficients. The polar and QR gauges both
  keep real problems real.
- Orientation conventions, for example whether `xiUPDO` is T₁₂ or T₂₁.
  Swapping them transposes Tₙ, which is a different model. These are pinned
  from the `nrg` source (phase 0).

Design:

- Canonical gauge for the core output (D3): Hermitian square-root (polar)
  gauge. V = (Θ/π)^{1/2} is Hermitian positive semidefinite, with Θ the sum
  of vᵢvᵢ† over the star coupling vectors vᵢ = √wᵢ·uᵢ (input normalization),
  and every Tₙ is Hermitian positive semidefinite. The gauge is unique, keeps
  symmetric problems in symmetric form (decoupled channels stay diagonal, real
  stays real), gives positive ξ for N = 1, and matches `nrginit`'s
  `band=manual` `pol2x2` route. The QR gauge stays available for comparison
  with Julia.

  The name comes from the polar decomposition A = W·P (W with orthonormal
  columns, P = (A†A)^{1/2} Hermitian positive semidefinite), the matrix
  analogue of z = e^{iθ}|z|. Each Lanczos step factors the residual block
  R = H·Qₙ − Qₙ·Eₙ − Qₙ₋₁·Tₙ₋₁† as R = Qₙ₊₁·Tₙ. QR takes Tₙ upper triangular;
  the polar gauge takes Tₙ = (R†R)^{1/2} and Qₙ₊₁ = R·(R†R)^{−1/2}, which is
  Löwdin symmetric orthonormalization: the orthonormal block closest to R.
  The first step does the same with the stacked star coupling vectors,
  giving V ∝ Θ^{1/2}. Since the factorization involves no preferred basis, a
  unitary rotation of the channels rotates every Eₙ, Tₙ and V the same way,
  which is why symmetric problems stay in symmetric form.
- General output in the `W` block layout, independent of any target.
- One writer per output format (D2). A writer transforms the canonical chain
  to its target form where needed, checks that every entry it drops has the
  value the form implies, and fails otherwise. First version:
  - uncoupled channels: `xi{ch}`, `zeta{ch}`, `theta{ch}` (needed for the
    N = 1 and decoupled-channel tests);
  - U1 `pol2x2` template: `xi1..4`, `zeta1..4`, `V{i}{j}{ch}`;
  - SPSU2 template: `xi`, `zeta`, `scdelta`, `sckappa`, `V`.

  Later, if needed: QS/QSZ rungs (`X` block).
- Nambu gauge rule: a separate work item that only the SPSU2 writer depends on.
  It must handle κ ≠ 0 and complex δ; the reference is Mathematica `tri=sc2` in
  `nrginit/wilson.m:511-638`.

### Numerical points

1. Eigenvalue crossings (D4). With sorted branches, a crossing inside an
   interval can make both representative energies pick the same eigenvector,
   so one direction gets both weights and the other none. Example: diagonal Γ
   whose ρ₁ and ρ₂ cross; the true chain is decoupled but the sorted-branch
   chain gets off-diagonal terms.

   Default: tracked branches.
   - Diagonalize Γ at each input node. Going outward from ω = 0 on each
     frequency branch, assign labels at node k+1 by the permutation that
     maximizes Σᵢ |⟨vᵢ(ωₖ)|v_π(i)(ωₖ₊₁)⟩|². N ≤ 4, so all N! permutations
     can be checked.
   - Eigenvalues degenerate within a tolerance at a node: rotate the
     degenerate subspace to best match the previous node's vectors
     (projection followed by Löwdin orthonormalization) before assigning.
   - Eigenvectors at a representative energy E: diagonalize Γ̃(E) and label
     the vectors by overlap with the tracked vectors at the bracketing node.
   - Branches whose representative energies coincide within a tolerance take
     their vectors from one diagonalization, so exactly degenerate branches
     (for example Γ ∝ identity) get an orthonormal set.
   - `sorted` stays available for comparison with Julia.

   Limitation: tracking removes label swaps at exact crossings, but not errors
   from eigenvectors that rotate quickly within an interval (avoided
   crossings, or a grid too coarse to follow the rotation). The method assumes
   eigenvectors change little within one interval.

   Diagnostics, in both modes:
   - per interval, ‖Σᵢ wᵢ uᵢuᵢ† − ∫Γ_in‖ / ‖∫Γ_in‖ over that interval, with the
     largest value reported and a warning above a threshold. This detects
     both the sorted-branch failure and fast eigenvector rotation;
   - the frequencies where tracked order differs from sorted order, reported
     as crossings.
2. Matrix sum rule. Over the whole star, the trace of Θ = Σᵢ vᵢvᵢ† equals
   the trace of ∫Γ_in exactly, but Θ = ∫Γ_in holds only approximately. Report
   the deviation (the analogue of `nrgchain`'s checksum) and compute V from Θ,
   the star actually built.
3. Precision. Late ξₙ ~ Λ^(−n/2) are below double-precision rounding error.
   Block Lanczos in GMP (already a required dependency), templated over the
   number type with a small complex wrapper and hand-written small-block linear
   algebra (N ≤ 4), full reorthogonalization, `preccpp` in bits as in
   `nrgchain`. Double instantiation for unit tests. Eigen in double precision
   for diagonalizing Γ.
4. Rank-deficient Γ. If Γ is rank-deficient over the whole band, Θ is singular
   and Lanczos breaks down at the first step; stop with a clear message.
5. Conventions. `bandrescale` exactly as in `nrgchain`. Θ is written exactly
   as integrated, like `theta.dat`. Writers that output V as a physical
   amplitude (the templates' `coefV`) apply √(1/π), as `matrix` does for
   scalar chains. `nrginit`'s `band=manual_V` expects √π times the physical
   amplitude.

## Validation

- N = 1 reproduces `adapt --integral` + `nrgchain` (to about 1e-12).
- Diagonal Γ gives one independent chain per channel, matching separate
  `nrgchain` runs.
- Γ rotated by a frequency-independent unitary gives the same chains, rotated
  (polar gauge), including a complex rotation.
- Full-length block Lanczos preserves the star spectrum.
- Real input gives exactly zero imaginary parts.
- Diagonal Γ with crossing ρ₁, ρ₂: tracked branches reproduce two
  independent `nrgchain` chains with zero off-diagonal terms; sorted branches
  trigger the sum-rule warning.
- Γ = diag(ρ, ρ): exactly degenerate branches still give an orthonormal set
  of coupling vectors.
- Γ rotated by a frequency-dependent angle (avoided crossing): the sum-rule
  diagnostic reports the deviation, and it shrinks as Λ decreases toward 1.
- End to end with `nrg`: a rotated decoupled problem through U1 `pol2x2`
  reproduces the decoupled QSZ results; an s-wave superconductor matches the
  Mathematica `tri=sc` path.
- Comparison with the Julia examples (`examples/*/E_T_matrices`) at a
  tolerance that reflects Julia's accuracy, not as a strict reference.

## Phases

0. Pin references: dump Julia stars and chains for the five examples plus a
   real 2×2 mixing case and a complex case; pin each target's index
   conventions (Tₙ vs Tₙᵀ for `xiR`, `xiUPDO`/`xiDOUP`, SPSU2 `conj_me`) from
   `c++/symmetry/sym-*-impl.hpp` and `tools/matrix/parser.cc`.
1. Refactor `adapt`'s integral method and mesh code into `tools/common`;
   existing `adapt` and `nrgchain` tests pass unchanged. Check why
   `nrgchain`'s `eps_pos`/`eps_neg` do not apply `adapt`'s `hardgap`
   rescaling.
2. Star library and `s` mode, with the N = 1, decoupled, and rotated tests.
3. Block Lanczos library and `l` mode, with unit tests in `test/unit/mixchain/`.
4. Command-line tool: options and `-v`/`-vv`/`-V`/`-h` per `tools/README.md`,
   `ConfigurationReport`, `W` output, uncoupled-channel writer,
   `test/tools/mixchain/` cases.
5. Output writers chosen in D2, with end-to-end `nrg` runs. The Nambu gauge
   rule is derived and tested here, before the SPSU2 writer.
6. Docs and packaging: tool README, `docs/docs/tools.md`, options table in
   `tools/README.md`, `ChangeLog`, tests in `recipe/meta.yaml`, EasyBuild
   sanity-check lists.

Related change outside the tool: `matrix -s` ties two things together.
`-s` is needed to read `V{i}{j}{ch}.dat` (`coefV`), but it also makes
`load_discretization_sc()` require `scdelta{ch}.dat` and `sckappa{ch}.dat`, and
`V{i}{j}{ch}.dat` for every channel even though `coefV` reads only channel 1.
The U1 template writes no `Z` block, so for U1 these files exist only to
satisfy the parser. Proposal: let `matrix` read the `V` files without
superconducting mode (or load `scdelta`/`sckappa` only in superconducting
mode), so the tool does not have to write placeholder files.

Later: adaptive mesh from a matrix weight (feed tr Γ or ‖Γ‖_F to `adapt` for
`GSOL`), `instantiate` integration, enabling `CHAIN_COEF` in `nrg`, Bauer's
discretization (`bauer.jl`), reconstruction of Γ from the chain.

## Design decisions

Status: open unless marked otherwise.

| # | Decision | Proposal | Status |
| --- | --- | --- | --- |
| D1 | Tool name | `mixchain`, in `tools/mixchain/` | decided |
| D2 | Output formats written in the first version | general `W` output, plus the files read by the SIAM/U1 (`pol2x2`) and SIAM/SPSU2 templates through their Perl `instantiate` script and `matrix -s` | decided |
| D3 | Canonical gauge of the core output | Hermitian square-root (polar); QR as option | decided |
| D4 | Branch ordering | tracked by eigenvector overlap by default; `sorted` as option; per-interval matrix sum-rule diagnostic in both modes | decided |
| D5 | Input | Γ (not Δ), normalized like `adapt`'s `dos` file; all N² components as two-column `Gamma_ij-re.dat` / `Gamma_ij-im.dat` files with Hermiticity and grid checks | decided |
| D6 | Scope of Bauer's method and reconstruction | both after the first version; Bauer's method as a second star generator feeding the `l` mode | decided |
| D7 | Coordination with upstream | discuss name and output format with the maintainer early | open |
