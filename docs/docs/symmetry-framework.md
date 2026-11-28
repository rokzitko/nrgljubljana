# Symmetry Framework

## Purpose

The symmetry layer determines how the Hilbert space is decomposed into invariant subspaces and how operators transform between them. In practice, it is one of the key extensibility points of the project.

The central abstraction is `Symmetry<S>` in `c++/symmetry.hpp`.

## What The Symmetry Backend Controls

A symmetry backend is responsible for:

- declaring the conserved quantum numbers and combination rules
- deciding which subspaces are allowed
- generating new subspaces from the previous iteration
- defining ancestor relationships for each new subspace
- constructing Hamiltonian matrix elements
- providing operator-recalculation hooks
- exposing multiplicity and prefactor functions used in observables

## Selection

The runtime chooses a concrete symmetry implementation through `mk_sym.hpp`, based on the configured `symtype` and channel count.

The symmetry types are listed in one registry, `c++/symmetries.def`. For each type it records the class
(`c++/symmetry/sym-<name>-impl.hpp`), the build set (`SYM_MORE`, `SYM_ALL`), and the properties that are needed
before a symmetry object is constructed: spin-resolved spectra (`isfield()`), the number of `<||f||>` blocks per
channel in `data`, the spin states per orbital, the supported numbers of channels, complex-only implementations, and
support for `substeps`, `polarized`, and pairing. The registry generates the factory functions and `get()` in
`mk_sym.hpp`, the traits table in `symmetry_traits.hpp`, and (in `c++/CMakeLists.txt`) one factory source per type
from `c++/sym-instance.cc.in`. `set_symmetry()` rejects a number of channels or a `substeps`/`polarized` setting that
the type does not support before its tables are loaded. `test/unit/mk_sym` checks the registry against the classes.

To add a symmetry type, write its class in `c++/symmetry/sym-<name>-impl.hpp`, add a line to `symmetries.def`, and
reconfigure.

## Interaction With The Main Flow

The symmetry backend is consulted throughout the iteration:

1. `InputData` constructs the symmetry object.
2. `SubspaceStructure` uses it to enumerate ancestors and valid new subspaces.
3. Hamiltonian construction uses it to fill block matrices.
4. Operator recalculation uses it to map old irreducible matrix elements to the new basis.
5. Thermodynamic and spectral prefactors depend on symmetry-specific multiplicity rules.

## Coefficient Tables

Each NRG step combines the states of the added Wilson-chain sites (one per channel)
with the invariant subspaces of the previous step. The coefficients that define the
new basis, the Hamiltonian blocks, and the recalculation of the operators come from
two sources:

- For the abelian types `NONE`, `P`, `U1`, `QSZ`, `SPU1`, `QSZ4`, and `SPZ4`, every combination is a
  single Fock state of the added sites. `AbelianShell::Shell` in `c++/abelian_shell.hpp`
  enumerates these states and computes all coefficients as fermion matrix elements,
  for any number of channels. `test/unit/symmetry/abelian_shell.cpp` checks that the
  result equals the Mathematica tables, which it reads from
  `c++/symmetry/coefnew/{none,p,u1,qsz,spu1}`, entry by entry. `QSZ4` and `SPZ4` have no
  Mathematica tables; they are checked through the cross-symmetry tests below.
- All other types use tables generated with Mathematica (`c++/symmetry/coefnew`,
  driven by `nrgcoef.m`) for each supported number of channels. The Hamiltonian tables
  (`*-offdiag.dat`, `*-diag.dat`, ...) and the tables of global operators are sequences
  of statements such as `OFFDIAG(i, j, ch, factor);`, whose factors are expressions in the
  quantum numbers of the subspace. Each statement is a macro defined once in
  `c++/table_terms.hpp`; it calls a method of a context object (`ShellTerms` in
  `make_matrix()`, `GlobalTerms` in `recalc_global()`) that selects the Wilson-chain
  coefficients and, with `substeps`, the channel of the added site. The recalculation
  tables (`*-doubletp.dat`, `*-spinupa.dat`, ...) are lists of rows included in the
  recalculation routines `c++/symmetry/nrg-recalc-*.hpp` through
  `Symmetry::recalc_table()` and `recalc_f_table()`; the rows are only evaluated if both
  subspaces hold states, since their factors may be undefined otherwise.

The conventions of the abelian tables are those of `nrgcoef.m`: modes are ordered by
channel with spin up before spin down; a Fock state is the product of the creation
operators of the occupied modes in increasing mode order; and in a combination
`|i, r>` of an added-site state `|i>` with a state `|r>` of the previous sites, the
added-site operators stand to the left. An operator `O` of the previous sites then
has `<i, r|O|j, r'> = delta_ij (-1)^(n_i |O|) <r|O|r'>`, with `n_i` the number of
electrons in `|i>` and `|O|` the fermion parity of `O`.

## Z<sub>4</sub> Spin Symmetry: `SPZ4` And `QSZ4`

The rhombic anisotropy $E(S_x^2-S_y^2) = (E/2)(S_+^2+S_-^2)$ breaks the conservation of
$S_z$ but changes $2S_z$ only by $\pm 4$. Terms such as s-wave pairing, $D S_z^2$, a
field $B S_z$, and isotropic or XXZ exchange conserve $2S_z$. In their presence with
$E \neq 0$, $2S_z$ modulo 4 is conserved. For spin-1/2 electrons this is a finer
symmetry than the fermion parity used by `P`:

| Type | Quantum numbers | Use instead of |
| --- | --- | --- |
| `SPZ4` | `SSZ` = $2S_z \bmod 4$ | `P` (superconducting bath with $E \neq 0$) |
| `QSZ4` | `Q` (charge), `SSZ` = $2S_z \bmod 4$ | `U1` (normal bath with $E \neq 0$) |

`f^dag_up` raises the label by 1 and `f^dag_down` by 3 (modulo 4). `nrginit` labels
the states by $2S_z+1$ before taking the modulus, as for `SPU1` and `QSZ`. The
modulo-4 quantum-number type is `mod4` in `c++/invar.hpp`, alongside `additive`,
`multiplicative`, and `mod3`. The C++ classes reuse the `SPU1` and `QSZ`
implementations with the spin label taken modulo 4.

Operators in `nrginit/operators.m` (with `spin` in `[extra]` for `model=KONDO`):

- Singlets (conserve $2S_z \bmod 4$): `SZd`, `SZ2d`, `SX2d`, `SY2d`, `Pu`, `Pd`, and
  `SZ`, `SZ2`, `SX2`, `SY2` for the impurity spin.
- Doublets: fermion operators such as `A_d`, which change $2S_z$ by $+1$ or $-1$
  (labels $+1$ and $+3$). The spectral functions are spin-resolved, as for `SPU1` and
  `QSZ`: `specd=A_d-A_d-u A_d-A_d-d`.
- Operators of type `t` (change $2S_z$ by $\pm 2$, i.e. by 2 modulo 4): `SXd`, `SYd`,
  `SPd`, `SMd`, and `SX`, `SY`, `SP`, `SM`. Their correlators are computed with
  `spect`, e.g. `spect=SX-SX` (output `spin_*_dens_SX-SX.dat`). With `P` or `U1`, the
  same operators are singlets and use `specs` (output `corr_*`).

The expectation values of `t` operators are zero by symmetry and are not written to
`custom`. `td` has no spin moments for `SPZ4` and `QSZ4`.

Example (spin-1 impurity exchange-coupled to a superconductor, with axial and rhombic
anisotropy):

```text
[extra]
spin=1
Jkondo=0.2
anD=0.05
anE=0.02
B=0.01

[param]
symtype=SPZ4
model=KONDO
variant=ANISONEW
bcsgap=0.1
...
ops=SZ SZ2 SX2 SY2 SX
spect=SX-SX
```

`test/c++sym_conventions/conv_aniso_*` checks this problem against `P` and `U1`.

## Conventions Shared By All Symmetry Types

The same Hamiltonian must give the same results for every symmetry type that can
describe it. The initial sites of the Wilson chain are set up by `nrginit`
(`HBANDonsite` and `HBANDhop` in `nrginit/hamiltonian.m`); the terms of the
sites added by the C++ iteration (`make_matrix()`) follow the same conventions.
For site $n$ of channel $\alpha$ with spin $\sigma\in\{\uparrow,\downarrow\}$:

- **On-site energy.** Without spin polarization the term is
  $\zeta_{n\alpha}(n_{n\alpha}-1)$; with `polarized=true` it is
  $\sum_\sigma \zeta_{n\alpha\sigma}(n_{n\alpha\sigma}-1/2)$. Both subtract the
  average occupancy, so the two forms coincide for
  $\zeta_{n\alpha\uparrow}=\zeta_{n\alpha\downarrow}$. Absolute energies (total
  energy, `E_fdm`, `F_fdm` in `tdfdm`, `dumpabs` levels) are therefore the same for
  all symmetry types and with or without `polarized`.
- **Pairing.** The on-site pairing term is
  $\Delta_{n\alpha} f^\dagger_{n\alpha\uparrow} f^\dagger_{n\alpha\downarrow} + \Delta^*_{n\alpha} f_{n\alpha\downarrow} f_{n\alpha\uparrow}$,
  with complex $\Delta_{n\alpha}$ allowed (`bcsgap=!0.1+0.1I`, `options=COMPLEX`).
  The phase is physical when it varies along the chain or between channels, for
  example in Josephson junctions. `SPSU2T` supports only real $\Delta$.
- **Unsupported terms are rejected.** Nonzero pairing (`delta`) coefficients are
  accepted only by `SPSU2`, `SPU1`, `SPU1LR`, `SPZ4`, `SPSU2LR`, `SPSU2T`, `SPSU2C3`, `P`,
  `PP`, and `NONE`; nonzero anomalous-hopping (`kappa`) coefficients only by
  `SPSU2`, `SPU1`, `SPU1LR`, `SPZ4`, `SPSU2LR`, and `SPSU2T`. `substeps=true` is accepted
  only by `QS`, `QSZ`, `SPSU2`, `SPU1`, `QSZ4`, `SPZ4`, `P`, `NONE`, and `U1`.

The regression cases in `test/c++sym_conventions/` solve the same problem with
several symmetry types and compare all of them against one shared reference
(see `test/c++sym_conventions/README.md`). A new symmetry type that supports
pairing or spin-polarized baths should be added there.

## Why It Matters For Contributors

Many algorithmic changes are symmetry-agnostic, but any work touching:

- new model support
- new operator types
- invariant-subspace logic
- recalculation rules
- symmetry-dependent observables

will eventually cross this layer.

## Documentation Gap To Fill Later

This page is currently a framework overview. A good future expansion would document one concrete symmetry family end to end, for example `QS`, as a worked example of:

- quantum number encoding
- ancestor generation
- Hamiltonian construction
- operator recalculation
