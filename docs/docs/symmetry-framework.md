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

The actual implementations are compiled from the `sym-*.cc` sources in `c++/` and the associated recalculation helpers in `c++/symmetry/`.

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

- For the abelian types `NONE`, `P`, `U1`, `QSZ`, and `SPU1`, every combination is a
  single Fock state of the added sites. `AbelianShell::Shell` in `c++/abelian_shell.hpp`
  enumerates these states and computes all coefficients as fermion matrix elements,
  for any number of channels. `test/unit/symmetry/abelian_shell.cpp` checks that the
  result equals the Mathematica tables, which it reads from
  `c++/symmetry/coefnew/{none,p,u1,qsz,spu1}`, entry by entry.
- All other types use tables generated with Mathematica (`c++/symmetry/coefnew`,
  driven by `nrgcoef.m`) for each supported number of channels, included through
  the m4 sources `c++/symmetry/nrg-recalc-*.hpp.m4`.

The conventions of the abelian tables are those of `nrgcoef.m`: modes are ordered by
channel with spin up before spin down; a Fock state is the product of the creation
operators of the occupied modes in increasing mode order; and in a combination
`|i, r>` of an added-site state `|i>` with a state `|r>` of the previous sites, the
added-site operators stand to the left. An operator `O` of the previous sites then
has `<i, r|O|j, r'> = delta_ij (-1)^(n_i |O|) <r|O|r'>`, with `n_i` the number of
electrons in `|i>` and `|O|` the fermion parity of `O`.

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
  accepted only by `SPSU2`, `SPU1`, `SPU1LR`, `SPSU2LR`, `SPSU2T`, `SPSU2C3`, `P`,
  `PP`, and `NONE`; nonzero anomalous-hopping (`kappa`) coefficients only by
  `SPSU2`, `SPU1`, `SPU1LR`, `SPSU2LR`, and `SPSU2T`. `substeps=true` is accepted
  only by `QS`, `QSZ`, `SPSU2`, `SPU1`, `P`, `NONE`, and `U1`.

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
