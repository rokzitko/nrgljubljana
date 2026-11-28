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
  only by `QS`, `QSZ`, `SPSU2`, and `SPU1`.

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
