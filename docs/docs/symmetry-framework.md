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
