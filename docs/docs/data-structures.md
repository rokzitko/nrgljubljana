# Data Structures

This page summarizes the C++ types that are most useful to understand before making nontrivial changes.

## How To Read This Page

The most useful mental split is:

- configuration and labels
- current-iteration state
- operator and measurement state
- cross-iteration stores
- phase-local helpers

The types below are grouped in that order.

## `Params`

Defined in `c++/params.hpp`.

`Params` is the central runtime configuration object. It owns:

- user-facing parameters parsed from `param`
- derived runtime values such as `Nmax` and `Nlen`
- optional `extra_params`
- the `Workdir`

This object influences nearly every major runtime decision: solver choice, truncation, algorithm selection, output behavior, and optional DM/FDM phases.

Typical ownership:

- owned by `NRG_calculation`
- referenced by nearly every subsystem
- partly mutable during startup because some values are derived from `data`

## `Invar` And `Twoinvar`

Defined in `c++/invar.hpp`.

`Invar` represents an invariant subspace label. It is the indexing language used throughout the codebase for:

- eigenspectra
- operator blocks
- subspace transitions
- recalculation routines

`Twoinvar` is the paired index used for block-structured matrices such as operator matrix elements between two invariant subspaces.

Why it matters:

- it is the common index space for block matrices
- many maps in `operators.hpp` and recalculation code are keyed by `Invar` or `Twoinvar`

## `Symmetry<S>`

Defined in `c++/symmetry.hpp`.

`Symmetry<S>` is the core abstraction that tells the solver how to:

- generate allowed new subspaces
- enumerate ancestor relations
- check coupling rules
- build Hamiltonian blocks
- recalculate operators
- provide multiplicity and spectral prefactors

Concrete symmetry implementations are selected through `mk_sym.hpp` and compiled from the `sym-*.cc` sources.

Typical ownership:

- constructed once during input loading
- shared across the runtime through `std::shared_ptr<Symmetry<S>>`
- consulted in nearly every iteration step

## `DiagInfo<S>` And `Eigen<S>`

Defined in `c++/eigen.hpp`.

These classes represent the eigenspectrum at one iteration.

- `DiagInfo<S>` is the map from `Invar` to `Eigen<S>`
- `Eigen<S>` stores one block's eigenvalues, eigenvectors, truncation metadata, and block-separated eigenvectors used for recalculation

Supporting types in the same file include:

- `Values<S>`: raw and derived eigenvalue views
- `Vectors<S>`: dense eigenvector storage
- `Blocks<S>`: eigenvectors split according to ancestor subspaces
- `StoredEigen<S>`: compact representation for archived thermodynamic storage

Lifecycle notes:

- `RawEigen<S>` is a short-lived LAPACK-facing result object
- `Eigen<S>` is the runtime object used during the current step
- `StoredEigen<S>` is the compact archived form kept in stores

## `Operators<S>` And `MatrixElements<S>`

Defined in `c++/operators.hpp`.

These structures store physical operators and density matrices in block form.

- `MatrixElements<S>` maps `Twoinvar` to a matrix block
- `Opch<S>` stores channel-resolved `f` operators
- `Operators<S>` groups all operator families used by the runtime
- `DensMatElements<S>` stores density-matrix blocks

The main operator categories in `Operators<S>` include singlet, parity-odd singlet, global singlet, doublet, triplet, quadruplet, and orbital-triplet operators.

Lifecycle notes:

- `InputData` seeds `Operators<S>` from the `data` file
- `Oprecalc<S>` determines which operator families actually need recalculation
- `core.hpp` trims operator matrices after each step to reduce memory footprint

## `SubspaceDimensions`, `SubspaceStructure`, And `TaskList`

Defined in `c++/subspaces.hpp`.

These types are the bridge between one iteration and the next.

- `SubspaceDimensions` describes how one new subspace is assembled from ancestor subspaces and what each ancestor contributes
- `SubspaceStructure` maps every new `Invar` to its `SubspaceDimensions`
- `TaskList` turns the structure into an ordered list of diagonalization jobs

If you want to understand Hamiltonian construction and block splitting, this is one of the first files to read.

In practice these types form the handoff boundary between:

- the previous iteration's eigenspectra
- the next set of diagonalization tasks
- the block structure used in recalculation

## `Coef<S>`

Defined in `c++/coef.hpp`.

`Coef<S>` contains the Wilson-chain and related discretization coefficients loaded from the `data` file. It is the numerical bridge between the initialization stage and the iterative Hamiltonian construction.

## `Stats<S>`

Defined in `c++/stats.hpp`.

`Stats<S>` tracks thermodynamic and per-phase aggregate quantities such as partition-function-like values, energies, and reporting data that are accumulated as the run proceeds.

It also owns the in-memory values later written into `td`, `tdfdm`, and optional HDF5 outputs.

## `ThermoStore<S>` And `BackiterStore`

Defined in `c++/store.hpp`.

These structures archive information across iterations.

- `ThermoStore<S>` keeps compact eigenspectrum information needed for thermodynamic and FDM calculations
- `BackiterStore` keeps per-iteration subspace metadata needed by backward-style passes

These stores become especially important once the runtime moves beyond the plain forward NRG sweep.

## Phase-Local Coordination Objects

Some important types are not long-lived stores but are central to the flow:

- `InputData<S>` in `c++/read-input.hpp`: reads the seed eigenspectrum, operators, and chain coefficients from `data`
- `Output<S>` in `c++/output.hpp`: owns text and HDF5 output sinks for the current phase
- `Oprecalc<S>` in `c++/oprecalc.hpp`: decides which operators and spectra need work and manages spectral algorithm instances
- `Step` in `c++/step.hpp`: packages iteration index, energy scale, runtime mode, and convenience helpers used throughout the loop

## Ownership Sketch

At runtime, the main ownership picture looks roughly like this:

```text
NRG_calculation
  ├── Params
  ├── InputData
  │    ├── Symmetry<S>
  │    ├── DiagInfo<S> (seed)
  │    ├── Operators<S> (seed)
  │    └── Coef<S>
  ├── Stats<S>
  ├── ThermoStore<S>
  └── BackiterStore

Per phase:
  run_phase(...)
    ├── Oprecalc<S>
    ├── Output<S>
    ├── Step
    └── current DiagInfo<S> / Operators<S>
```

## How They Fit Together

At a high level:

1. `Params` and `InputData` set up the run.
2. `Symmetry<S>` defines how invariant subspaces combine.
3. `DiagInfo<S>` stores the current eigenspectra.
4. `SubspaceStructure` and `TaskList` define the next diagonalization work.
5. `Operators<S>` are recalculated in the new basis.
6. `Stats<S>`, `ThermoStore<S>`, and `BackiterStore` preserve information needed for reporting and later phases.

If you want the function-level sequencing that moves these objects around, see [Iteration engine](iteration-engine.md).
