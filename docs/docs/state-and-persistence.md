# State and Persistence

This page describes which data lives only inside the current phase, which data is archived across iterations, and which data is serialized to disk.

## Workdir-Owned Temporary State

`Workdir` in `c++/workdir.hpp` creates a temporary directory for iteration artifacts.

Important file families stored there include:

- unitary/eigenspectrum snapshots via `unitaryfn(...)`
- shell density matrices via `rhofn(...)`

These files are used to bridge phases such as the second DM-NRG or FDM sweep.

## In-Memory Current-Step State

During a phase, the main mutable state is:

- current `DiagInfo<S>`
- current `Operators<S>`
- current `Stats<S>`
- phase-local `Output<S>`
- phase-local `Oprecalc<S>`

This is the state actively transformed by `iterate(...)` and `after_diag(...)`.

## Cross-Iteration Stores

Two long-lived stores are updated after each iteration:

### `ThermoStore<S>`

Defined in `c++/store.hpp`.

This stores compact eigenspectrum information needed for:

- thermodynamics
- FDM weights
- later backward-style processing

Each iteration contributes a `ThermoSubs<S>` object that stores compact `StoredEigen<S>` snapshots for all relevant invariant subspaces.

### `BackiterStore`

Also defined in `c++/store.hpp`.

This stores `BackiterSubs`, which capture:

- subspace dimensions
- numbers of kept states
- dimensions of eigenspectra

It is used by backward density-matrix propagation and operator-related postprocessing.

## Serialization Points

Several types serialize themselves directly:

- `DiagInfo<S>::save/load(...)` in `c++/eigen.hpp`
- `DensMatElements<S>::save/load(...)` in `c++/operators.hpp`
- HDF5 save helpers across `Values`, `DiagInfo`, `Operators`, `Stats`, and stores

The code consistently uses temporary files followed by rename for important binary outputs, which helps protect against partial writes.

## Output Front-End: `Output<S>`

`Output<S>` in `c++/output.hpp` owns the output sinks for one phase.

It may open and manage:

- `annotated.dat`
- `energies.nrg`
- `states.nrg`
- `report.nrg`
- `custom`
- `customfdm`
- `raw.h5` or `raw-dm.h5`

This makes `Output<S>` the phase-local front-end for persistence, while `ThermoStore` and `BackiterStore` are the cross-phase in-memory archive.

## Density-Matrix Persistence

The DM-NRG machinery in `c++/dmnrg.hpp` relies on saved per-shell data.

The broad pattern is:

1. initialize `rho` or `rhoFDM` at the last shell
2. save it to the workdir
3. backpropagate shell by shell using `BackiterStore` and stored eigenspectra
4. reload those files during the second sweep when needed

That is why the forward sweep must archive enough information even when the second phase is deferred.

## Persistence Boundaries In The Main Loop

Within `core.hpp`, the important persistence boundaries are:

- `persist_nrg_iteration_outputs(...)`
- `archive_iteration_state(...)`
- `diag.save(...)` for later DM/FDM use
- optional HDF5 dumps for structure, Hamiltonians, eigenspectra, and operators

## Practical Reading Order

If you are debugging file lifecycle or state reuse, read in this order:

1. `c++/workdir.hpp`
2. `c++/store.hpp`
3. `c++/output.hpp`
4. `c++/eigen.hpp`
5. `c++/operators.hpp`
6. `c++/dmnrg.hpp`
