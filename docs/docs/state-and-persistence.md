# State and Persistence

This page describes which data lives only inside the current phase, which data is archived across iterations, and which data is serialized to disk.

## Workdir-Owned State

`Workdir` in `c++/workdir.hpp` normally creates a unique temporary directory
for iteration artifacts. `nrg --checkpoint-dir DIR` instead creates or reopens
exactly `DIR`, requires `resume=true`, holds an exclusive process lock, and
does not recursively remove that caller-owned directory.

Important file families stored there include:

- unitary/eigenspectrum snapshots via `unitaryfn(...)`
- shell density matrices via `rhofn(...)`

With `resume=false` (the default), DM-NRG and FDM serialize and reload these
files between sweeps and check their structure on read. This same-run use
requires no input hashing, provenance scan, prior-chain discovery, or random IDs.

With `resume=true`, discovery runs after `data` establishes `Nmax` and accepts
only a contiguous unitary prefix from `Ninit`. Before opening result files,
the mandatory preflight reads each discovered unitary and checks that its
predecessor UUID matches the previous shell's generation UUID (nil predecessor
at `Ninit`). Preflight rejects a splice between independent chains even if
`param` and `data` are identical. Normal loads also check the generation IDs
against those expected for the shell. These comparisons use the existing
reads, not a separate pass to hash checkpoint contents. They do not detect
in-place payload edits that leave the IDs intact.

Cached shells are replayed without diagonalization to rebuild statistics,
stores, operators, and result files. If the complete unitary chain is reused,
each requested density family (`rho` or `rhofdm`) can be reused only if all
its required files are compatible and name the final unitary's UUID (the
chain tip). Otherwise that family is rebuilt. A partial restart computes the
missing shells with new generation IDs and rebuilds all requested densities.
Resumable archives remain in the checkpoint directory until the user removes them.

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

Workdir archives are written to temporary files and renamed after successful
serialization. Unitary files store the shell ground-state shift needed for
replay. See [Workdir files](output-formats.md#workdir-files) for the required
resumable headers and their compatibility limits.
Persistent result files such as `td`, spectra, and
HDF5 output are generally opened directly and can remain partial after an
interrupted run.

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

The current user-visible layouts and their compatibility limitations are
documented in the [output format reference](output-formats.md).

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
