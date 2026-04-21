# Runtime Flow

## Inputs

The runtime side centers around two input files in the working directory:

- `param`: runtime parameters and workflow switches
- `data`: generated initial eigenspectrum, operators, and chain coefficients

`param` is parsed by `Params` in `c++/params.hpp`.

`data` is parsed by `InputData` in `c++/read-input.hpp`.

## Entry Point

Execution starts in `c++/nrg.cc`:

1. initialize MPI
2. report OpenMP and environment information
3. validate that `param` and `data` exist
4. set up the temporary workdir
5. run master or slave behavior depending on MPI rank

The executable-facing orchestration then moves into `NRG_calculation` in `c++/nrg-general.hpp`.

`nrg.cc` is intentionally thin. Nearly all interesting control flow begins once `NRG_calculation` is constructed.

## Initialization In `NRG_calculation`

`NRG_calculation` constructs and wires together the core runtime objects:

1. `Params`
2. `InputData`
3. `Symmetry<S>`
4. `Stats<S>`
5. `ThermoStore<S>`
6. `BackiterStore`
7. diagonalization engine selection

Important methods in `NRG_calculation`:

- `select_diag_engine()`: choose MPI / OpenMP / serial diagonalization backend
- `run_phase(...)`: build per-phase helpers and execute the iteration loop
- `prepare_rho(...)`: initialize and backpropagate shell density matrices
- `prepare_rhoFDM()`: initialize and backpropagate FDM density matrices
- `run_dm_phase(...)`: launch the second sweep when enabled

`InputData` performs these steps:

1. read the `data` header
2. instantiate the symmetry backend
3. load the seed eigenspectrum and seed operators
4. read operator blocks and coefficient tables
5. finalize `Nmax` and related derived parameters

See [Input and configuration](input-and-configuration.md) for the input-side details.

## Main NRG Phase

The main run is started through `run_phase(RUNTYPE::NRG, ...)`, which enters `nrg_loop(...)` in `c++/core.hpp`.

At a high level, each iteration does the following:

1. build the next-shell `SubspaceStructure`
2. convert that into a `TaskList`
3. construct the Hamiltonian blocks for each invariant subspace
4. diagonalize each block through the selected `DiagEngine`
5. establish energy references and truncation criteria
6. split eigenvectors into ancestor blocks when needed
7. update statistics and outputs
8. recalculate operators in the new basis
9. compute observables and spectral data for the current phase
10. archive state into thermodynamic and backward-sweep stores
11. truncate the eigenspectra
12. recalculate irreducible operators for the next step

The most important coordination logic lives in `c++/core.hpp`.

See [Iteration engine](iteration-engine.md) for the function-level breakdown of this phase.

## Optional DM-NRG And FDM Phases

After the first sweep, `NRG_calculation` may run density-matrix-related phases depending on `Params`:

- prepare `rho`
- prepare `rhoFDM`
- run `RUNTYPE::DMNRG`

This is still driven through the same broad iteration machinery, but with different runtime mode and output semantics.

## Outputs

The runtime produces several forms of output:

- spectral and thermodynamic files in the working directory
- temporary stored eigenspectra and density matrices inside the generated workdir
- optional HDF5 output via `Output` and `h5save` helpers
- logs and diagnostic dumps when enabled by parameters

The workdir is represented by `Workdir` in `c++/workdir.hpp` and is used for transient iteration-state files such as unitary matrices and density matrices.

See [State and persistence](state-and-persistence.md) for the storage model and serialization boundaries.

## Mathematica Side

The C++ executable does not construct the initial problem description by itself. The Mathematica-side `nrginit` pipeline is responsible for:

- model definition
- basis construction
- initial Hamiltonian diagonalization
- seed operator generation
- `data` file emission

That is why understanding the full project flow requires keeping both `nrginit/` and `c++/` in mind, even though the iterative heavy lifting happens in C++.

See [nrginit workflow](nrginit-workflow.md) for the Mathematica-side entry path.
