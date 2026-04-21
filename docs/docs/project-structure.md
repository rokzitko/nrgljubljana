# Project Structure

## Top-Level Layout

- `c++/`: core solver, runtime orchestration, diagonalization wrappers, operator recalculation, stores, output handling, and symmetry abstraction
- `tools/`: standalone preprocessing and postprocessing executables such as `adapt`, `nrgchain`, `broaden`, `hilb`, `resample`, and `unitary`
- `nrginit/`: Mathematica-side initialization code used to generate the initial `data` file and related model artifacts
- `test/`: unit tests, regression tests, tool tests, and Mathematica-driven integration suites
- `share/`: installed package metadata and auxiliary CMake modules used by the build and downstream consumers
- `scripts/`: helper scripts for inspecting or aggregating output files
- `doc/`: legacy Sphinx documentation tree
- `docs/`: new MkDocs documentation tree
- `templates/`: model and workflow templates used by higher-level generation or tests
- `floquet/`: repository-local example and experiment directories for Floquet-related workflows

## Core Native Code In `c++/`

The `c++/` directory is the heart of the runtime. The most important entry and coordination files are:

- `nrg.cc`: executable entry point
- `nrg-general.hpp`: high-level orchestration through `NRG_calculation`
- `core.hpp`: main iteration logic and postprocessing helpers
- `read-input.hpp`: parses the generated `data` file into in-memory structures
- `params.hpp`: runtime parameter model and parsing of `param`
- `symmetry.hpp`: abstract symmetry interface
- `mk_sym.hpp`: symmetry backend selection
- `diag.hpp`: LAPACK wrappers and eigensolver selection
- `eigen.hpp`: eigenspectrum storage and serialization
- `operators.hpp`: operator matrix containers
- `subspaces.hpp`: subspace layout and task generation
- `store.hpp`: thermodynamic and backward-sweep storage

If you are trying to orient yourself quickly, the most useful reading order is:

1. `nrg.cc`
2. `nrg-general.hpp`
3. `read-input.hpp`
4. `core.hpp`
5. `eigen.hpp`, `operators.hpp`, `store.hpp`, `subspaces.hpp`
6. `symmetry.hpp` and `mk_sym.hpp`

Concrete symmetry backends are compiled from many `sym-*.cc` files and corresponding headers.

See also:

- [Input and configuration](input-and-configuration.md)
- [Iteration engine](iteration-engine.md)
- [Data structures](data-structures.md)

## Standalone Tools In `tools/`

Each tool is built as its own executable and linked against the project's common utility code and core numerical dependencies. The tool set covers several categories:

- discretization and chain generation: `adapt`, `nrgchain`
- matrix and numerical transforms: `diag`, `hilb`, `kk`, `unitary`
- spectral postprocessing: `broaden`, `bw`, `resample`, `specmoments`
- averages and integration helpers: `binavg`, `intavg`, `integ`, `tdavg`, `mats`
- file and format utilities: `h5write`, `matrix`

## Test Layout In `test/`

The test tree is organized by style rather than by component:

- `test/unit/`: focused C++ unit tests for core data structures and utilities
- `test/c++/`: main regression runs for the C++ executable
- `test/c++sym_basic/`, `test/c++sym_more/`, `test/c++sym_all/`: regression suites across increasing symmetry coverage
- `test/tools/`: standalone tool tests
- `test/nrginit*`, `test/models*`, `test/templates/`: Mathematica-dependent end-to-end workflows
- `test/complex/`, `test/test_long/`: specialized or long-running suites

## Documentation During Migration

For now there are two documentation trees:

- `doc/` is the existing Sphinx content
- `docs/` is the new MkDocs content

The intent is to move contributor-oriented documentation into `docs/` first, while keeping the existing material available during the transition.
