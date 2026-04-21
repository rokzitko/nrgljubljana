# NRG Ljubljana

NRG Ljubljana is a hybrid Mathematica and C++ framework for numerical renormalization group calculations on quantum impurity models. The repository contains the low-level iterative solver, the model initialization pipeline, a large collection of analysis tools, and a substantial regression and unit test suite.

This documentation set is aimed primarily at contributors and advanced users who need to understand how the repository is organized and how data flows through the code.

## The Big Picture

The project has three major layers:

1. `nrginit/` prepares the initial Hamiltonian, basis, operators, and Wilson chain inputs.
2. `c++/` runs the iterative NRG calculation, truncation, operator recalculation, and optional DM-NRG or FDM workflows.
3. `tools/` provides standalone utilities for preprocessing and postprocessing data.

The input side is centered around the `param` and `data` files. The runtime side revolves around a parameter object, a symmetry backend, eigenspectra for each invariant subspace, operator containers, and iteration history stores.

## Start Here

- [Getting started](getting-started.md): build, test, and run the code locally
- [Project structure](project-structure.md): top-level repository map
- [Input and configuration](input-and-configuration.md): how `param`, `data`, symmetry selection, and seed-state loading fit together
- [Runtime flow](runtime-flow.md): how execution moves from input to NRG phases and outputs
- [Iteration engine](iteration-engine.md): one iteration step in terms of the actual coordination functions in `core.hpp`
- [Data structures](data-structures.md): the central C++ types worth understanding first
- [State and persistence](state-and-persistence.md): what is stored on disk, what stays in memory, and how phases reuse state
- [nrginit workflow](nrginit-workflow.md): where Mathematica-side initialization fits into the overall pipeline
- [Developer guide](developer-guide.md): where to make changes for common tasks

## Documentation Scope

The goal of this MkDocs tree is to document the codebase as it exists in the repository:

- subsystem boundaries
- execution flow
- main data structures
- testing and maintenance workflow

It is not yet intended as a full scientific manual or end-user tutorial replacement for all legacy materials in `doc/`.
