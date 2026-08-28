# NRG Ljubljana

NRG Ljubljana is a hybrid Mathematica and C++ framework for numerical renormalization group calculations on quantum impurity models. The repository contains the low-level iterative solver, the model initialization pipeline, a large collection of analysis tools, and a substantial regression and unit test suite.

This documentation set provides an end-user path from installation through a
first calculation, references the supported configuration and current output
formats, and documents the architecture for contributors and advanced users.

## The Big Picture

The project has three major layers:

1. `nrginit/` prepares the initial Hamiltonian, basis, operators, and Wilson chain inputs.
2. `c++/` runs the iterative NRG calculation, truncation, operator recalculation, and optional DM-NRG or FDM workflows.
3. `tools/` provides standalone utilities for preprocessing and postprocessing data.

The input side is centered around the `param` and `data` files. The runtime side revolves around a parameter object, a symmetry backend, eigenspectra for each invariant subspace, operator containers, and iteration history stores.

## Start Here

- [Getting started](getting-started.md): build, test, and run the code locally
- [Parameter reference](parameter-reference.md): `param` syntax, ownership,
  defaults, constraints, and interactions
- [Output format reference](output-formats.md): current result filenames,
  columns, units, lifecycle, and implementation-dependent formats
- [Project structure](project-structure.md): top-level repository map
- [Input and configuration](input-and-configuration.md): how `param`, `data`, symmetry selection, and seed-state loading fit together
- [Runtime flow](runtime-flow.md): how execution moves from input to NRG phases and outputs
- [Iteration engine](iteration-engine.md): one iteration step in terms of the actual coordination functions in `core.hpp`
- [Data structures](data-structures.md): the central C++ types worth understanding first
- [State and persistence](state-and-persistence.md): what is stored on disk, what stays in memory, and how phases reuse state
- [nrginit workflow](nrginit-workflow.md): where Mathematica-side initialization fits into the overall pipeline
- [Developer guide](developer-guide.md): where to make changes for common tasks

## Documentation Scope

The goal of this MkDocs tree is to document both the user workflow and the
codebase as it exists in the repository:

- installation and the first end-to-end calculation
- user-facing parameters and current output layouts
- subsystem boundaries
- execution flow
- main data structures
- testing and maintenance workflow

It is not a full introduction to NRG theory. The older `doc/` tree remains
available for legacy scientific material that has not yet moved into this
site.
