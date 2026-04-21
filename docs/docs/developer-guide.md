# Developer Guide

## Local Workflow

The shortest normal edit cycle is:

1. configure a local build
2. build the affected targets
3. run focused tests first
4. run broader CTest coverage if the change touches shared code paths

Typical commands are collected in `CONTRIBUTING.md`.

## Where To Start Reading

For most C++ changes, start with:

- `c++/nrg-general.hpp`
- `c++/core.hpp`
- `c++/params.hpp`
- `c++/read-input.hpp`
- `c++/eigen.hpp`
- `c++/operators.hpp`
- `c++/store.hpp`

These files define most of the runtime vocabulary and the main orchestration flow.

Useful companion pages in this documentation set:

- [Input and configuration](input-and-configuration.md)
- [Runtime flow](runtime-flow.md)
- [Iteration engine](iteration-engine.md)
- [Data structures](data-structures.md)
- [State and persistence](state-and-persistence.md)

## Common Change Paths

### Runtime Flow Changes

Look at:

- `c++/nrg-general.hpp`
- `c++/core.hpp`
- `c++/step.hpp`

### Parameter Or Input Changes

Look at:

- `c++/params.hpp`
- `c++/read-input.hpp`
- `nrginit/`

### Symmetry Work

Look at:

- `c++/symmetry.hpp`
- `c++/mk_sym.hpp`
- `c++/sym-*.cc`
- `c++/symmetry/`

### Operator Recalculation Or Measurement Changes

Look at:

- `c++/operators.hpp`
- `c++/oprecalc.hpp`
- `c++/measurements.hpp`
- `c++/spectral.hpp`
- `c++/algo*.hpp`

### Tool Changes

Look at:

- `tools/CMakeLists.txt`
- the corresponding subdirectory under `tools/`
- the matching tests in `test/tools/`

## Testing Strategy

Prefer targeted tests first. The repository has enough breadth that broad runs are more useful after the local change is already stable.

If you touch:

- workdir or filesystem behavior: run unit tests and at least one `ctest -j` run
- `Params` or input parsing: run unit tests and one or more regression directories under `test/c++/`
- symmetry or recalculation logic: run both unit tests and regression suites
- tool code: run the matching tool tests under `test/tools/`

## Legacy And New Docs

The `doc/` tree is still present during migration. New contributor-oriented documentation should go into `docs/` unless there is a strong reason to extend the legacy Sphinx content instead.
