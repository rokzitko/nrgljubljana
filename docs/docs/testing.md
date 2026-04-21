# Testing

The project uses a mix of unit tests, executable regression tests, tool tests, and Mathematica-dependent integration suites.

## Unit Tests

`test/unit/` contains focused C++ tests for:

- parameter parsing
- eigenspectrum containers
- operator storage
- workdir behavior
- subspace and symmetry helpers
- numerical utilities
- file I/O and HDF5 persistence

These tests are built as individual binaries and run through CTest.

## Regression Tests For `nrg`

The `test/c++/`, `test/c++sym_basic/`, `test/c++sym_more/`, `test/c++sym_all/`, and related directories exercise the full executable on prepared `param` and `data` inputs and compare outputs against reference results.

This is where end-to-end behavior across many symmetry types is checked.

## Tool Tests

`test/tools/` contains standalone tool coverage. These tests are useful when changing parsing, I/O, or command-line behavior in `tools/`.

## Mathematica-Dependent Suites

Several suites depend on Mathematica being available during configuration:

- `test/nrginit/`
- `test/nrginit_spsu2/`
- `test/nrginit+nrgrun/`
- `test/models/`
- `test/templates/`

These cover the higher-level initialization and generated-input workflow.

## Long Tests

Long-running suites are gated behind `-DTEST_LONG=ON`.

## Local Commands

Build and run the default test suite:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Useful focused runs from `CONTRIBUTING.md`:

```sh
ctest --test-dir build --output-on-failure -R '^(store|test_clean|test0_clean)$'
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)$'
ctest --test-dir build --output-on-failure -R '^(test_dmnrg_only|test_fdm_only|test65_algorithms_mats)$'
```

## CI Coverage

The current CI exercises:

- Linux and macOS build matrices
- sanitizer builds
- static analysis builds
- documentation builds

The CI setup is strict enough that it also catches test harness issues such as shared workdir races and platform-specific exception or filesystem behavior.
