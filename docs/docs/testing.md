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
The `chain-fixed-seed` cases select a runtime reconstruction backend while
retaining the same committed seed; they do not test matched-backend seed
generation. Separate `chain-preparation` cases regenerate `data` before running
the solver and require Mathematica.

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

All four configuration conditions must be satisfied for these suites to be
registered: `Build_Tests=ON`, `SYM_ALL=ON`,
`NRGLJUBLJANA_ENABLE_MATHEMATICA=ON`, and successful Mathematica discovery.
The Conda recipe additionally defaults its
`nrgljubljana_enable_mathematica` and `nrgljubljana_build_tests` variants to
`OFF`; CMake kernel hints do not override those recipe variants.

## Long Tests

Long-running suites are gated behind `-DTEST_LONG=ON`.

## Independent Chain Backends

Scalar-chain producers use separate `base_legacy` and `base_rkpw` test names
and workdirs, sharing only immutable source fixtures and references. For
example, use `nrginit0_minimal_rkpw`, `test46_adapt_prepare_rkpw`, or the
dedicated `nrginit_pipeline_rkpw` (formerly `nrginit_rkpw_pipeline`).

`TEST_CHAIN_LEGACY` and `TEST_CHAIN_RKPW` both default to `ON`. Configure with
`-DTEST_CHAIN_LEGACY=OFF -DTEST_CHAIN_RKPW=ON` to retain only independent RKPW
coverage. A `chain-rkpw` test does not invoke the production legacy backend or
consume another test's generated outputs; deleting legacy registrations must
retain shared fixtures, references, and staging helpers. `chain-legacy` has the
corresponding independent meaning. Production selectors still default to
`tri=old` and `tridiag_method=lanczos`.

Optional `TEST_CHAIN_CROSSCHECK=ON` requires both backends enabled. These tests
carry `chain-crosscheck`, `chain-uses-legacy`, and `chain-uses-rkpw`, not the
independent backend labels. Neutral default-policy tests inspect parsing or
dispatch without executing a chain backend.

```sh
ctest --test-dir build -N -L '^chain-rkpw$'
ctest --test-dir build -L '^chain-rkpw$' --output-on-failure --no-tests=error
# Only after explicitly enabling crosschecks and both backends:
ctest --test-dir build -L '^chain-crosscheck$' --output-on-failure --no-tests=error
```

Use `ctest -N` to check coverage after applying the Mathematica, symmetry, long,
and scientific gates, rather than relying on a fixed total. The repository's
[test guide](https://github.com/rokzitko/nrgljubljana/blob/master/test/README.md#independent-chain-backends)
records the source-case inventory and the `ChainBackendTests.cmake` /
`chain-backend.pl` input-only staging contract. Backend tests do not regenerate
golden references.

## Scientific Validation

`TEST_SCIENTIFIC=ON` adds Python/NumPy ED checks. The 12 prepared-data NRG
comparisons and Python unit tests are license-free consumers, not generation
coverage. With Mathematica and `SYM_ALL`, three additional fixture-generation
tests per enabled backend run all four temperature/mode checks against the same
independent ED oracle. They use fresh candidates and `--check`, leaving
committed data and historical provenance unchanged.

```sh
# License-free scientific coverage, even if regeneration tests are registered:
ctest --test-dir build -L '^scientific$' -LE '^chain-generation$' \
  --output-on-failure --no-tests=error
# Independently regenerated RKPW coverage, when its gates are enabled:
ctest --test-dir build -R '^scientific_prepare_.*_rkpw$' \
  --output-on-failure --no-tests=error
```

See the [scientific suite guide](https://github.com/rokzitko/nrgljubljana/blob/master/test/scientific/README.md)
for Python requirements, model conventions, and validation limits.

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
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)'
ctest --test-dir build --output-on-failure -R '^(test_dmnrg_only|test_fdm_only|test65_algorithms_mats)$'
```

## CI Coverage

The current CI exercises:

- Linux and macOS build matrices
- sanitizer builds
- static analysis builds
- documentation builds

The CI setup is strict enough that it also catches test harness issues such as shared workdir races and platform-specific exception or filesystem behavior.
