# Getting Started

## Dependencies

The C++ code depends on:

- BLAS and LAPACK
- MPI
- OpenMP
- Boost
- GSL
- GMP
- HDF5

The `nrginit` side additionally depends on Mathematica. In normal workflows Mathematica is used to generate the `data` file consumed by the C++ executable.

## Build

Configure and build with CMake:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
cmake --install build
```

For a debug build:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DCMAKE_BUILD_TYPE=Debug
```

Useful configure options:

- `-DTEST_LONG=ON` enables long-running tests
- `-DASAN=ON -DUBSAN=ON` enables sanitizer builds
- `-DANALYZE_SOURCES=ON` turns on static analysis hooks
- `-DBuild_Documentation=ON -DSphinx_Only=ON` builds the legacy Sphinx docs
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

## Test

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure
```

Some tests depend on Mathematica and are only configured when Mathematica is available.

Useful focused runs from `CONTRIBUTING.md`:

```sh
ctest --test-dir build --output-on-failure -R '^(store|test_clean|test0_clean)$'
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)$'
ctest --test-dir build --output-on-failure -R '^(test_dmnrg_only|test_fdm_only|test65_algorithms_mats)$'
```

## Running The Code

The runtime executable is `nrg`. It expects `param` and `data` in the working directory.

Typical workflow:

1. generate `data` with the Mathematica-side initialization path
2. prepare `param`
3. run `nrg`
4. inspect the generated workdir outputs and tool postprocessing results

The executable also creates a temporary workdir for iteration artifacts and density-matrix files.

## Legacy Documentation

The older Sphinx documentation remains in `doc/` while this MkDocs tree is being expanded. Use this new tree for contributor-facing structure and code-orientation material.
