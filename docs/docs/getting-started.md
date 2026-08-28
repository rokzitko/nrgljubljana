# Getting Started

## Dependencies

Source builds require CMake 3.25 or newer and C and C++ compilers with
C++20 support. These dependencies must be installed on the system regardless
of how the remaining C++ dependencies are resolved:

- threaded BLAS and LAPACK, typically MKL or OpenBLAS
- MPI with C++ support
- Boost 1.53 or newer, including the MPI and serialization components
- GSL
- GMP, including the C++ library
- HDF5 with the C and high-level components

OpenMP is required only when application-level OpenMP regions are enabled or
when the selected BLAS/LAPACK threading implementation requires it. It is not
an unconditional dependency. The `nrginit` workflow additionally requires
Wolfram Mathematica to generate the `data` file consumed by the C++ executable.
Optional CUDA support requires a CUDA Toolkit that provides cuSOLVER, cuBLAS,
and cuDART.

### Dependency resolution

The default `NRGLJUBLJANA_USE_SYSTEM_DEPS=OFF` configuration downloads
CPM.cmake and uses it to obtain:

- GoogleTest 1.15.2 when tests are enabled (the default for a top-level build)
- Eigen 3.4.0 when a compatible system Eigen package is unavailable
- HighFive 3.1.0
- fmt 12.2.0
- range-v3 0.12.0

The first configure requires Git and network access unless these sources are
already available in the CPM source cache.

For a network-free build, install CMake config packages for Eigen 3.4 or newer,
HighFive, fmt 11 or newer, and range-v3, then configure with
`-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON`. GoogleTest must also be installed when
tests are enabled. The test suite additionally uses Perl and the HDF5
command-line tools; configure with `-DBuild_Tests=OFF` when tests are not
needed.

## Build

CMake requires an explicit absolute install prefix. Configure, build, install,
and load the installed environment with:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
cmake --install build
source "$HOME/nrgljubljana/share/nrgljubljanavars.sh"
```

Source `nrgljubljanavars.sh` in each new shell that will use the installation.
It sets the executable, library, include, and CMake package search paths. With
a nondefault `CMAKE_INSTALL_DATADIR`, the script is installed under that
directory instead of `share`; CMake prints the exact path during configuration.

For a debug build:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DCMAKE_BUILD_TYPE=Debug
```

Useful configure options:

- `-DBuild_Tests=ON|OFF` controls the test build (default: `ON` for a top-level build)
- `-DTEST_LONG=ON` enables long-running tests
- `-DASAN=ON -DUBSAN=ON` enables sanitizer builds
- `-DANALYZE_SOURCES=ON` turns on static analysis hooks
- `-DBuild_Documentation=ON -DSphinx_Only=ON` builds the legacy Sphinx docs
- `-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON|OFF` enables application-level OpenMP regions (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_CUDA=ON|OFF` requests CUDA/cuSOLVER support (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

## Test

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure --timeout 3600 --no-tests=error
```

Increase `--timeout` for slow machines or debug builds.

Some tests depend on Mathematica and are only configured when Mathematica is available.

Useful focused runs from `CONTRIBUTING.md`:

```sh
ctest --test-dir build --output-on-failure -R '^(store|test_clean|test0_clean)$'
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)'
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
