[![build](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml)
[![conda](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml)

# NRG Ljubljana

NRG Ljubljana is a framework for numerical renormalization group (NRG) calculations for quantum impurity problems. It combines a Mathematica-based initialization layer 
(using [sneg](https://github.com/rokzitko/sneg)) with a C++20 runtime for the iterative NRG calculation, density-matrix variants, and a collection of standalone analysis tools.

NRG Ljubljana is used in published research on quantum dots, Kondo systems,
superconducting impurities, molecular magnets, DMFT, and correlated materials.
See [publications using NRG Ljubljana](docs/docs/publications.md).

If you publish work using NRG Ljubljana, please add your paper via a
[pull request](https://github.com/rokzitko/nrgljubljana/edit/master/docs/docs/publications.md)
or [open an issue](https://github.com/rokzitko/nrgljubljana/issues).

## What It Covers

- multiple symmetry backends, including QS, QSZ, ISO, ISOSZ, SPSU2, and extended symmetry sets
- standard NRG iteration together with CFS, DM-NRG, and FDM workflows
- thermodynamic quantities, expectation values, spectral functions, Matsubara quantities, and conductance-related calculations
- preprocessing and postprocessing tools for discretization, chain generation, broadening, Hilbert transforms, resampling, and file conversion
- structured output in text, binary, and HDF5 formats

## Applications

- dynamical mean-field theory (DMFT); see a simple sample code for [Kondo lattice
model](https://github.com/rokzitko/DMFT_NRG_KLM)

## Installation

NRG Ljubljana is available from conda-forge:

   conda install -c conda-forge nrgljubljana

Packages are available for Linux x86-64, Linux aarch64, macOS x86-64,
and macOS arm64. Source builds remain appropriate when CUDA support or
nonstandard numerical-library configurations are required.

## Source builds

### Requirements

Source builds require CMake 3.25 or newer and C and C++ compilers with
C++20 support. The following dependencies must be provided by the system in
all dependency modes:

- threaded BLAS and LAPACK, typically MKL or OpenBLAS
- MPI with C++ support
- Boost 1.68 or newer, including the MPI and serialization components
- GSL 2.0 or newer
- GMP, including the C++ library
- HDF5 with the C and high-level components

OpenMP is required only when application-level OpenMP regions are enabled or
when the selected BLAS/LAPACK threading implementation requires it. Wolfram
Mathematica is required for the `nrginit` workflow, which prepares the initial
Hamiltonian, basis, and operator data used by the C++ executable. CUDA support
is optional and requires a CUDA Toolkit that provides cuSOLVER, cuBLAS, and
cuDART.

By default, `NRGLJUBLJANA_USE_SYSTEM_DEPS` is `OFF`. Configuration downloads
CPM.cmake and uses it to obtain GoogleTest 1.15.2 when tests are enabled,
Eigen 3.4.0 when a compatible system Eigen package is unavailable, HighFive
3.1.0, fmt 12.2.0, and range-v3 0.12.0. The first configure therefore requires
Git and network access unless these sources are already present in the CPM
source cache.

For a network-free build, install CMake config packages for Eigen 3.4 or newer,
HighFive, fmt 11 or newer, and range-v3, then configure with
`-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON`. GoogleTest must also be installed when the
default test build is enabled. The test suite additionally uses Perl and the
HDF5 command-line tools; use `-DBuild_Tests=OFF` when tests are not needed.

### Build and test

CMake requires an explicit absolute install prefix. Configure, build, install,
and load the installed environment with:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
cmake --install build
source "$HOME/nrgljubljana/share/nrgljubljanavars.sh"
```

Source `nrgljubljanavars.sh` in each new shell that will use the installation.
For a nondefault `CMAKE_INSTALL_DATADIR`, the script is under that directory
instead of `share`; CMake prints its exact path during configuration. The
script configures the executable, library, include, and CMake package search
paths.

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure --timeout 3600 --no-tests=error
```

Increase `--timeout` for slow machines or debug builds.

Useful developer options:

- `-DCMAKE_BUILD_TYPE=Debug`
- `-DBuild_Tests=ON|OFF` controls the test build (default: `ON` for a top-level build)
- `-DTEST_LONG=ON`
- `-DASAN=ON -DUBSAN=ON`
- `-DANALYZE_SOURCES=ON`
- `-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON|OFF` enables application-level OpenMP regions such as simultaneous diagonalisation scheduling (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_CUDA=ON|OFF` requests CUDA/cuSOLVER support (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_BLAS_ILP64=ON|OFF` selects the 64-bit or 32-bit BLAS/LAPACK integer ABI (default: `OFF`)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

### BLAS/LAPACK integer ABI

The default LP64 BLAS/LAPACK interface uses 32-bit integers even on a 64-bit
system. `NRGLJUBLJANA_BLAS_ILP64=ON` changes the project declarations and
library search to the 64-bit-integer ILP64 interface. BLAS, LAPACK, and the
project must use the same ABI; mixing LP64 and ILP64 can link but is unsafe.
Use a matching implementation, such as `BLA_VENDOR=Intel10_64lp` for MKL LP64
or `BLA_VENDOR=Intel10_64ilp` for MKL ILP64. An ILP64 OpenBLAS build must also
export a compatible symbol convention. Apple Accelerate and the standard
Conda variants are LP64. See the
[getting-started guide](docs/docs/getting-started.md#blaslapack-integer-abi) for
the full compatibility constraints.

## Parallelism Model

NRG Ljubljana's default performance model is BLAS/LAPACK-internal threading. The executable should normally have one numerical threading backend in the process: threaded MKL or threaded OpenBLAS. Application-level OpenMP regions are disabled by default so the code does not accidentally link a second OpenMP runtime such as GNU `libgomp` together with Intel `libiomp5`.

Use `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and scheduler CPU binding to control numerical kernel threading. When running with MPI, choose the BLAS/LAPACK thread count together with the rank count; `mpi_ranks * blas_threads` should usually not exceed the CPUs allocated to the job.

For MKL builds that use the `mkl_rt` dispatcher, set `-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU`, `INTEL`, or `LLVM` when you need an explicit threading backend. This links the matching compiler OpenMP runtime through CMake's `OpenMP::OpenMP_CXX` target while keeping application-level OpenMP regions disabled unless `NRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is also set.

`-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is an expert option for simultaneous diagonalisation scheduling (`diag_mode=OpenMP`, `diagth>1`) and a few non-BLAS loops. It can create nested parallelism when BLAS/LAPACK is also threaded, so CMake checks the visible link line for mixed OpenMP runtime families and the executable prints startup diagnostics and warnings about the detected MKL/OpenBLAS/OpenMP/MPI threading configuration.

## Repository Map

- `c++/`: core NRG engine, runtime flow, diagonalization, symmetry framework, operators, stores, and numerical utilities
- `tools/`: standalone preprocessing and postprocessing executables
- `nrginit/`: Mathematica-side model initialization and input generation
- `nrgspawn/`: examples and reference outputs for prepared runs generated by `tools/instantiate/`
- `test/`: unit tests, regression suites, tool tests, and Mathematica-driven integration tests
- `share/`: installed auxiliary CMake files and runtime assets
- `scripts/`: small helper scripts for inspecting and postprocessing outputs
- `doc/`: legacy Sphinx documentation
- `docs/`: new MkDocs documentation tree

## Documentation

The in-tree documentation refresh is being migrated to MkDocs under `docs/`, while the legacy Sphinx content in `doc/` remains available.

MkDocs documentation: http://auger.ijs.si/nrgljubljana/site/

## Contributing

See `CONTRIBUTING.md` for the local build, test, sanitizer, analysis, and documentation commands used in development.

## License

NRG Ljubljana is distributed under the GNU General Public License. See `COPYING` for the full license text.

## Contact

- project home page: http://nrgljubljana.ijs.si/
- Rok Zitko, "Jozef Stefan" Institute, Ljubljana, Slovenia
- rok.zitko@ijs.si

## Acknowledgements

NRG Ljubljana started during Rok Zitko's PhD work at the University of Ljubljana and the "Jozef Stefan" Institute. The codebase reflects collaboration and discussions with multiple researchers in the NRG community and contributions from collaborators over many years.
