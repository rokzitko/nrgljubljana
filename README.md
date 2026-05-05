[![build](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml)
[![conda](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml)

# NRG Ljubljana

NRG Ljubljana is a framework for numerical renormalization group (NRG) calculations for quantum impurity problems. It combines a Mathematica-based initialization layer with a C++20 runtime for the iterative NRG calculation, density-matrix variants, and a collection of standalone analysis tools.

## What It Covers

- multiple symmetry backends, including QS, QSZ, ISO, ISOSZ, SPSU2, and extended symmetry sets
- standard NRG iteration together with CFS, DM-NRG, and FDM workflows
- thermodynamic quantities, expectation values, spectral functions, Matsubara quantities, and conductance-related calculations
- preprocessing and postprocessing tools for discretization, chain generation, broadening, Hilbert transforms, resampling, and file conversion
- structured output in text, binary, and HDF5 formats

## Quick Start

Build locally with CMake:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
cmake --install build
```

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure --timeout 3600 --no-tests=error
```

Increase `--timeout` for slow machines or debug builds.

Useful developer options:

- `-DCMAKE_BUILD_TYPE=Debug`
- `-DTEST_LONG=ON`
- `-DASAN=ON -DUBSAN=ON`
- `-DANALYZE_SOURCES=ON`
- `-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON|OFF` enables application-level OpenMP regions such as simultaneous diagonalisation scheduling (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

Core native dependencies:

- threaded BLAS and LAPACK, typically MKL or OpenBLAS
- MPI
- OpenMP only when enabled for application-level regions, or when required by the selected BLAS/LAPACK implementation
- Boost
- GSL
- GMP
- HDF5

Wolfram Mathematica is required for the `nrginit` side of the workflow, which prepares the initial Hamiltonian, basis, and operator data used by the C++ executable.

## Parallelism Model

NRG Ljubljana's default performance model is BLAS/LAPACK-internal threading. The executable should normally have one numerical threading backend in the process: threaded MKL or threaded OpenBLAS. Application-level OpenMP regions are disabled by default so the code does not accidentally link a second OpenMP runtime such as GNU `libgomp` together with Intel `libiomp5`.

Use `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and scheduler CPU binding to control numerical kernel threading. When running with MPI, choose the BLAS/LAPACK thread count together with the rank count; `mpi_ranks * blas_threads` should usually not exceed the CPUs allocated to the job.

For MKL builds that use the `mkl_rt` dispatcher, set `-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU`, `INTEL`, or `LLVM` when you need an explicit threading backend. This links the matching compiler OpenMP runtime through CMake's `OpenMP::OpenMP_CXX` target while keeping application-level OpenMP regions disabled unless `NRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is also set.

`-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is an expert option for simultaneous diagonalisation scheduling (`diag_mode=OpenMP`, `diagth>1`) and a few non-BLAS loops. It can create nested parallelism when BLAS/LAPACK is also threaded, so CMake checks the visible link line for mixed OpenMP runtime families and the executable prints startup diagnostics and warnings about the detected MKL/OpenBLAS/OpenMP/MPI threading configuration.

## Repository Map

- `c++/`: core NRG engine, runtime flow, diagonalization, symmetry framework, operators, stores, and numerical utilities
- `tools/`: standalone preprocessing and postprocessing executables
- `nrginit/`: Mathematica-side model initialization and input generation
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
