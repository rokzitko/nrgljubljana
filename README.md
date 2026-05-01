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
ctest --test-dir build --output-on-failure
```

Useful developer options:

- `-DCMAKE_BUILD_TYPE=Debug`
- `-DTEST_LONG=ON`
- `-DASAN=ON -DUBSAN=ON`
- `-DANALYZE_SOURCES=ON`
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

Core native dependencies:

- BLAS and LAPACK
- MPI
- OpenMP
- Boost
- GSL
- GMP
- HDF5

Wolfram Mathematica is required for the `nrginit` side of the workflow, which prepares the initial Hamiltonian, basis, and operator data used by the C++ executable.

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
