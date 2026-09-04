[![build](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/build.yml)
[![conda](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml/badge.svg)](https://github.com/rokzitko/nrgljubljana/actions/workflows/conda.yml)

# NRG Ljubljana

NRG Ljubljana implements Wilson's numerical renormalization group (NRG) for
equilibrium quantum impurity models. The continuum bath is logarithmically
discretized and mapped to a Wilson chain, whose sites are added successively to
the impurity. The Hamiltonian is diagonalized in symmetry sectors after each
addition and high-energy states are discarded. The resulting flow of rescaled
many-body levels and the associated thermodynamic and dynamical observables
span many decades of energy and temperature. For the single-impurity Anderson
model, this resolves the free-orbital and local-moment regimes, Kondo screening,
and the strong-coupling fixed point.

The impurity Hamiltonian, symmetry-adapted basis, and local operators are
constructed in Mathematica using [sneg](https://github.com/rokzitko/sneg). The
numerically intensive iterative diagonalization and density-matrix calculations
are performed in C++20. Additional programs carry out logarithmic bath
discretization, Wilson-chain construction, Hilbert transforms, and spectral
broadening.

## Physical Quantities and NRG Schemes

- entropy, susceptibility, heat capacity, free energy, and local expectation values
- many-body level flows, crossover scales, and fixed-point spectra
- real-frequency spectral densities and Matsubara-frequency correlation functions
- linear-response conductance, transport moments, and finite-temperature response functions
- conventional NRG, complete-Fock-space (CFS), density-matrix NRG (DM-NRG), and full-density-matrix NRG (FDM-NRG)
- Abelian and non-Abelian symmetry sectors for charge, spin, axial charge (isospin), and orbital degrees of freedom
- logarithmic discretization, Wilson-chain coefficients, Hilbert transforms, $z$-averaging, and spectral broadening

## Installation

Install NRG Ljubljana from conda-forge in a dedicated environment with one
command:

```sh
conda create -n nrg -c conda-forge nrgljubljana
```

Activate it with `conda activate nrg`.

Packages are available for Linux x86-64, Linux aarch64, macOS x86-64,
and macOS arm64. See the [full Conda installation guide](README.conda-forge.md)
for Miniforge setup, updating, and troubleshooting.

The conda-forge distribution provides all programs needed for bath
discretization, Wilson-chain construction, iterative diagonalization, and
spectral analysis.
Mathematica is needed only to derive the symmetry-adapted Hamiltonian and
operator matrix elements for a new impurity model with `nrginit`. The tutorial
calculations below contain these matrix elements already and do not require
Mathematica.

## First Calculation

The [NRG Ljubljana SIAM tutorial](https://github.com/rokzitko/NRG_SIAM) is the
recommended introduction. Its first calculation considers the
particle-hole-symmetric single-impurity Anderson model with bath half-bandwidth
$D$, interaction $U=0.1D$, and constant hybridization spectrum
$\Gamma(\omega)=0.01D$ for $|\omega|\leq D$. It logarithmically discretizes the
hybridization spectrum, constructs the Wilson chain, and iteratively
diagonalizes the many-body Hamiltonian. With NRG Ljubljana installed, the
calculation can be completed in about ten minutes and requires Git, Bash, and
Make:

```sh
git clone https://github.com/rokzitko/NRG_SIAM.git
cd NRG_SIAM/01_minimal
make
```

All required Hamiltonian and operator matrix elements are supplied in
symmetry-adapted form, so Mathematica is not needed. The calculation produces:

- `run/td`: shell-by-shell thermodynamic quantities and thermal moments of total spin and charge
- `run/custom`: the impurity occupancy `n_d` and double occupancy `n_d_ud`

Inspect the first results with:

```sh
cat run/custom
cat run/td
```

At particle-hole symmetry, `n_d` should remain one within numerical precision.
The subsequent lessons extract the impurity entropy and susceptibility by
subtracting the corresponding impurity-free system, follow the crossover from
a local moment to the screened strong-coupling fixed point, display the
many-body level flow, and calculate the FDM-NRG impurity spectral density with
$z$-averaging.

## Example Result

[Lesson 5](https://github.com/rokzitko/NRG_SIAM/tree/main/05_z_averaging)
calculates the per-spin local impurity spectral function
$A(\omega)=-\mathrm{Im}\,G_{d\sigma}^R(\omega)/\pi$ with FDM-NRG. The
unbroadened, frequency-binned weights from four shifted logarithmic grids,
$z=0.25, 0.50, 0.75, 1.00$, are averaged before a common broadening is applied.
The figure shows $\pi\Gamma A(\omega)$: the narrow Kondo resonance is centered
at the Fermi level, while the Hubbard satellites occur near the atomic charge
excitation energies $\omega=\pm U/2$.

[![Four-twist z-averaged FDM-NRG impurity spectral function](https://raw.githubusercontent.com/rokzitko/NRG_SIAM/main/05_z_averaging/figures/z_averaged_spectrum.svg)](https://github.com/rokzitko/NRG_SIAM/tree/main/05_z_averaging)

No rescaling is used to impose the $T=0$ Friedel-sum-rule limit
$\pi\Gamma A(0)=1$. At $T=10^{-5}D$ and the smallest nonzero mesh frequency,
the four-grid average is $0.9986$. For each $z$, the unbroadened weights satisfy
the zeroth-moment sum rule and the odd moments vanish by particle-hole symmetry;
the averaged second moment is checked against the exact flat-band result. The
lesson notes that $z$-averaging does not replace convergence checks in
$\Lambda$, retained-state truncation, temperature, or broadening width.

## Applications

NRG Ljubljana is used for research on:

- quantum dots, nanostructures, Kondo screening, and quantum phase transitions
- superconducting impurities and Yu-Shiba-Rusinov states
- molecular magnets, surface impurities, and multiorbital impurity models
- impurities in Dirac, Weyl, altermagnetic, and other structured hosts
- dynamical mean-field theory (DMFT) and correlated materials
- equilibrium reference calculations and benchmarks for other many-body methods

For research-grade DMFT(NRG) calculations for strongly correlated electron
systems, see
[DMFT_NRG_KLM](https://github.com/rokzitko/DMFT_NRG_KLM), a Kondo-lattice-model
implementation in which NRG is the DMFT impurity solver. It uses the improved
logarithmic discretization and self-energy estimator, accepts a user-specified
tabulated noninteracting density of states, adjusts the chemical potential to
the desired filling, and calculates single-particle spectra and transport
coefficients.

## Publications Using NRG Ljubljana

The [research bibliography](docs/docs/publications.md) collects a
representative, non-exhaustive selection of applications by independent groups
and project collaborators in quantum impurity physics, superconductivity,
nanostructures, DMFT, and correlated materials.

If you publish work using NRG Ljubljana, please add your paper via a
[pull request](https://github.com/rokzitko/nrgljubljana/edit/master/docs/docs/publications.md)
or [open an issue](https://github.com/rokzitko/nrgljubljana/issues).

## Documentation

The [current documentation](http://auger.ijs.si/nrgljubljana/site/) describes
the definition of an impurity problem, numerical parameters, observables,
normalization conventions, and advanced numerical considerations. Useful
starting points are:

- [specifying the impurity problem and numerical parameters](docs/docs/input-and-configuration.md)
- [complete parameter reference](docs/docs/parameter-reference.md)
- [physical quantities and output conventions](docs/docs/output-formats.md)
- [bath discretization and spectral-analysis programs](docs/docs/tools.md)
- [construction of Hamiltonian and operator matrix elements with `nrginit`](docs/docs/nrginit-workflow.md)
- [single-frequency Floquet formalism, model construction, and truncation](docs/docs/floquet-formalism.md)
- [compilation from source](docs/docs/getting-started.md)
- [parallel execution and numerical-library threading](docs/docs/parallelism.md)

Older scientific notes remain available under `doc/`.

## Citation

If you use results produced with NRG Ljubljana, please cite the method paper
and the archived software release:

- R. Zitko and T. Pruschke, "Energy resolution and discretization artifacts in the numerical renormalization group," *Physical Review B* **79**, 085106 (2009), [doi:10.1103/PhysRevB.79.085106](https://doi.org/10.1103/PhysRevB.79.085106).
- R. Zitko, "NRG Ljubljana" (version 8f90ac4), Zenodo (2021), [doi:10.5281/zenodo.4841076](https://doi.org/10.5281/zenodo.4841076).

[`CITATION.cff`](CITATION.cff) gives the software citation in a standard
citation format.

## Questions and Contributions

- Ask questions about setting up or interpreting a calculation in
  [GitHub Discussions](https://github.com/rokzitko/nrgljubljana/discussions/categories/q-a).
  Reproducible numerical problems can be submitted with the
  [bug report](https://github.com/rokzitko/nrgljubljana/issues/new?template=bug.md).
- See [`CONTRIBUTING.md`](CONTRIBUTING.md) for ways to add a published
  application, example, model, benchmark, documentation improvement, or code
  change.
- Visit the [project home page](http://nrgljubljana.ijs.si/) for additional background and examples.
- Contact Rok Zitko at [rok.zitko@ijs.si](mailto:rok.zitko@ijs.si) for project and research inquiries.

## Advanced: Compilation from Source

Most users should install the conda-forge distribution. Compile from source
when modifying the NRG implementation, enabling CUDA diagonalization, or using
a specific BLAS/LAPACK installation.

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
Mathematica is required for `nrginit`, which constructs the initial Hamiltonian,
symmetry-adapted basis, and operator matrix elements. CUDA support is optional
and requires a CUDA Toolkit that provides cuSOLVER, cuBLAS, and cuDART.

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

### Compilation and Tests

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

Useful CMake options:

- `-DCMAKE_BUILD_TYPE=Debug`
- `-DBuild_Tests=ON|OFF` controls the test build (default: `ON` for a top-level build)
- `-DTEST_LONG=ON`
- `-DASAN=ON -DUBSAN=ON`
- `-DANALYZE_SOURCES=ON`
- `-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON|OFF` enables application-level OpenMP regions such as simultaneous diagonalization scheduling (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_CUDA=ON|OFF` requests CUDA/cuSOLVER support (default: `OFF`)
- `-DNRGLJUBLJANA_ENABLE_MATHEMATICA=ON|OFF` controls `FindMathematica` (default: `OFF` on `aarch64`, `ON` otherwise)
- `-DNRGLJUBLJANA_BLAS_ILP64=ON|OFF` selects the 64-bit or 32-bit BLAS/LAPACK integer ABI (default: `OFF`)
- `-DNRGLJUBLJANA_INSTALL_NRGINIT=ON|OFF` controls installation of the `nrginit` scripts (default: `ON`)
- `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON|OFF` uses preinstalled dependencies instead of CPM downloads (default: `OFF`)

### BLAS/LAPACK Integer ABI

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

### Parallelism Model

NRG Ljubljana normally obtains parallelism from threaded BLAS/LAPACK. A
calculation should use one numerical threading library, either threaded MKL or
threaded OpenBLAS. Application-level OpenMP regions are disabled by default so
that GNU `libgomp` and Intel `libiomp5` are not inadvertently loaded into the
same process.

Use `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and scheduler CPU binding to control numerical kernel threading. When running with MPI, choose the BLAS/LAPACK thread count together with the rank count; `mpi_ranks * blas_threads` should usually not exceed the CPUs allocated to the job.

For MKL builds that use the `mkl_rt` dispatcher, set `-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU`, `INTEL`, or `LLVM` when an explicit threading layer is needed. This links the matching compiler OpenMP runtime through CMake's `OpenMP::OpenMP_CXX` target while keeping application-level OpenMP regions disabled unless `NRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is also set.

`-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON` is an expert option for simultaneous diagonalization scheduling (`diag_mode=OpenMP`, `diagth>1`) and a few non-BLAS loops. It can create nested parallelism when BLAS/LAPACK is also threaded, so CMake checks the visible link line for mixed OpenMP runtime families and the executable prints startup diagnostics and warnings about the detected MKL/OpenBLAS/OpenMP/MPI threading configuration.

## Contents of the Distribution

- `c++/`: iterative diagonalization, truncation, symmetries, operators, density matrices, and numerical routines
- `tools/`: bath discretization, Wilson-chain construction, spectral analysis, and numerical transformations
- `nrginit/`: Mathematica construction of the impurity Hamiltonian, basis, and operators
- `nrgspawn/`: example inputs and reference results prepared with `tools/instantiate/`
- `test/`: numerical checks, reference calculations, and tests involving Mathematica
- `share/`: CMake and environment files installed with NRG Ljubljana
- `scripts/`: utilities for inspecting and combining numerical results
- `doc/`: older scientific documentation
- `docs/`: current documentation

## Contributing

See `CONTRIBUTING.md` for research, documentation, example, benchmark, and code
contributions. Most contributions do not require compiling the source.

## License

NRG Ljubljana is distributed under the GNU General Public License. See `COPYING` for the full license text.

## Acknowledgements

NRG Ljubljana started during Rok Zitko's PhD work at the University of
Ljubljana and the "Jozef Stefan" Institute. The implementation reflects
collaboration and discussions with researchers in the NRG community and
contributions from collaborators over many years.
