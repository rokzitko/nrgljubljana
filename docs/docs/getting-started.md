# Getting Started

## Dependencies

Source builds require CMake 3.25 or newer and C and C++ compilers with
C++20 support. These dependencies must be installed on the system regardless
of how the remaining C++ dependencies are resolved:

- threaded BLAS and LAPACK, typically MKL or OpenBLAS
- MPI with C++ support
- Boost 1.68 or newer, including the MPI and serialization components
- GSL 2.0 or newer
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

### BLAS/LAPACK integer ABI

The BLAS integer ABI is independent of the operating system's pointer size.
The usual LP64 interface passes 32-bit integers to BLAS and LAPACK, while the
ILP64 interface passes 64-bit integers. NRG Ljubljana defaults to LP64 through
`NRGLJUBLJANA_BLAS_ILP64=OFF`; enabling the option selects ILP64. The option
also fixes CMake's `BLA_SIZEOF_INTEGER` to `4` or `8`, and configuration rejects
an explicitly contradictory value.

The selected BLAS and LAPACK libraries, the project declarations, and any
dispatcher interface must all use the same integer width and symbol naming
convention. An ABI mismatch may still link successfully but passes arguments
with the wrong layout and is unsafe. Typical explicit selections are:

```sh
# LP64 OpenBLAS
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ \
  -DNRGLJUBLJANA_BLAS_ILP64=OFF -DBLA_VENDOR=OpenBLAS

# LP64 Intel MKL
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ \
  -DNRGLJUBLJANA_BLAS_ILP64=OFF -DBLA_VENDOR=Intel10_64lp

# ILP64 Intel MKL
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ \
  -DNRGLJUBLJANA_BLAS_ILP64=ON -DBLA_VENDOR=Intel10_64ilp
```

For ILP64, prefer MKL's direct `Intel10_64ilp` interface. An `mkl_rt`
dispatcher build is only compatible when its runtime interface layer is also
set to ILP64; NRG Ljubljana does not select that layer. OpenBLAS is compatible
with ILP64 only when the installed library was built for 64-bit integers and
exports a symbol convention detected by this project. Apple Accelerate is an
LP64 choice with the supported CMake versions. ILP64 declarations currently
use a 64-bit C++ `long`, so the option is not supported on LLP64 platforms such
as Windows. The ILP64 path is covered in CI on Linux x86-64 with Intel MKL.

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
- `-DNRGLJUBLJANA_BLAS_ILP64=ON|OFF` selects the 64-bit or 32-bit BLAS/LAPACK integer ABI (default: `OFF`)
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

## First End-To-End Calculation

The following calculation is the same minimal single-impurity Anderson model
used by the `nrginit0_minimal` regression fixture. It exercises both installed
entry points and requires a working Mathematica or Wolfram Engine kernel.

Start in a new directory. `nrginit` and `nrg` both read and write the current
directory, and optional files left by an older calculation are not cleaned up
automatically.

```sh
run_dir="$(mktemp -d "$HOME/nrg-minimal.XXXXXX")"
cd "$run_dir"
cat > param <<'EOF'
[param]
symtype=QS
model=SIAM

U=1.0
Gamma=0.1
delta=0.1
band=flat

discretization=Z
Lambda=2.0

Nmax=4
keep=100

ops=n_d
EOF
```

Generate the seed Hamiltonian, operators, and Wilson-chain coefficients:

```sh
nrginit && test -s data
```

The launcher uses the command `math` by default. If the kernel has another
name or is not on `PATH`, pass it as the only positional argument:

```sh
nrginit /path/to/WolframKernel
```

A successful initializer exits with status zero, prints `Success!`, and leaves
`data` and `mmalog` in the calculation directory. `data` is the generated input
for `nrg`; `mmalog` is a diagnostic Mathematica transcript.

Run the iterative solver only after `nrginit` succeeds:

```sh
nrg && test -f DONE
```

For this input, the main persistent results are:

- `td`: one row of thermodynamic quantities for the seed problem and each NRG
  shell
- `custom`: expectation values for requested scalar operators, including
  `n_d`
- `DONE`: an empty completion marker created at the end of a successful fresh
  run when `done=true`

Inspect the headers and first data rows with:

```sh
sed -n '1,2p' td
sed -n '1,3p' custom
```

The exact floating-point values can vary with the eigensolver and numerical
libraries. Use the headers rather than fixed character offsets when parsing the
tables. See [Output format reference](output-formats.md) for column meanings and
file lifecycle details.

The runtime also creates a uniquely named temporary workdir for eigenspectra
and density matrices used between phases. `nrg -w DIR` or `NRG_WORKDIR=DIR`
selects the parent of that temporary directory; neither setting relocates
`td`, `custom`, spectra, or other persistent results.

For restartable calculations, set `resume=true` in `param` before the initial
run and invoke `nrg --checkpoint-dir DIR`. Unlike `-w`, this option creates or
reopens exactly `DIR`, locks it against concurrent solver processes, and
preserves the directory after exit. Repeat the same command after an
interruption; completed diagonalizations are loaded while all shells are
replayed to reconstruct complete output files. The checkpoint files are native
same-build artifacts and require unchanged `param` and `data` inputs.

It is normal for `nrg` to report initializer-only keys such as `model`, `U`,
`Gamma`, `band`, and `Nmax` under `Unused settings`. Review the list because a
misspelled runtime key appears in the same place. The complete ownership and
syntax rules are in the [parameter reference](parameter-reference.md).

## Legacy Documentation

The older Sphinx documentation remains in `doc/` while this MkDocs tree is
being expanded. Use this MkDocs site for the current build, workflow, parameter,
and output references as well as contributor-facing code orientation.
