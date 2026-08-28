# Conda build notes

This document describes how to build `nrgljubljana` with conda using dependencies from conda-forge. The `nrginit` scripts are installed by default; using them requires Mathematica at runtime.

## Recipe layout

The conda recipe lives in:

- `recipe/meta.yaml`
- `recipe/build.sh`

Build with:

```bash
conda build -c conda-forge recipe
```

## Prerequisites

1. Conda tooling:

```bash
conda install -n base -c conda-forge conda-build
```

2. Optional: a local Mathematica installation on the machine, needed to run
   `nrginit` and Mathematica-dependent tests.

## Optional: enable Mathematica detection

The recipe always installs the `nrginit` scripts, but configure-time
Mathematica detection defaults to `OFF` on every platform through the
`nrgljubljana_enable_mathematica` variant. CMake root or kernel hints only
control where `FindMathematica` searches; they do not enable detection. The
recipe passes its explicit `NRGLJUBLJANA_ENABLE_MATHEMATICA` value after
`CMAKE_ARGS`, so setting that option in `CMAKE_ARGS` is also insufficient.

Enable detection with the recipe variant:

```bash
conda build -c conda-forge \
  --variants "{nrgljubljana_enable_mathematica: ['ON']}" \
  recipe
```

CTest is disabled in Conda builds by default. To detect Mathematica and run
the Mathematica-dependent tests, enable both variants:

```bash
conda build -c conda-forge \
  --variants "{nrgljubljana_enable_mathematica: ['ON'], nrgljubljana_build_tests: ['ON']}" \
  recipe
```

The finder recognizes a standard installation through the `math` command. For
a nonstandard native installation, pass the host root and kernel location in
addition to enabling the variant:

```bash
export CMAKE_ARGS="${CMAKE_ARGS:-} \
  -DMathematica_HOST_ROOT_DIR=/path/to/Wolfram \
  -DMathematica_ROOT_DIR=/path/to/Wolfram \
  -DMathematica_KERNEL_EXECUTABLE=/full/path/to/WolframKernel"
```

Build-time detection does not store the kernel path in the installed launcher.
At runtime, `nrginit` invokes `math` by default; pass a different executable as
its only argument, for example `nrginit /full/path/to/WolframKernel`.

## Build the package

From repository root:

```bash
conda build -c conda-forge recipe
```

Artifacts are placed in your local conda-bld directory (for example under `~/miniconda3/conda-bld/`).

By default, compilation and CTest use conda-build's `CPU_COUNT`. To override the number of build jobs:

```bash
conda build -c conda-forge --variants "{nrgljubljana_build_jobs: ['64']}" recipe
```

CTest is disabled during conda builds by default. To build and run the non-long CTest suite during packaging:

```bash
conda build -c conda-forge --variants "{nrgljubljana_build_tests: ['ON']}" recipe
```

To build with one job count and test with another:

```bash
conda build -c conda-forge --variants "{nrgljubljana_build_tests: ['ON'], nrgljubljana_build_jobs: ['64'], nrgljubljana_test_jobs: ['8']}" recipe
```

CTest uses a one-hour per-test timeout by default. To use a longer timeout for slow or debug builds:

```bash
conda build -c conda-forge --variants "{nrgljubljana_build_tests: ['ON'], nrgljubljana_test_timeout: ['14400']}" recipe
```

Linux `aarch64` conda builds should use Clang:

```bash
conda build -c conda-forge --variants "{c_compiler: ['clang'], cxx_compiler: ['clangxx']}" recipe
```

On NVIDIA Grace systems, build the Linux `aarch64` package with NVPL BLAS/LAPACK:

```bash
conda build -c conda-forge --no-anaconda-upload \
  --variants "{c_compiler: ['clang'], cxx_compiler: ['clangxx'], blas_impl: ['nvpl'], nrgljubljana_build_tests: ['ON']}" \
  recipe
```

To include long tests as well:

```bash
conda build -c conda-forge --variants "{nrgljubljana_build_tests: ['ON'], nrgljubljana_test_long: ['ON']}" recipe
```

## Install from local build

```bash
conda create -n nrg-local -c local -c conda-forge nrgljubljana
conda activate nrg-local
```

To reinstall a rebuilt local package into an existing environment:

```bash
conda install -n nrg-local -c local -c conda-forge --force-reinstall nrgljubljana
```

## Verify install

Basic binary/library checks:

```bash
command -v nrg
command -v adapt
test -f "${CONDA_PREFIX}/lib/cmake/nrgljubljana/nrgljubljana-config.cmake"
```

Activation hook:

```bash
echo "${NRGLJUBLJANA_ROOT}"
printf '%s\n' "${PATH}" | tr ':' '\n' | grep -Fx "${CONDA_PREFIX}/nrginit"
```

It should print the active conda prefix and show that `nrginit/` is on `PATH`.

## Verify nrginit installation

The `nrginit` scripts are installed by default. Verify the launcher is present and on `PATH` after activation:

```bash
test -x "${CONDA_PREFIX}/nrginit/nrginit"
command -v nrginit
```

If Mathematica detection was enabled and succeeded, the logs contain a line
similar to `Mathematica executable ...`. Mathematica-dependent tests are only
registered when tests and the full symmetry set are enabled and Mathematica is
found.

## Notes

- Current recipe skips Windows.
- The conda recipe uses Intel MKL on `linux-64`, Accelerate on macOS, and OpenBLAS on generic Linux `aarch64`; NVPL is available as an explicit Linux `aarch64` variant.
- The provided Conda variants use the LP64 BLAS/LAPACK integer ABI. `blas_impl` selects an implementation, not ILP64; the recipe does not currently define or test an ILP64 package variant.
- Mathematica is not provided by conda dependencies; using `nrginit` relies on a local system installation.
- Mathematica detection is disabled by default; `nrginit` is still installed, but Mathematica-dependent tests are not registered.
