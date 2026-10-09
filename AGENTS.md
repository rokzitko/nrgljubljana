# Repository Guide

## Build And Run

Run from the repository root; CMake rejects a missing/default or relative install prefix:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/nrgljubljana/"
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --timeout 3600 --no-tests=error
```

- Requires CMake >=3.25 and C++20. See `README.md` for system prerequisites; CPM does not supply MPI, BLAS/LAPACK, Boost, GSL, GMP, or HDF5. Default configuration downloads other dependencies; `-DNRGLJUBLJANA_USE_SYSTEM_DEPS=ON` uses preinstalled packages instead.
- Tests use the build-tree library; no installation is needed. Use `build/c++/nrg`, not the root `./nrg` script, which hard-codes a developer's `~/repos/.../build-x86` paths.
- Threaded BLAS/LAPACK supplies normal parallelism; application OpenMP defaults OFF. Keep BLAS/LAPACK and `NRGLJUBLJANA_BLAS_ILP64` on the same integer ABI, and do not mix GNU, Intel, and LLVM OpenMP runtimes. See `docs/docs/parallelism.md`.
- `local_conda_build_test.sh` builds committed HEAD, excluding uncommitted edits; it is not a check of the working tree.

## Code Paths

- Runtime flow: `c++/nrg.cc` -> `c++/nrg-lib.cc` -> `NRG_calculation` in `c++/nrg-general.hpp` -> iteration engine in `c++/core.hpp`. The shared library target is `nrgljubljana_c`; most solver implementation is in headers. Shared numerical changes must support both `double` and `std::complex<double>` instantiations.
- `nrginit/` constructs the initial model data with Mathematica; C++ reads `param` and `data` from the calculation directory. Prepared data allows solver regressions without Mathematica. Parameter handling spans `c++/params.hpp`, `c++/read-input.hpp`, and the initializer, not just the CLI.

## Focused Verification

Build the affected targets before CTest. For example, a single unit-test executable and matching CTest entry:

```sh
cmake --build build --target store --parallel 2
ctest --test-dir build -R '^store$' --output-on-failure --no-tests=error
```

- Unit tests under `test/unit/` get target/test names from their `.cpp` basenames. CMake globs them at configure time: reconfigure after adding files. Regression suites such as `test/c++/` register cases through `simple_tests` and `mpi_tests`; merely adding a directory does not register a test. Inputs are copied into the build tree during configuration.
- Tool-focused run: `ctest --test-dir build -R '^(adapt|nrgchain)' --output-on-failure --no-tests=error`. More focused groups and sanitizer/static-analysis commands are in `CONTRIBUTING.md`.
- `Build_Tests` defaults ON for a top-level build. Long suites require `-DTEST_LONG=ON` and `SYM_ALL`; Mathematica suites require detected Mathematica and `SYM_ALL`. Check enabled coverage rather than interpreting absent tests as passes.
- Independent SIAM validation is opt-in: `-DTEST_SCIENTIFIC=ON`, Python >=3.10 and NumPy >=1.26,<3 in the selected `Python3_EXECUTABLE`. Run `ctest --test-dir build -L '^scientific$' --output-on-failure --no-tests=error`. See `test/scientific/README.md`; prepared fixtures need no Mathematica, but regeneration does.

## Regression Contracts

- Follow `test/README.md`: solver runners clean physical outputs and `DONE`, run in the build-tree case directory, require a fresh `DONE`, then invoke `compare.pl --strict`. Running directly in source fixtures risks overwriting inputs/references or using stale results.
- Do not add golden `.bin`/`.h5` regression outputs; declare semantic validation in `ref/.physical-outputs`. HDF5 validation requires `h5dump` in `PATH`. New physical-output names require coordinated changes to `test/PhysicalOutput.pm` and cleanup/comparator contract tests.
- Keep the standard comparison tolerances unless a demonstrated numerical requirement justifies a narrow `FILENAME.tol` override. Run comparator contracts without a CMake build via `perl test/scripts/comparators.t test`.

## Generated Sources

- Edit `c++/symmetry/nrg-recalc-*.hpp.m4` and shared `recalc-macros.m4`, then run `sh recalc-macros.update` with working directory `c++/symmetry/`. Commit regenerated `.hpp` files with their sources; ordinary CMake builds do not regenerate them.
- For the matrix tool, edit `tools/matrix/parser.yy` or `matrix.ll`, then run `bash update` with working directory `tools/matrix/` (requires Bison/Flex). CMake compiles the checked-in generated parser/scanner rather than regenerating them.
- Numerical reference data under `test/reference/generated/` is generated with exactly `mpmath==1.3.0`; follow `test/reference/README.md`. `python3 test/reference/generate.py --check` checks committed bytes without rewriting. Normal CTest uses committed data and needs no generator environment.

## Style And Docs

- `.clang-format` uses 150 columns, two-space indentation, no include sorting, and no comment reflow. Its `Standard: Cpp11` setting does not describe the actual C++20 build requirement.
- Current docs belong in `docs/docs/`; `doc/` is legacy Sphinx. With `mkdocs` and `pymdown-extensions` installed, validate using `python3 -m mkdocs build --strict -f docs/mkdocs.yml`. `make -C docs upload` deploys via SCP; it is not a validation command.
- Follow `docs/docs/developer-guide.md`: parameter/default/parser/ownership changes update `docs/docs/parameter-reference.md` in the same change; generation-locked semantic changes also update first-run and initializer docs. Output naming, triggers, layouts, units, precision, or binary/HDF5 changes update `docs/docs/output-formats.md`.
- Record feature removals in root-level `REMOVED_FEATURES.md`, newest first, with the removal date (`YYYY-MM-DD`), affected features, and relevant replacement guidance. Keep removal notices there, not in user guides or release summaries; update user guides only to describe current supported behavior.
