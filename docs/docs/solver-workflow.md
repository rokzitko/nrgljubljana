# Native Solver Workflow

The C++20 library `nrgljubljana_workflow` provides file-backed preparation and a
sequential solver entry point. Its only public header is
`<nrgljubljana/solver_workflow.hpp>`, in namespace `NRG::Workflow`. It invokes the
existing discretizer, initializer, chain generator and `NRG::run_nrg` directly;
it does not invoke executables, simulate `main`, or parse diagnostic output.
The standalone tools and their numerical formulas/defaults remain unchanged.

## Consumer Setup

```cmake
find_package(nrgljubljana_workflow CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE nrgljubljana::workflow)
```

The source-tree alias is also `nrgljubljana::workflow`. MPI, solver internals,
BLAS integer-ABI definitions, Eigen and other implementation headers remain
private. A consumer needs only the installed workflow header and target.

Build in the upstream root using the same configuration as the core, especially
the BLAS/LAPACK integer ABI and MKL/OpenMP settings. An existing configured build
can be used without replacing these settings:

```sh
cmake -S . -B build-x86
cmake --build build-x86 --target nrgljubljana_workflow --parallel 2
cmake --install build-x86 --component workflow --prefix /absolute/temporary/prefix
```

Set `nrgljubljana_workflow_DIR` to
`/absolute/temporary/prefix/lib/cmake/nrgljubljana_workflow`, or the corresponding
`lib64` directory on platforms using that layout. The workflow component installs
the matching core binary in `lib[64]/nrgljubljana/workflow`; it does not import a
possibly older installed core package. The package supplies the private core's
link-time dependency path. On Linux, the workflow library uses an origin-relative
RPATH so this core precedes older same-SONAME libraries in `LD_LIBRARY_PATH`.
External runtime dependencies still need the environment used by the upstream
build. Do not link or preload a different NRG core into the same process: an
already loaded library with the same SONAME cannot be replaced by RPATH.

## Preparation

```cpp
#include <nrgljubljana/solver_workflow.hpp>

using namespace NRG::Workflow;
discretize({"/work/case"});
discretize({"/work/case", "param", Branch::negative});
const auto initial = instantiate({"/work/case", "param", "template", 68});
solve("/work/case", "/work/scratch");
```

`calculation_directory` must be nonempty and already exist. A relative
calculation directory is resolved from the caller's CWD. Relative
`parameter_file` and `template_directory` paths are resolved against the
calculation directory, not the parameter file's parent. Absolute paths are
accepted. Other input paths in the files, such as `dos`, retain the tools'
calculation-relative interpretation. Template matrix lookup retains the CLI's
template-directory-first, calculation-directory-fallback behavior.

Discretization requires a regular parameter file (symlinks to regular files are
accepted). Seek/read failures throw before any result is published; they are not
interpreted as missing parameters. An empty readable file or an absent exact
`[param]` block still selects the original tool defaults.

`DiscretizationOptions` defaults to `parameter_file="param"`, positive branch,
and `integral=true`. The latter forces the integral method even if `f_method=ode`
is present, just like `adapt --integral`. With `integral=false`, the parameter
file selects the method and the original ODE default applies if it is absent.
Results are `FSOL.dat` or `FSOLNEG.dat`, and `GSOL.dat` or `GSOLNEG.dat` when the
adaptive mesh is calculated. Both branches must be prepared before initializing
a file-backed bath. Analytic `band=flat` initialization does not require FSOL.

`InitializerOptions` defaults to `parameter_file="param"` and
`template_directory="template"`. **The caller must supply `nmax` in [1,998].**
Its default zero is deliberately invalid. Resolve local `Nmax`, `Tmin` or other
DMFT policies in the caller; the library does not implement that resolver.
The supplied value overrides a copy of the chain parameter map before chain
configuration (including the default `mMAX=2*nmax`) is resolved. It also sets
the generated data header and coefficient lengths. An explicit `mMAX` is
preserved. The original solver parameter file is never edited.

Full initialization currently supports real, one-channel, unpolarized QS
templates. It evaluates seed Hamiltonians, diagonalizes them, transforms the
operator blocks, writes the ground energy and Wilson coefficients, and validates
the generated text with the solver's real `InputData<double>` reader before
publication. Existing linear/Steffen densities, Lanczos/RKPW chain backends,
analytic flat bands, fixed/adaptive meshes, supported hard gaps, and chain-table
load/save settings retain their existing tool semantics. Ordinary ungapped
FSOL interpolation continues to hold its last value beyond the tabulated
`xmax`; it is not truncated at that point. Existing hard-gap extent checks are
unchanged.

Unsupported template symmetries, complex or multichannel templates, polarized
or 2x2 spin coefficients, rungs, substeps, Floquet mode, `tri=cpp`, rescaled
Wilson coefficients and disabled chain tridiagonalization are rejected with
exceptions. The CLI's other existing modes remain available through the CLI.
`InitializerResult.data_file` is the absolute published `data` path;
`ground_energy` is returned directly from seed diagonalization in the same
absolute units as the generated `e` block, not extracted from stdout.

## State And Errors

All three workflow entry points share one serialization mutex and a saved-state
guard, including quantum-number metadata restoration around `solve`. Preparation
restores the caller's CWD, GMP default precision, `std::cout`/`std::cerr`
formatting (flags, precision, width and fill), and quantum-number metadata on
success and exceptions. Scoped GSL handler changes are restored. Each chain
invocation resets its calculation state. Under the guard, `solve` calls the
hardened `NRG::run_nrg`, which provides CWD restoration and fresh `DONE` handling.

This is **not general thread safety**. The host must not overlap other solver
APIs, GSL operations, GMP-default changes, stream-format changes, CWD changes,
or CWD-dependent I/O while a workflow call runs. This includes calls through
other libraries with their own locks. The original CWD must remain accessible
until restoration. Diagnostics can be emitted; the parent may capture logs.

Preparation uses unique staging directories inside the explicit calculation
directory. Intermediate chains and matrices are removed on return or exception.
Completed outputs replace destination symlinks rather than writing through them.
Numerical, template or validation errors leave previously published `data` and
FSOL outputs intact. Each output is published by rename; multiple adaptive
outputs are not a filesystem-wide transaction, so an I/O error during publication
may leave an already-published GSOL alongside the old FSOL. Callers must treat
any exception as a failed preparation and not proceed to solve.

For `solve`, the scratch directory must be absolute and already exist. A unique
temporary child is removed on success and failure. Solver physical outputs
remain in the calculation directory and are not transactional; inspect success
and the fresh `DONE` marker rather than assuming files imply completion.

## Verification

```sh
cmake --build build-x86 --target workflow_contract adapt instantiate nrgchain adapt_high --parallel 2
ctest --test-dir build-x86 -R '^(solver_workflow_contract|adapt_high|adapt[0-9]|instantiate[0-9]|nrgchain[0-9])' --output-on-failure
```

The workflow contract uses disposable copies of the QS template, compares both
FSOL branches and complete initializer output byte-for-byte against the CLIs,
exercises A/B/A and failed calls, checks CWD/GMP/GSL/stream/invariant restoration
and staging cleanup, and runs a small native solver calculation. It covers
directory parameter input, seek/read failures, and a solver exception after
symmetry initialization. Both chain backends
and an ungapped FSOL tail beyond `xmax` are exercised. To test another current
template without touching it, run `tools/solver_workflow/tests/contract.py` with
`--driver`, `--tools`, `--template` and `--temporary-parent` paths.
`tools/solver_workflow/tests/consumer` is a minimal standalone installed-package
consumer requiring no internal headers or explicit MPI dependency.
Pass its executable as `--consumer` to the contract script to include the
installed-library solve check. Adding `--current-input /path/to/calculation
--nmax 68` also tests preparation using disposable copies of that directory's
`param` and `Delta.dat` plus the selected template. The supplied chain length is
explicit, not resolved by the test or library; this check does not run the full
production solver calculation.
