# Native Solver Workflow

The C++20 library `nrgljubljana_workflow` provides file-backed preparation and a
sequential solver entry point through `<nrgljubljana/solver_workflow.hpp>`, in
namespace `NRG::Workflow`. In-memory preparation helpers are declared separately
in `<nrgljubljana/solver_preparation.hpp>`. The file-backed workflow invokes the
existing discretizer, initializer, chain generator and `NRG::run_nrg` directly;
it does not invoke executables, simulate `main`, or parse diagnostic output.
The native and standalone paths share the same numerical formulas and defaults.

## Consumer Setup

```cmake
find_package(nrgljubljana_workflow CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE nrgljubljana::workflow)
```

The source-tree alias is also `nrgljubljana::workflow`. MPI, solver internals,
BLAS integer-ABI definitions, Eigen and other implementation headers remain
private. A consumer needs only the installed workflow header and target.

The additive API level is independent of the date-based package/library version
(`2026.09`). Require the diagnostic-routing overloads at compile time with:

```cpp
#include <nrgljubljana/solver_workflow.hpp>
static_assert(NRG::Workflow::workflow_api_version >= 2);
```

Level 3 adds `read_parameter_sections` and, in `<nrgljubljana/solver_chain.hpp>`,
`discretize_hybridization`, `diagonal_coefficients`, `nrgchain_coefficients`, an
`instantiate` overload and `validate_preparation_parameters` with `ChainSource`.

The build-tree and exported `nrgljubljana::workflow` targets also expose the
custom `WORKFLOW_API_VERSION` property, currently `3`, for configure-time checks
with `get_target_property`. Existing option/result layouts and original function
symbols remain unchanged; consumers must link a library providing the new
overloads when using them. The project version has not been renumbered.

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
Its default zero is deliberately invalid. Resolve `Nmax` before invoking the
initializer; the separate `resolve_chain_length` in-memory helper accepts the
selected `[param]` entries. DMFT iteration policy remains with the caller.
The supplied value overrides a copy of the chain parameter map before chain
configuration (including the default `mMAX=2*nmax`) is resolved. It also sets
the generated data header and coefficient lengths. An explicit `mMAX` is
preserved. The original solver parameter file is never edited.

`read_parameter_sections(text)` in the preparation header reads parameter-file
text into its sections, for example to select `[param]` for the helpers above.
It is strict: besides blank lines and whole-line `#` comments, only header lines
that are exactly `[name]` and `key=value` lines inside a section are accepted, and
a section or key given twice is an error naming the line. What it accepts, `nrg`
and the tools read the same way.

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

Both initializer paths fix each real seed eigenvector's sign by making its
first largest-magnitude component positive and normalize exact zeros to positive
zero. All operator blocks use this same basis. This removes arbitrary LAPACK
sign differences, not backend rounding or rotations within degenerate subspaces.

Unsupported template symmetries, complex or multichannel templates, polarized
or 2x2 spin coefficients, rungs, substeps, Floquet mode, `tri=cpp`, rescaled
Wilson coefficients and disabled chain tridiagonalization are rejected with
exceptions. The CLI's other existing modes remain available through the CLI.
`InitializerResult.data_file` is the absolute published `data` path;
`ground_energy` is returned directly from seed diagonalization in the same
absolute units as the generated `e` block, not extracted from stdout.

## Discretization From Gamma In Memory

`discretize_hybridization` in `<nrgljubljana/solver_chain.hpp>` builds the
Wilson chain of a hybridization matrix Gamma(omega) for each twist, with
`mixchain`'s star and chain stages and no files. Gamma is given as
`HybridizationFunction` (row-major N x N per node, real or complex); the result
is one `WilsonChain` per twist, in the units of `chain.dat` (E and T times
`bandrescale`), with the diagnostics of its header as fields. The `[param]`
entries are read with `mixchain`'s names and defaults; `channels`, `dos_prefix`
and `z` are ignored, `Nmax` is replaced by the `nmax` argument (one length for
all twists, `mMAX` defaults to `2*nmax`), and `discretization` must be Z. Gamma
must be positive semidefinite: clipping or projecting it is the caller's. The
call runs under the workflow guard without a calculation directory; a
`diagnostic_file` receives `mixchain`'s log and is resolved against the caller's
working directory, which is not changed. The test `solver_workflow_chain`
compares the result exactly with `chain.dat` from `mixchain --Nz 2`.

`diagonal_coefficients` turns a chain whose channels do not couple into
`ChainCoefficients`: per channel k, `theta` = V_kk^2, `xi` = T_n(k,k) and `zeta` =
E_n(k,k), the coefficients `nrgchain` computes and the initializer reads. Every
element left out (off-diagonal elements, imaginary parts) must be exactly zero,
otherwise it is an error naming the element; a diagonal Gamma discretized with
`split_blocks=true` (the default) gives such a chain. The gauge must be polar.

The `instantiate` overload in `<nrgljubljana/solver_chain.hpp>` takes
`ChainCoefficients` instead of computing the chain with `nrgchain`; everything
else is shared with the `instantiate` above, including its scope (real,
one-channel, unpolarized QS). The coefficients must cover the sites `0..nmax`
(a longer chain is cut), match the coefficient channels of the template, and have
the z, Lambda and bandrescale of the parameter file, to 1e-12 relative. The chain
keys of the file (`dos`, `band`, `adapt`, `tridiag_method`, ...) are not read.
`nrgchain_coefficients(options)` returns the coefficients the `instantiate` above
computes, from the `FSOL` files that `discretize` writes, so both routes give
coefficients to the same `instantiate`; with them it writes the same `data` byte
for byte (test `solver_workflow_instantiate_chain`). `instantiate(options)` is
short for `instantiate(options, nrgchain_coefficients(options))`.

```cpp
// adapt and nrgchain: Delta.dat in each directory
discretize({directory, "param", Branch::positive});
discretize({directory, "param", Branch::negative});
auto coefficients = nrgchain_coefficients({directory, "param", "template", nmax});

// mixchain: Gamma in memory, every twist at once
auto chains       = discretize_hybridization(gamma, param, twists, nmax_all);
auto coefficients = diagonal_coefficients(chains[i]);

instantiate({directory, "param", "template", nmax}, coefficients);
solve(directory, scratch);
```

`validate_preparation_parameters(section, ChainSource::mixchain)` is the static
check of a `[param]` section for this route: `discretization` Z, the chain length
(`Nmax`, or `Tmin` and its keys), z, Lambda, bandrescale, `tri` other than `cpp`,
and the keys `discretize_hybridization` reads, as it reads them. The keys of
`adapt` and `nrgchain` are not read. `ChainSource::nrgchain` is the
one-argument check.

For one channel the chain agrees with `adapt --integral` followed by `nrgchain`
(test `solver_workflow_vs_nrgchain`: flat bands, tabulated densities, a hard gap, a
density on a logarithmic grid). It agrees to rounding when every `x = j + z`
is a node of `FSOL` (z a multiple of `outputstep`, as for Nz = 2 and 4 at the
default 1/64). Otherwise `nrgchain` interpolates `FSOL` linearly, and the two
differ by that interpolation error, quadratic in `outputstep` (up to 3e-5 at
1/64 in the test, at the first sites).

## Diagnostic Routing

Each file-backed entry point has an additional overload accepting a diagnostic
filename, without adding fields to the existing option or result structures:

```cpp
discretize({"/work/case"}, "discretize.log");
const auto initial = instantiate({"/work/case", "param", "template", 68}, "initialize.log");
solve("/work/case", "/work/scratch", "nrg.log");
```

An empty filename is exactly the original stream-routing behavior, as are calls
to the original overloads. Otherwise, the library creates one combined log for
C `stdout`/`stderr` and the descriptor-backed C++ `std::cout`/`std::cerr` streams.
Relative log paths are resolved against `calculation_directory`, not the caller's
CWD or the parameter file's parent. Resolution happens before any workflow CWD
change, including when the calculation directory itself is relative. Absolute
log paths are accepted unchanged.

The log's parent must already exist. Creation uses `O_EXCL` and mode `0600`
(subject to the host umask): an existing file, directory, or symlink, including a
dangling symlink, is rejected without truncating or following it. Use a fresh
filename for each call. A log created before a failed calculation is retained;
it contains diagnostics emitted up to failure, not an automatic rendering of
the thrown exception. Open/redirection failures throw before the backend runs.
No log-routing parameter, environment setting, or callback is introduced, and
physical output files and logging-verbosity parameters are unaffected.

Capture flushes C and C++ output before switching descriptors and again before
restoration, including during exception unwinding. It preserves the caller's
stdout/stderr destinations and descriptor flags. Normal stream buffering still
determines the ordering of interleaved diagnostic messages. Redirection is
**process-wide**, inside the existing workflow serialization and saved-state
guard. It is not per-thread logging or permission to run concurrent application
I/O; unrelated host output during the call can also enter the log. Custom C++
stream buffers installed by the host retain their own routing.

## State And Errors

All three file-backed entry points share one serialization mutex and a saved-state
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
until restoration. Diagnostics can be emitted; use the optional diagnostic-file
overloads to capture them without moving redirection outside the workflow guard.

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

`solve` makes no MPI calls: it runs the solver serially, with threaded
BLAS/LAPACK. A host may initialize MPI and run independent calculations on
different ranks at the same time, each in its own calculation directory; a shared
scratch directory is fine, since each call uses its own child. Forward-only runs
(`fdm=false`) are supported and write no FDM output.

## Verification

```sh
cmake --build build-x86 --target workflow_diagnostics workflow_preparation --parallel 2
ctest --test-dir build-x86 -R '^solver_workflow_(diagnostics(_unsynced)?|preparation)$' --output-on-failure --no-tests=error
cmake --build build-x86 --target workflow_contract workflow_solve_mpi adapt instantiate nrgchain adapt_high --parallel 2
ctest --test-dir build-x86 -R '^(solver_workflow_contract|solver_workflow_solve|adapt_high|adapt[0-9]|instantiate[0-9]|nrgchain[0-9])' --output-on-failure
```

Reconfigure after adding test sources. The focused diagnostics target compiles
the production routing and discretization sources with test-only initializer and
solver stand-ins; it neither links nor runs the full solver. It checks a real
discretization failure after logging, combined C/C++ output, buffered tails,
relative/absolute paths, exclusive creation and symlink rejection, empty-path
and original-overload behavior, and exact descriptor/stream-state restoration
on success and exceptions. Separate processes cover synchronized and
`sync_with_stdio(false)` C++ streams. The preparation target is also solver-free.

The workflow contract uses disposable copies of the QS template, compares both
FSOL branches and complete initializer output byte-for-byte against the CLIs,
exercises A/B/A and failed calls, checks CWD/GMP/GSL/stream/invariant restoration
and staging cleanup, and runs a small native solver calculation. It covers
directory parameter input, seek/read failures, and a solver exception after
symmetry initialization. Both chain backends
and an ungapped FSOL tail beyond `xmax` are exercised. To test another current
template without touching it, run `tools/solver_workflow/tests/contract.py` with
`--driver`, `--tools`, `--template` and `--temporary-parent` paths.
`solver_workflow_solve` and `solver_workflow_solve_np2` solve two cases of
`test/c++` (`mpitest1_qs` with FDM, `legacy_discretization_d1` forward-only)
from a driver that initializes MPI, serially and with one case per rank under
`mpiexec -n 2`, and compare them strictly with the references of the `nrg` tests.
`tools/solver_workflow/tests/consumer` is a minimal standalone installed-package
consumer requiring no internal headers or explicit MPI dependency.
Pass its executable as `--consumer` to the contract script to include the
installed-library solve check. Adding `--current-input /path/to/calculation
--nmax 68` also tests preparation using disposable copies of that directory's
`param` and `Delta.dat` plus the selected template. The supplied chain length is
explicit, not resolved by the test or library; this check does not run the full
production solver calculation.
