# Native DMFT Numerics

`nrgljubljana_numerics` is a small C++20 library, independent of the scientific
solver. Its public header is `include/nrgljubljana/numerics.hpp`; the namespace
is `NRG::Numerics`. Only standard-library types cross the public interface.

## API

```cpp
namespace NRG::Numerics {
struct Table { std::vector<double> x; std::vector<double> y; };
struct Lattice { Table real; Table spectral; };

double integrate(const Table &);
double integrate_abs(const Table &);
double integrate_fermi(const Table &, double temperature,
                       double epsabs = 1e-10, double epsrel = 1e-9);
Table resample(const Table &, const std::vector<double> &grid);
Table kramers_kronig(const Table &);
Lattice lattice(const Table &dos, const Table &re_sigma, const Table &im_sigma,
                double mu, double clip_sigma, unsigned moment = 0);
}
```

All interpolants are Steffen piecewise cubics. Integrals, absolute integrals,
KK transforms, and lattice transforms use the existing analytic polynomial
kernels; finite-temperature Fermi integration uses the existing interval-wise
QAG implementation in `tools/integ/integ.hpp`. There is no new quadrature or
transform algorithm. The implementation uses the library interface of
`tools/kk/kk.hpp` and the analytic tabulated path and safety validation from
`tools/hilb/hilb.hpp` / `tools/common/`. No CLI constructors or parsers run.

## Numerical Contract

- Input coordinates and values must be finite, grids strictly increasing,
  domain widths finite, and x/y lengths equal. Inputs are never sorted or
  modified. Interpolation requires at least three points.
- `integrate_abs` integrates the absolute value of the represented cubic,
  partitioned at its roots. It does not interpolate the absolute sample values.
- `integrate_fermi` uses `f(x)=1/(exp(x/T)+1)` with Fermi level zero. `T=0`
  uses the analytic negative-energy integral. Finite `T>=0` and valid finite
  nonnegative GSL tolerances are required. If absolute tolerance is zero,
  relative tolerance must be at least `max(50*epsilon(double), 0.5e-28)`.
  The inherited kernel checks the summed error against the global target;
  failure throws rather than printing a warning or accepting an inaccurate sum.
- `resample` requires a nonempty finite strictly increasing target grid.
  Extrapolation is rejected beyond
  `32*epsilon(double)*max(abs(a),abs(b),b-a)` of the source domain `[a,b]`.
  Smaller excursions use the endpoint value, retaining the requested x label.
- KK retains the CLI convention `PV integral S(E)/(E-w) dE / pi`, converting
  imaginary input to real output. The grid must be even-sized and symmetric
  about zero, using the CLI's `gsl_fcmp` symmetry tolerance `1e-8`. Function
  values need not be even or odd. Endpoints return the CLI's finite subtracted
  integral `integral (S(E)-S(w))/(E-w) dE / pi`, not a divergent one-sided limit.
- `lattice` computes `G_n(w) = integral E^n DOS(E)/(z(w)-E) dE`, with
  `z(w) = w + mu - ReSigma(w) + i*max(-ImSigma(w),clip_sigma)`. The DOS is
  nonnegative, need not be normalized or symmetric, and is zero outside its
  tabulated support. The self-energy tables need at least one point, must have
  exactly matching frequency grids, and are not interpolated.
- **`Lattice::real` is raw `Re G_n`; `Lattice::spectral` is `-Im G_n/pi`.**
  This deliberately differs from default CLI `hilb`, which divides both output
  columns by `-pi`. Output grids are the self-energy grid. For `n>0` these are
  energy-weighted correlators, and odd moments can have signed spectral weight.
  Energy weighting multiplies the DOS polynomial before transformation; in
  particular `n=1` never uses the cancellation-prone `z*G_0-mass` identity.
- `mu` must be finite and `clip_sigma >= sqrt(min normal double)` must be
  finite, matching `hilb`'s safety floor (about `1.4916681462400413e-154`). The
  CLI clipping rule intentionally also clips positive imaginary self-energies.
  Moments must fit in `int`; high degrees may overflow coefficients or exhaust
  memory. No normalization or artificial positivity is imposed on results.
- Signed data is accepted by integration, Fermi integration, resampling, and
  KK. Negative DOS samples are rejected. Zero DOS is allowed.

Errors propagate as standard C++ exceptions, including invalid arguments,
out-of-domain requests, numerical overflow, and quadrature failure. Existing
KK validation sometimes throws `std::runtime_error` for malformed meshes.
Nonfinite numerical outputs are rejected. The library does not print, exit,
read environment worker settings, or perform file I/O.

GSL operations use the shared `NRG::Tools::GslErrorHandlerGuard`, which restores
the previous process-global handler on success and exceptions, including the
default aborting handler (`nullptr`) and the distinct disabled handler. Its recursive
mutex serializes cooperating GSL regions; independent callers must also use
that shared mechanism rather than concurrently changing GSL's global handler.
Analytic transforms need no GSL handler changes after interpolation setup.
Each call prepares its own interpolant; there is currently no public cache or
`PreparedDos` object.

## Build And Test

From the repository root, with an out-of-source temporary build directory:

```sh
cmake -S tools/numerics -B /tmp/nrgljubljana-numerics-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/tmp/nrgljubljana-numerics-install \
  -DNRGLJUBLJANA_NUMERICS_BUILD_TESTS=ON
cmake --build /tmp/nrgljubljana-numerics-build --parallel 2
ctest --test-dir /tmp/nrgljubljana-numerics-build --output-on-failure --no-tests=error
cmake --install /tmp/nrgljubljana-numerics-build
```

Requires CMake >=3.25, C++20, GSL >=2.1, Boost headers, and platform threads.
No Eigen, MPI, BLAS/LAPACK solver stack, dependency downloads, or test framework
are needed. Standalone here means building this subtree within a repository
checkout: the sibling shared kernels are still needed at build time. The
installed package does not depend on the checkout or Boost headers. Existing
private CLI headers use POSIX headers; Linux is the supported/tested platform.
`BUILD_SHARED_LIBS=ON` selects a shared library instead of the default static
library. Strict floating-point flags are applied privately to the library.

The six CTest groups cover flat/linear and cubic integration oracles, absolute
interpolant roots, finite/zero-temperature Fermi integrals, Steffen resampling,
KK sign/parity/endpoints, lattice moments and far-field cancellation, malformed
input, clipping, and GSL handler restoration including concurrent calls.

## Consume

Installed package:

```cmake
find_package(nrgljubljana_numerics 1 CONFIG REQUIRED)
target_link_libraries(dmft PRIVATE nrgljubljana::numerics)
```

Set `CMAKE_PREFIX_PATH` to the installation prefix. The standalone config
package is named `nrgljubljana_numerics`, not the solver's `nrgljubljana` package.
GSL and Threads are resolved by the config for static-library linking; Boost
is a build-only header dependency.

The normal upstream build includes this library through `tools/CMakeLists.txt`
and installs the same `nrgljubljana_numerics` package. Its export set remains
separate from `nrgljubljana-targets`, with no duplicate installed target names.

Alternatively, an external parent can wire just the subtree directly:

```cmake
add_subdirectory(path/to/nrgljubljana/tools/numerics numerics-build)
target_link_libraries(dmft PRIVATE nrgljubljana::numerics)
```

The build target is `nrgljubljana_numerics`; its alias and installed imported
target are both `nrgljubljana::numerics`. In the upstream build the focused tests
default on when the root tests are enabled. In other parent projects they default
off; enable `NRGLJUBLJANA_NUMERICS_BUILD_TESTS` explicitly and have the parent
enable CTest if desired.

The minimal consumer under `tests/consumer` checks both integration modes:

```sh
cmake -S tools/numerics/tests/consumer -B /tmp/numerics-package-consumer \
  -DCMAKE_PREFIX_PATH=/tmp/nrgljubljana-numerics-install
cmake --build /tmp/numerics-package-consumer --parallel 2
ctest --test-dir /tmp/numerics-package-consumer --output-on-failure

cmake -S tools/numerics/tests/consumer -B /tmp/numerics-subdir-consumer \
  -DNRGLJUBLJANA_NUMERICS_SOURCE_DIR="$PWD/tools/numerics"
cmake --build /tmp/numerics-subdir-consumer --parallel 2
ctest --test-dir /tmp/numerics-subdir-consumer --output-on-failure
```

To check coexistence with an installed solver package, add its prefix to
`CMAKE_PREFIX_PATH` and configure the package consumer with
`-DNRGLJUBLJANA_NUMERICS_CHECK_SOLVER_PACKAGE=ON`. The consumer still links only
the numerical library and does not execute the solver.
