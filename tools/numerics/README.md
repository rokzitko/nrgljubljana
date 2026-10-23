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

struct SpectrumBin { double frequency; double weight; };
using SpectrumBins = std::vector<SpectrumBin>;
enum class SpectrumLayout { real_pairs, complex_real, complex_imaginary };
enum class BroadeningKernel {
  hybrid_output_frequency, hybrid_peak_frequency, gaussian
};
struct BroadeningOptions {
  double alpha;
  double temperature;
  double omega0_ratio = 1e-9;
  BroadeningKernel kernel = BroadeningKernel::hybrid_output_frequency;
  double accumulation = 0;
  double filter_low = 0;
  double filter_high = std::numeric_limits<double>::max();
  bool keep_nonnegative = true;
  bool keep_negative = true;
  std::optional<double> final_gaussian_gamma;
  std::optional<double> final_derfd_gamma;
  std::size_t jobs = 1;
};
SpectrumBins load_spectrum_bins(const std::string &,
                               SpectrumLayout = SpectrumLayout::real_pairs);
Table broaden(const std::vector<SpectrumBins> &, const std::vector<double> &mesh,
              const BroadeningOptions &);
Table broaden_files(const std::vector<std::string> &, const std::vector<double> &mesh,
                    const BroadeningOptions &,
                    SpectrumLayout = SpectrumLayout::real_pairs);
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

- Table coordinates and values must be finite, grids strictly increasing,
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
read environment worker settings, launch subprocesses, parse CLI arguments, or
change the working directory. **`load_spectrum_bins` and `broaden_files` are the
explicit exceptions to the no-file-I/O contract:** they only read the named
binary inputs. All other entry points, including `broaden`, work in memory.

GSL operations use the shared `NRG::Tools::GslErrorHandlerGuard`, which restores
the previous process-global handler on success and exceptions, including the
default aborting handler (`nullptr`) and the distinct disabled handler. Its recursive
mutex serializes cooperating GSL regions; independent callers must also use
that shared mechanism rather than concurrently changing GSL's global handler.
Analytic transforms need no GSL handler changes after interpolation setup.
Each call prepares its own interpolant; there is currently no public cache or
`PreparedDos` object.

## Broadening Contract

Broadening shares `tools/broaden/broadening_core.hpp` with the CLI, including
exact-key averaging, filtering, `std::transform_reduce` accumulation and final
convolution. Primary kernel formulas remain in `c++/broadening.hpp`, also used
by the solver. The native implementation does not include the CLI header and
never accesses or changes the GSL error handler. Calls own all mutable state;
independent calls and `jobs > 1` are supported, with the same per-destination
arithmetic regardless of worker count. `jobs` defaults to 1, not an environment
setting. The CLI retains its existing options, defaults and reporting.

- Binary files contain native-endian, native-representation `double` records,
  without a header: pairs `(frequency, weight)` for `real_pairs`, or triplets
  `(frequency, real, imaginary)` for either complex layout. The selected
  component becomes the signed weight. Record lengths and **all components**,
  including the unselected component, are validated. Truncation, failed I/O,
  nonfinite values and invalid layouts throw; empty files are valid.
- Loading preserves input order and duplicates. Broadening accepts unsorted
  bins, signed weights and empty individual contributions, but the outer list
  of spectra or filenames must be nonempty. In list order, and then record
  order, weights are added to a `std::map<double,double>` keyed by exact
  frequency. Each accumulated value is divided **once** by the number of
  contributions, including empty ones. No per-file normalization, tolerance
  coalescing, positive-weight clipping or total-weight normalization occurs.
- The caller supplies a nonempty, finite, strictly increasing destination mesh,
  retained exactly in `Table::x`. One- and two-point meshes are allowed.
  Hybrid kernels reject zero destination frequencies instead of producing a
  log-Gaussian NaN; Gaussian broadening accepts zero. Zero-frequency input bins
  are valid for every kernel.
- `alpha`, `temperature`, `omega0_ratio`, and each supplied final gamma must be
  finite and positive, as must the products `omega0_ratio*temperature` and
  `gamma*temperature`. `accumulation` must be finite, `jobs >= 1`, and filters
  must satisfy finite `0 <= filter_low <= filter_high`. Nonfinite accumulated
  weights or numerical results throw rather than returning NaNs or infinities.
- Filters zero weights with `abs(frequency) < filter_low` or
  `abs(frequency) > filter_high`, retaining both boundary values. Sign filters
  operate on the **frequency**, never the weight; zero is nonnegative. Setting
  both sign flags false is allowed and yields zero on ordinary valid meshes.
- `hybrid_output_frequency` is the default CLI crossover. The alternative
  `hybrid_peak_frequency` corresponds to `-n`; it is a kernel choice, not a
  post-hoc normalization. Both use `omega0 = omega0_ratio*temperature` and the
  CLI accumulation-point shift only in the logarithmic term. `gaussian` matches
  `-g`, using width `alpha` (standard deviation `alpha/sqrt(2)`) and ignoring
  accumulation in its kernel.
- Final Gaussian (`-f`) and derivative-of-Fermi (`-x`) passes use
  `sigma = gamma*temperature`. The Gaussian's standard deviation is `sigma`;
  the other kernel is `1 / (2*sigma*(1+cosh((x-y)/sigma)))`. When both are set,
  Gaussian runs first. Each pass snapshots its input, sums **interior** source
  points using the exact legacy cell width
  `(mesh[j+1]+mesh[j])/2 - (mesh[j-1]+mesh[j])/2`, and does **not** normalize the
  kernel sum. Endpoints remain untouched. Only destinations satisfying strict
  `abs(w) < 100*sigma` are convolved; the source points have no cutoff.

Equivalent to the existing script command
`broaden -x gamma -L mesh basename Nz alpha T 1e-9`, with an already loaded mesh:

```cpp
NRG::Numerics::BroadeningOptions options{.alpha = alpha, .temperature = T};
options.final_derfd_gamma = gamma;
auto spectrum = NRG::Numerics::broaden_files(filenames_in_z_order, mesh, options);
```

Pass paths explicitly, for example `{"1/spec.bin", "2/spec.bin"}`. The API
neither constructs numbered directories nor loads a mesh from a text file.

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

The original six CTest groups cover flat/linear and cubic integration oracles, absolute
interpolant roots, finite/zero-temperature Fermi integrals, Steffen resampling,
KK sign/parity/endpoints, lattice moments and far-field cancellation, malformed
input, clipping, and GSL handler restoration including concurrent calls.
Six additional broadening groups cover independent kernel/convolution oracles,
signed symmetry, exact-key averaging and its rounding order, all binary layouts,
malformed records/options/meshes, empty contributions, cutoff boundaries,
unmodified endpoints, final-pass order, and repeated/concurrent calls.

Optionally set
`-DNRGLJUBLJANA_NUMERICS_BROADEN_CLI=/absolute/path/to/broaden` to add CLI parity
testing (requires Python 3). The harness alone launches subprocesses, in fresh
temporary directories, and compares ten signed cases at 1, 3 and 20 workers.
An existing pre-extraction executable can serve as an independent oracle;
the numerical library itself never launches a CLI. No root CMake linkage
change is needed for `broaden`, which uses the shared header directly.

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
