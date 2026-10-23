#ifndef NRGLJUBLJANA_NUMERICS_HPP
#define NRGLJUBLJANA_NUMERICS_HPP

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace NRG::Numerics {

struct Table {
  std::vector<double> x;
  std::vector<double> y;
};

struct Lattice {
  Table real;     // Raw Re G_n, NOT -Re G_n/pi.
  Table spectral; // A_n = -Im G_n/pi; signed for odd moments.
};

// Interpolated tables require at least three finite points with matching x/y
// sizes, strictly increasing x, and finite domain width. Interpolation is
// always Steffen. Signed ordinates are supported except for the physical DOS.
// Failures throw standard C++ exceptions; no CLI parsing or exits. Only the
// explicitly named spectrum-loading functions below perform file I/O.

// Integrate the represented interpolant over its full tabulated domain.
double integrate(const Table &table);

// Integrate over [lower,upper], clipped to the represented table domain.
double integrate_range(const Table &table, double lower, double upper);

// Weight the represented interpolant, not the tabulated knots, by x^moment.
double integrate_moment(const Table &table, unsigned moment);

// Integrate |S(x)|, splitting at interpolant roots, NOT a spline through |y|.
double integrate_abs(const Table &table);

// Integrate S(x)/(exp(x/temperature)+1), with the Fermi level at x=0.
// Temperature must be finite and nonnegative; zero uses the exact step limit.
// Tolerances must be finite and nonnegative. If epsabs=0, epsrel must be at
// least max(50*epsilon(double), 0.5e-28). Finite-T quadrature errors are checked
// against the global max(epsabs, epsrel*abs(result)) target.
double integrate_fermi(const Table &table, double temperature, double epsabs = 1e-10, double epsrel = 1e-9);

// The target grid must be nonempty, finite, and strictly increasing. Reject
// extrapolation except within 32*epsilon(double)*max(|a|,|b|,b-a) of [a,b].
// Such roundoff excursions use the endpoint value; output x preserves grid.
Table resample(const Table &table, const std::vector<double> &grid);

// CLI convention: K(w) = PV integral S(E)/(E-w) dE / pi, so Im input yields Re.
// Requires an even-sized mesh symmetric about zero (CLI relative tolerance
// 1e-8); ordinates need not have parity. Output uses the input mesh, including
// endpoints, where K_sub(w) = integral (S(E)-S(w))/(E-w) dE / pi is returned.
Table kramers_kronig(const Table &table);

// G_n(w) = integral E^moment * DOS(E)/(z(w)-E) dE over the DOS domain,
// z(w) = w + mu - ReSigma(w) + i*max(-ImSigma(w), clip_sigma).
// DOS must be nonnegative, but need not be normalized or symmetric. Moments
// use direct polynomial energy weighting, including n=1 (never z*G_0-mass).
// Self-energy tables need at least one point and exactly matching grids;
// neither is interpolated. Both output tables use that frequency grid.
// mu must be finite; clip_sigma must be finite and >= sqrt(min normal double).
// As in hilb, positive ImSigma is deliberately clipped too. The moment must
// fit in int, as required by hilb; large degrees may exhaust memory or overflow.
Lattice lattice(const Table &dos, const Table &re_sigma, const Table &im_sigma, double mu, double clip_sigma,
                unsigned moment = 0);

struct SpectrumBin { double frequency; double weight; };
using SpectrumBins = std::vector<SpectrumBin>;
enum class SpectrumLayout { real_pairs, complex_real, complex_imaginary };
enum class BroadeningKernel { hybrid_output_frequency, hybrid_peak_frequency, gaussian };

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

// Native-endian binary doubles: (frequency,weight) or (frequency,real,imag).
// Validate whole records and all components, selecting the requested weight.
// Empty files are valid. No sorting, coalescing, or normalization on loading.
SpectrumBins load_spectrum_bins(const std::string &, SpectrumLayout = SpectrumLayout::real_pairs);

// Exact-key map averaging in input order, divided once by the outer size
// (which must be nonzero, including empty contributions). Signed weights and
// unsorted duplicate bins are supported. The explicit mesh must be nonempty,
// finite and strictly increasing; hybrid kernels reject zero destinations.
// All options must be finite; alpha, temperature, omega0_ratio, optional gamma
// and their width products must be positive. Require 0<=filter_low<=filter_high
// and jobs>=1. Filters include both boundaries and classify zero as nonnegative.
// Final passes use sigma=gamma*temperature, Gaussian then derivative of Fermi,
// interior source cell widths, unchanged endpoints, and |destination|<100*sigma.
// No normalization, clipping, environment access, logging, cwd or GSL changes.
Table broaden(const std::vector<SpectrumBins> &, const std::vector<double> &mesh, const BroadeningOptions &);
Table broaden_files(const std::vector<std::string> &, const std::vector<double> &mesh, const BroadeningOptions &,
                    SpectrumLayout = SpectrumLayout::real_pairs);

} // namespace NRG::Numerics

#endif
