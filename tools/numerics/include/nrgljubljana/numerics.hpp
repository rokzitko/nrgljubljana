#ifndef NRGLJUBLJANA_NUMERICS_HPP
#define NRGLJUBLJANA_NUMERICS_HPP

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
// Failures throw standard C++ exceptions; no CLI parsing, file I/O, or exits.

// Integrate the represented interpolant over its full tabulated domain.
double integrate(const Table &table);

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

} // namespace NRG::Numerics

#endif
