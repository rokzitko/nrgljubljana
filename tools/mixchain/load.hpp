// Channel-mixing discretization for NRG
// ** Loading of the matrix of the spectral function Gamma(omega) of the hybridisation function

#ifndef _mixchain_load_hpp_
#define _mixchain_load_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../common/io.hpp"
#include "../common/tabulated.hpp"
#include "linint.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// The input file names carry the matrix indices as single digits
constexpr int max_channels = 9;

// The node added at omega=0, as add_zero_point() does in adapt and nrgchain.
constexpr double zero_point = 1e-99;

inline auto gamma_filename(const std::string &prefix, const int i, const int j, const bool imaginary) {
  return prefix + "_" + std::to_string(i) + std::to_string(j) + (imaginary ? "-im.dat" : "-re.dat");
}

struct GammaOptions {
  std::string prefix{"Gamma"};
  int channels{1};
  double bandrescale{1.0};
  double hermiticity_tolerance{1e-8};
};

// One frequency branch of the input: |omega| strictly increasing, with the Hermitian matrix Gamma at each node.
template<typename S> struct GammaBranch {
  std::vector<double> omega;
  std::vector<Matrix<S>> gamma;
  // The smallest |omega| that was actually tabulated, before the node at omega=0 was added. Below it the density is
  // the constant continuation of the input; zero means that nothing was added.
  double innermost{};
  auto size() const { return omega.size(); }
};

// Gamma(omega) on the input grid, split into the positive and the negative frequency branch.
template<typename S> struct GammaInput {
  int channels{};
  GammaBranch<S> pos, neg;
};

namespace detail {

using Pairs = std::vector<std::pair<double, double>>;

inline void validate_options(const GammaOptions &options) {
  if (options.channels < 1 || options.channels > max_channels)
    throw std::invalid_argument("channels must be between 1 and " + std::to_string(max_channels) + ".");
  if (!(std::isfinite(options.bandrescale) && options.bandrescale > 0.0))
    throw std::invalid_argument("bandrescale must be a positive finite number.");
  if (!(std::isfinite(options.hermiticity_tolerance) && options.hermiticity_tolerance > 0.0))
    throw std::invalid_argument("The Hermiticity tolerance must be a positive finite number.");
}

// Read one two-column component file. Unlike load_rho() in adapt, the values are signed: an off-diagonal element of
// Gamma is not a density.
inline auto read_column_file(const std::string &filename) {
  std::ifstream F;
  NRG::Tools::open_input(F, filename);
  auto pairs = NRG::Tools::read_strict_pairs(F, filename);
  if (pairs.size() < 2) throw std::runtime_error(filename + ": at least two data points are required.");
  return pairs;
}

inline void check_strictly_increasing(const Pairs &pairs, const std::string &filename) {
  for (std::size_t k = 1; k < pairs.size(); k++) {
    if (!(pairs[k].first > pairs[k - 1].first))
      throw std::runtime_error(filename + ": frequencies must be strictly increasing; found "
                               + std::to_string(pairs[k].first) + " after " + std::to_string(pairs[k - 1].first) + ".");
  }
}

// Every component must be tabulated on one and the same grid: this tool does not reinterpolate the input. The
// comparison is not bitwise, so that the same nodes written at different output precision are still accepted.
inline void check_same_grid(const Pairs &reference, const std::string &reference_name, const Pairs &other,
                            const std::string &name) {
  if (other.size() != reference.size())
    throw std::runtime_error(name + ": expected the same frequency grid as " + reference_name + "; found "
                             + std::to_string(other.size()) + " points instead of " + std::to_string(reference.size())
                             + ".");
  for (std::size_t k = 0; k < reference.size(); k++) {
    const auto expected = reference[k].first;
    const auto found    = other[k].first;
    if (std::abs(found - expected) > 1e-12 * std::max(1.0, std::abs(expected)))
      throw std::runtime_error(name + ": expected the same frequency grid as " + reference_name + "; point "
                               + std::to_string(k + 1) + " is " + std::to_string(found) + " instead of "
                               + std::to_string(expected) + ".");
  }
}

// Which imaginary parts are tabulated. The off-diagonal ones are required either all or none; a diagonal one is
// optional and must vanish. 'complex' reports whether the calculation has to be done in complex arithmetic.
struct ImaginaryParts {
  std::vector<bool> present;
  bool complex{};
};

inline auto scan_imaginary_parts(const GammaOptions &options) {
  const auto n = options.channels;
  ImaginaryParts result;
  result.present.assign(static_cast<std::size_t>(n) * n, false);
  int offdiag_present = 0;
  std::string missing;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      const auto index      = static_cast<std::size_t>((i - 1) * n + (j - 1));
      const auto filename   = gamma_filename(options.prefix, i, j, true);
      const auto found      = std::filesystem::exists(filename);
      result.present[index] = found;
      if (i == j) continue; // a diagonal imaginary part must vanish, so it does not decide the arithmetic
      if (found)
        offdiag_present++;
      else if (missing.empty())
        missing = filename;
    }
  }
  const int offdiag_total = n * (n - 1);
  if (offdiag_present != 0 && offdiag_present != offdiag_total)
    throw std::runtime_error("The off-diagonal imaginary parts of Gamma must be given either all or none: "
                             + std::to_string(offdiag_present) + " of " + std::to_string(offdiag_total)
                             + " files are present, " + missing + " is missing.");
  result.complex = offdiag_present > 0;
  return result;
}

// The largest |Gamma_ij| in the input, used as the scale of the relative tolerance tests. Only the real parts are
// scanned: for a positive semidefinite Gamma the diagonal, which is real, dominates every element, and a spuriously
// large imaginary part must not inflate the tolerance that is meant to catch it.
inline auto magnitude(const std::vector<Pairs> &components) {
  double scale = 0.0;
  for (const auto &pairs : components)
    for (const auto &[omega, value] : pairs) scale = std::max(scale, std::abs(value));
  return scale;
}

// Gamma - Gamma^dagger measures both the mismatch of Gamma_ji against conj(Gamma_ij) and a nonzero imaginary part of
// a diagonal element, the latter as 2i Im Gamma_ii. For a real S it is the deviation from a symmetric matrix, and the
// tabulated imaginary parts are checked separately by the caller.
template<typename S>
void check_hermiticity(const std::vector<double> &omega, const std::vector<Matrix<S>> &gamma, const double scale,
                       const double tolerance) {
  double worst           = 0.0;
  std::size_t worst_node = 0;
  Eigen::Index worst_i   = 0;
  Eigen::Index worst_j   = 0;
  for (std::size_t k = 0; k < gamma.size(); k++) {
    const Matrix<double> deviation = (gamma[k] - gamma[k].adjoint()).cwiseAbs();
    Eigen::Index i                 = 0;
    Eigen::Index j                 = 0;
    const auto value               = deviation.maxCoeff(&i, &j);
    if (value > worst) {
      worst      = value;
      worst_node = k;
      worst_i    = i;
      worst_j    = j;
    }
  }
  if (worst > tolerance * scale)
    throw std::runtime_error("Gamma is not Hermitian: at omega=" + std::to_string(omega[worst_node]) + " element ("
                             + std::to_string(worst_i + 1) + "," + std::to_string(worst_j + 1) + ") deviates by "
                             + std::to_string(worst) + ", which exceeds " + std::to_string(tolerance)
                             + " times the magnitude " + std::to_string(scale) + " of Gamma.");
}

// Keep the nodes of one sign, indexed by |omega|. A node at exactly omega=0 belongs to neither branch, as in the sign
// predicate of load_rho() in adapt.
template<typename S>
auto select_branch(const std::vector<double> &omega, const std::vector<Matrix<S>> &gamma, const Sign sign) {
  GammaBranch<S> branch;
  for (std::size_t k = 0; k < omega.size(); k++) {
    const auto accept = sign == Sign::POS ? omega[k] > 0.0 : omega[k] < 0.0;
    if (!accept) continue;
    branch.omega.push_back(std::abs(omega[k]));
    branch.gamma.push_back(gamma[k]);
  }
  if (branch.omega.empty()) throw std::runtime_error("No data found at " + sign_name(sign) + " frequencies.");
  if (sign == Sign::NEG) { // |omega| must increase
    std::reverse(branch.omega.begin(), branch.omega.end());
    std::reverse(branch.gamma.begin(), branch.gamma.end());
  }
  branch.innermost = branch.omega.front();
  if (branch.omega.front() > zero_point) { // constant extrapolation to omega=0, as add_zero_point() does
    branch.omega.insert(branch.omega.begin(), zero_point);
    branch.gamma.insert(branch.gamma.begin(), branch.gamma.front());
  }
  return branch;
}

// The integral of one diagonal element of Gamma over the nodes with lower <= omega <= upper, taking the input
// tabulation as piecewise linear. Only used for the coverage report below.
template<typename S>
auto integrate_diagonal(const std::vector<double> &omega, const std::vector<Matrix<S>> &gamma, const int i,
                        const double lower, const double upper) {
  Vec table;
  table.reserve(omega.size());
  for (std::size_t k = 0; k < omega.size(); k++) table.emplace_back(omega[k], std::real(gamma[k](i, i)));
  return piecewise_linear_integral(table, lower, upper);
}

// The discretization mesh never reaches beyond |omega|=1, so any weight tabulated further out is discarded; and where
// the input stops short of the band edge, TabulatedDensity continues the density at its last tabulated value, which
// adds weight that is not in the input. Both are reported per diagonal element: for a positive semidefinite Gamma
// these are real and nonnegative, and |Gamma_ij|^2 <= Gamma_ii Gamma_jj bounds the off-diagonal elements by them.
template<typename S>
void report_band_coverage(const std::vector<double> &omega, const std::vector<Matrix<S>> &gamma, std::ostream &out) {
  const auto channels = static_cast<int>(gamma.front().rows());
  const auto lowest   = omega.front();
  const auto highest  = omega.back();
  for (int i = 0; i < channels; i++) {
    const auto total = integrate_diagonal(omega, gamma, i, lowest, highest);
    const auto name  = "# Gamma_" + std::to_string(i + 1) + std::to_string(i + 1) + ": ";
    const auto fraction = [total](const double weight) {
      return total > 0.0 ? " (" + std::to_string(100.0 * weight / total) + "%)" : std::string();
    };

    const auto discarded = integrate_diagonal(omega, gamma, i, 1.0, highest)
                           + integrate_diagonal(omega, gamma, i, lowest, -1.0);
    if (discarded > 0.0)
      out << name << discarded << " of " << total << fraction(discarded)
          << " of the weight lies beyond the band edge and is discarded" << std::endl;

    double extrapolated = 0.0;
    if (highest < 1.0) extrapolated += std::real(gamma.back()(i, i)) * (1.0 - highest);
    if (lowest > -1.0) extrapolated += std::real(gamma.front()(i, i)) * (1.0 + lowest);
    if (extrapolated > 0.0)
      out << name << extrapolated << " of " << total << fraction(extrapolated)
          << " of the weight is added by extrapolation to the band edge" << std::endl;
  }
}

template<typename S> void report_branch(const GammaBranch<S> &branch, const Sign sign, std::ostream &out) {
  out << "# Gamma - " << sign_name(sign) << " - " << branch.size() << " nodes - interval [ " << branch.omega.front()
      << " : " << branch.omega.back() << " ]" << std::endl;
}

} // namespace detail

// Whether the input requires complex arithmetic. Call this first and instantiate load_gamma() and everything
// downstream for the corresponding scalar type.
inline auto gamma_is_complex(const GammaOptions &options) {
  detail::validate_options(options);
  return detail::scan_imaginary_parts(options).complex;
}

// Read the channels^2 components of Gamma. All -re files are required. The off-diagonal -im files are required either
// all or none; none means that Gamma is real. A diagonal -im file is optional and must vanish within the Hermiticity
// tolerance, as must the deviation of Gamma_ji from conj(Gamma_ij); Gamma is then symmetrized to (Gamma+Gamma^dag)/2.
template<typename S> auto load_gamma(const GammaOptions &options, std::ostream &out = std::cout) {
  detail::validate_options(options);
  const auto n         = options.channels;
  const auto imaginary = detail::scan_imaginary_parts(options);
  if (imaginary.complex && !is_complex_v<S>)
    throw std::logic_error("The off-diagonal imaginary parts of Gamma are tabulated, so the discretization must be "
                           "performed in complex arithmetic.");

  const auto index = [n](const int i, const int j) { return static_cast<std::size_t>((i - 1) * n + (j - 1)); };

  std::vector<detail::Pairs> re(static_cast<std::size_t>(n) * n);
  std::vector<detail::Pairs> im(static_cast<std::size_t>(n) * n);
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      re[index(i, j)] = detail::read_column_file(gamma_filename(options.prefix, i, j, false));
      if (imaginary.present[index(i, j)])
        im[index(i, j)] = detail::read_column_file(gamma_filename(options.prefix, i, j, true));
    }
  }

  const auto reference_name = gamma_filename(options.prefix, 1, 1, false);
  const auto &grid          = re[index(1, 1)];
  detail::check_strictly_increasing(grid, reference_name);
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      if (!(i == 1 && j == 1))
        detail::check_same_grid(grid, reference_name, re[index(i, j)], gamma_filename(options.prefix, i, j, false));
      if (imaginary.present[index(i, j)])
        detail::check_same_grid(grid, reference_name, im[index(i, j)], gamma_filename(options.prefix, i, j, true));
    }
  }

  const auto scale = detail::magnitude(re);
  if (!(scale > 0.0)) throw std::runtime_error("Gamma vanishes identically.");

  // In real arithmetic the tabulated imaginary parts never reach the matrix, so they are checked here rather than
  // through the Hermiticity test. Only the diagonal ones can be present at this point.
  if constexpr (!is_complex_v<S>) {
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        if (!imaginary.present[index(i, j)]) continue;
        for (const auto &[frequency, value] : im[index(i, j)]) {
          if (std::abs(value) > options.hermiticity_tolerance * scale)
            throw std::runtime_error(gamma_filename(options.prefix, i, j, true)
                                     + ": the imaginary part of a diagonal element of Gamma must vanish; at omega="
                                     + std::to_string(frequency) + " it is " + std::to_string(value) + ".");
        }
      }
    }
  }

  const auto nodes = grid.size();
  std::vector<double> omega(nodes);
  std::vector<Matrix<S>> gamma(nodes, Matrix<S>::Zero(n, n));
  for (std::size_t k = 0; k < nodes; k++) {
    omega[k] = grid[k].first;
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        const auto real_part      = re[index(i, j)][k].second;
        const auto imaginary_part = imaginary.present[index(i, j)] ? im[index(i, j)][k].second : 0.0;
        gamma[k](i - 1, j - 1)    = make_scalar<S>(real_part, imaginary_part);
      }
    }
  }

  detail::check_hermiticity<S>(omega, gamma, scale, options.hermiticity_tolerance);
  for (auto &m : gamma) m = (0.5 * (m + m.adjoint())).eval();

  // omega -> omega/bandrescale and Gamma -> Gamma*bandrescale puts the band edge at 1 and leaves the integral of
  // Gamma unchanged. This is what rescalevecxy() does in adapt and nrgchain.
  for (std::size_t k = 0; k < nodes; k++) {
    omega[k] /= options.bandrescale;
    gamma[k] *= options.bandrescale;
  }

  out << "# Gamma: channels=" << n << " complex=" << (is_complex_v<S> ? 1 : 0) << " nodes=" << nodes << " interval [ "
      << omega.front() << " : " << omega.back() << " ]" << std::endl;
  detail::report_band_coverage<S>(omega, gamma, out);

  GammaInput<S> input;
  input.channels = n;
  input.pos      = detail::select_branch<S>(omega, gamma, Sign::POS);
  input.neg      = detail::select_branch<S>(omega, gamma, Sign::NEG);
  detail::report_branch(input.pos, Sign::POS, out);
  detail::report_branch(input.neg, Sign::NEG, out);
  return input;
}

} // namespace NRG::MixChain

#endif
