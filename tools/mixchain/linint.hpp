// Channel-mixing discretization for NRG
// ** Interpolation code

#ifndef _mixchain_linint_hpp_
#define _mixchain_linint_hpp_

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../common/linint.hpp"

namespace NRG::MixChain {

// Structures for storing tabulated data, such as one element of Gamma(omega) or the weight function of the mesh.
using Pair = std::pair<double, double>;
using Vec  = std::vector<Pair>;

struct ThrowError {
  [[noreturn]] void operator()(const std::string &message) const { throw std::runtime_error(message); }
};

using LinInt = NRG::Tools::LinIntBase<Vec, ThrowError>;

// The exact integral of a piecewise-linear table over [lower,upper], restricted to the tabulated range. The values
// may have either sign, so this is not a job for TabulatedDensity.
inline auto piecewise_linear_integral(const Vec &table, const double lower, const double upper) {
  double sum = 0.0;
  for (std::size_t k = 1; k < table.size(); k++) {
    const auto [x0, y0] = table[k - 1];
    const auto [x1, y1] = table[k];
    const auto left     = std::max(x0, lower);
    const auto right    = std::min(x1, upper);
    if (right <= left) continue;
    const auto slope = (y1 - y0) / (x1 - x0);
    // Captured by name: structured bindings cannot be captured implicitly before Clang 16.
    const auto at    = [x0 = x0, y0 = y0, slope](const double x) { return y0 + slope * (x - x0); };
    sum += (right - left) * (at(left) + at(right)) / 2.0;
  }
  return sum;
}

} // namespace NRG::MixChain

#endif
