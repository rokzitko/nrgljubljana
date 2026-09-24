#ifndef NRG_STAR_TO_CHAIN_HPP
#define NRG_STAR_TO_CHAIN_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace NRG {

struct StarPoint {
  double energy;
  double amplitude; // Square root of the spectral weight; need not be normalized.
};

struct ScalarChain {
  std::vector<double> zeta;
  std::vector<double> xi;
  bool terminated;
};

// Gragg and Harrod, Numer. Math. 44 (1984), 317-335, equations (7)-(14).
// The unsquared rotation form avoids underflow of weights and squared hoppings.
// Preserve insertion order: Wilson stars should interleave positive/negative
// shells from high to low energy. Grouping branches can spoil tail accuracy.
// Return count onsite energies and hoppings, xi[n] coupling sites n and n+1.
// Only exact support exhaustion produces a terminal zero; no tolerance deflation.
inline ScalarChain scalar_star_to_chain(std::span<const StarPoint> star, const std::size_t count) {
  if (count == 0) throw std::invalid_argument("rkpw: coefficient count must be positive.");
  double largest_amplitude = 0.0;
  for (const auto &point : star) {
    if (!std::isfinite(point.energy) || !std::isfinite(point.amplitude) || point.amplitude < 0.0)
      throw std::invalid_argument("rkpw: energies must be finite and amplitudes finite and nonnegative.");
    largest_amplitude = std::max(largest_amplitude, point.amplitude);
  }
  int exponent = 0;
  std::frexp(largest_amplitude, &exponent);
  std::vector<StarPoint> support;
  std::map<double, std::size_t> positions;
  for (const auto &point : star) {
    if (point.amplitude == 0.0) continue;
    // A common power-of-two scale is immaterial to the measure. Normalize BEFORE
    // merging/norms: hypot(denorm_min,denorm_min) loses the rotation's direction.
    const auto scaled = std::scalbn(point.amplitude, 1 - exponent);
    if (scaled == 0.0) throw std::runtime_error("rkpw: amplitude normalization underflow.");
    const auto [it, inserted] = positions.emplace(point.energy, support.size());
    if (inserted) {
      support.push_back({point.energy, scaled});
    } else {
      auto &amplitude = support[it->second].amplitude;
      amplitude = std::hypot(amplitude, scaled);
      if (!std::isfinite(amplitude)) throw std::runtime_error("rkpw: merged amplitude is not representable.");
    }
  }
  const auto rank = support.size();
  if (count > rank) throw std::invalid_argument("rkpw: requested coefficient count exceeds effective nonzero support.");
  const auto cap = count < rank ? count + 1 : rank;
  std::vector<double> a(cap, 0.0), b(cap, 0.0);
  a[0] = support[0].energy;
  b[0] = support[0].amplitude; // b[0] carries the total amplitude, not xi[0].
  for (std::size_t i = 1; i < rank; ++i) {
    const auto node = support[i].energy;
    if (i < cap) a[i] = node;
    double pi = support[i].amplitude, cold = 1.0, sold = 0.0, tauold = 0.0;
    // All poles contribute even when only a short leading chain is requested.
    for (std::size_t k = 0; k < std::min(i + 1, cap); ++k) {
      const auto oldb = b[k];
      const auto rho = std::hypot(oldb, pi);
      const auto nextb = cold * rho;
      auto c = rho == 0.0 ? 1.0 : oldb / rho;
      auto s = rho == 0.0 ? 0.0 : -pi / rho;
      if (rho > 0.0 && rho < std::numeric_limits<double>::min()) {
        // The rounded subnormal norm may have very few bits. Recover the
        // rotation direction from locally scaled operands, not that norm.
        const auto scale = std::max(std::abs(oldb), std::abs(pi));
        const auto x = oldb / scale, y = pi / scale;
        const auto norm = std::hypot(x, y);
        c = x / norm;
        s = -y / norm;
      }
      const auto nextpi = s * (a[k] - node) + c * sold * oldb;
      const auto nexttau = s * nextpi;
      a[k] -= nexttau - tauold;
      b[k] = nextb;
      if (!std::isfinite(a[k]) || !std::isfinite(nextb) || !std::isfinite(nextpi))
        throw std::runtime_error("rkpw: numerical breakdown or nonrepresentable intermediate coefficient.");
      pi = nextpi;
      cold = c;
      sold = s;
      tauold = nexttau;
    }
  }

  ScalarChain result{{a.begin(), a.begin() + count}, {b.begin() + 1, b.end()}, count == rank};
  for (const auto hopping : result.xi)
    if (!(hopping > 0.0)) throw std::runtime_error("rkpw: nonpositive or underflowed hopping before support exhaustion.");
  if (result.terminated) result.xi.push_back(0.0);
  return result;
}

} // namespace NRG

#endif
