#ifndef _tools_common_cumulative_weight_hpp_
#define _tools_common_cumulative_weight_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tabulated_density.hpp"

namespace NRG::Tools {

// Cumulative hybridisation weight on one frequency branch, normalized to its value at omega=1:
//
//   W(omega) = int_0^omega rho(u) du / int_0^1 rho(u) du,
//
// together with its generalized inverse. Where rho vanishes over an interval, W has a plateau and is not invertible;
// the upper edge of the plateau is returned in that case.
//
// The density is held by reference: the caller owns the TabulatedDensity. Nothing here is const, because evaluating
// the density updates its internal caches.
class CumulativeWeight {
 private:
  TabulatedDensity *rho_{};
  double total_{};
  std::vector<std::pair<double, double>> plateaus_; // (weight at the plateau, upper edge)

 public:
  CumulativeWeight() = default;

  // 'samples' are the raw tabulated (omega, rho) pairs backing the density; they are scanned for zero-density
  // plateaus. The density and the samples must describe the same data.
  template<typename VecT>
  CumulativeWeight(TabulatedDensity &rho, const VecT &samples) : rho_(&rho) {
    total_  = rho_->integral(0.0, 1.0);
    if (!(std::isfinite(total_) && total_ > 0.0)) {
      throw std::runtime_error("Integral method requires positive spectral weight in [0,1].");
    }
    plateaus_.clear();
    for (std::size_t i = 0; i + 1 < samples.size();) {
      if (samples[i].second == 0.0 && samples[i + 1].second == 0.0) {
        std::size_t last = i + 1;
        while (last + 1 < samples.size() && samples[last + 1].second == 0.0) { last++; }
        const double lower = std::max(0.0, samples[i].first);
        double upper = std::min(1.0, samples[last].first);
        if (last + 1 == samples.size() && upper < 1.0) upper = 1.0;
        if (lower < upper) {
          const double midpoint = lower + (upper - lower) / 2.0;
          plateaus_.emplace_back(normalized(midpoint), upper);
        }
        i = last;
      } else {
        i++;
      }
    }
  }

  [[nodiscard]] auto total() const { return total_; }

  double normalized(const double omega) {
    return rho_->integral(0.0, omega) / total_;
  }

  // Generalized inverse of W. For a zero-density plateau, return its upper edge.
  auto inverse(const double weight) {
    for (const auto &[plateau_weight, upper_edge] : plateaus_) {
      if (weight == plateau_weight) return upper_edge;
    }
    double lower = 0.0;
    double upper = 1.0;
    while (true) {
      const double midpoint = lower + (upper - lower) / 2.0;
      if (midpoint == lower || midpoint == upper) break;
      if (normalized(midpoint) <= weight) {
        lower = midpoint;
      } else {
        upper = midpoint;
      }
    }
    return lower + (upper - lower) / 2.0;
  }
};

} // namespace NRG::Tools

#endif
