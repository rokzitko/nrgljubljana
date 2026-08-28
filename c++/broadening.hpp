// Broadening kernels shared by the NRG runtime and the broaden tool.
// Copyright (C) 2005-2026 Rok Zitko

#ifndef _broadening_hpp_
#define _broadening_hpp_

#include <cassert>
#include <cmath>

namespace NRG::Broadening {

inline constexpr double sqrt_pi = 1.7724538509055160273;

enum class CrossoverMode { output_frequency, peak_frequency };

// Modified log-Gaussian kernel from A. Weichselbaum and J. von Delft,
// Phys. Rev. Lett. 99, 076402 (2007), with gamma=alpha/4.
inline double log_gaussian(const double output, const double peak, const double alpha) {
  if ((output < 0.0 && peak > 0.0) || (output > 0.0 && peak < 0.0)) return 0.0;
  if (peak == 0.0) return 0.0;
  const double gamma = alpha / 4.0;
  return std::exp(-std::pow(std::log(output / peak) / alpha - gamma, 2))
         / (alpha * std::abs(output) * sqrt_pi);
}

// Shift the logarithmic accumulation points to +/-accumulation. The Gaussian
// and crossover terms continue to use the original frequencies.
inline double log_gaussian(const double output, const double peak, const double alpha,
                           const double accumulation) {
  auto shifted_output = output;
  auto shifted_peak = peak;
  if (shifted_output > accumulation && shifted_peak > accumulation) {
    shifted_output -= accumulation;
    shifted_peak -= accumulation;
  }
  if (shifted_output < -accumulation && shifted_peak < -accumulation) {
    const auto shift = -accumulation;
    shifted_output -= shift;
    shifted_peak -= shift;
  }
  return log_gaussian(shifted_output, shifted_peak, alpha);
}

// Unit-area Gaussian with the width convention used by the primary NRG
// broadening kernel. Its standard deviation is width/sqrt(2).
inline double gaussian(const double output, const double peak, const double width) {
  return std::exp(-std::pow((output - peak) / width, 2)) / (width * sqrt_pi);
}

inline double crossover(const double frequency, const double alpha, const double omega0) {
  const auto absolute_frequency = std::abs(frequency);
  return absolute_frequency > omega0
           ? 1.0
           : std::exp(-std::pow(std::log(absolute_frequency / omega0) / alpha, 2));
}

inline double crossover(const double output, const double peak, const double alpha,
                        const double omega0, const CrossoverMode mode) {
  const auto frequency = mode == CrossoverMode::peak_frequency ? peak : output;
  return crossover(frequency, alpha, omega0);
}

inline double blend(const double log_part, const double gaussian_part, const double h) {
  return log_part * h + gaussian_part * (1.0 - h);
}

inline double hybrid_kernel(const double output, const double peak, const double alpha,
                            const double omega0, const CrossoverMode mode,
                            const double accumulation = 0.0) {
  const auto log_part = log_gaussian(output, peak, alpha, accumulation);
  const auto gaussian_part = gaussian(output, peak, omega0);
  const auto h = crossover(output, peak, alpha, omega0, mode);
  assert(h >= 0.0 && h <= 1.0);
  return blend(log_part, gaussian_part, h);
}

} // namespace NRG::Broadening

#endif
