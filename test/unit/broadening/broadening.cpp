#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include <broadening.hpp>

namespace {

template<typename Function>
double midpoint_integral(Function function, const double lower, const double upper,
                         const std::size_t intervals) {
  const auto width = (upper - lower) / static_cast<double>(intervals);
  auto sum = 0.0;
  for (std::size_t i = 0; i < intervals; ++i) {
    const auto frequency = lower + (static_cast<double>(i) + 0.5) * width;
    sum += function(frequency);
  }
  return sum * width;
}

} // namespace

TEST(Broadening, base_kernels) { // NOLINT
  using namespace NRG::Broadening;

  constexpr double alpha = 0.3;
  constexpr double omega0 = 0.1;
  constexpr double peak = 0.2;

  EXPECT_EQ(log_gaussian(-0.1, peak, alpha), 0.0);
  EXPECT_EQ(log_gaussian(0.1, 0.0, alpha), 0.0);
  EXPECT_DOUBLE_EQ(gaussian(peak, peak, omega0), 1.0 / (omega0 * sqrt_pi));
  EXPECT_DOUBLE_EQ(crossover(omega0, alpha, omega0), 1.0);
  EXPECT_DOUBLE_EQ(crossover(2.0 * omega0, alpha, omega0), 1.0);
}

TEST(Broadening, crossover_mode_selects_frequency) { // NOLINT
  using namespace NRG::Broadening;

  constexpr double output = 0.03;
  constexpr double peak = 0.2;
  constexpr double alpha = 0.3;
  constexpr double omega0 = 0.1;

  const auto log_part = log_gaussian(output, peak, alpha);
  const auto gaussian_part = gaussian(output, peak, omega0);
  const auto output_h = crossover(output, alpha, omega0);
  const auto peak_h = crossover(peak, alpha, omega0);

  EXPECT_DOUBLE_EQ(hybrid_kernel(output, peak, alpha, omega0, CrossoverMode::output_frequency),
                   blend(log_part, gaussian_part, output_h));
  EXPECT_DOUBLE_EQ(hybrid_kernel(output, peak, alpha, omega0, CrossoverMode::peak_frequency),
                   blend(log_part, gaussian_part, peak_h));
  EXPECT_NE(output_h, peak_h);
}

TEST(Broadening, shifted_accumulation_coordinates) { // NOLINT
  using namespace NRG::Broadening;

  constexpr double alpha = 0.3;
  constexpr double accumulation = 0.2;

  EXPECT_DOUBLE_EQ(log_gaussian(0.5, 0.4, alpha, accumulation),
                   log_gaussian(0.3, 0.2, alpha));
  EXPECT_DOUBLE_EQ(log_gaussian(-0.5, -0.4, alpha, accumulation),
                   log_gaussian(-0.3, -0.2, alpha));
  EXPECT_DOUBLE_EQ(log_gaussian(0.1, 0.4, alpha, accumulation),
                   log_gaussian(0.1, 0.4, alpha));
}

TEST(Broadening, peak_frequency_mode_conserves_full_axis_weight) { // NOLINT
  using namespace NRG::Broadening;

  constexpr double peak = 0.05;
  constexpr double alpha = 0.3;
  constexpr double omega0 = 0.1;
  constexpr std::size_t intervals = 200000;

  const auto normalized = midpoint_integral(
    [](const double output) {
      return hybrid_kernel(output, peak, alpha, omega0, CrossoverMode::peak_frequency);
    },
    -2.0, 2.0, intervals);
  const auto output_weighted = midpoint_integral(
    [](const double output) {
      return hybrid_kernel(output, peak, alpha, omega0, CrossoverMode::output_frequency);
    },
    -2.0, 2.0, intervals);

  EXPECT_NEAR(normalized, 1.0, 1e-10);
  EXPECT_GT(std::abs(output_weighted - 1.0), 0.1);
}
