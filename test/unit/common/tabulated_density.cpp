#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <common/gsl_piecewise_polynomial.hpp>
#include <common/tabulated_density.hpp>

namespace {

auto knots(const NRG::Tools::DensityTable &samples) {
  std::vector<double> result;
  for (const auto &[energy, value] : samples) {
    (void)value;
    result.push_back(energy);
  }
  return result;
}

auto values(const NRG::Tools::DensityTable &samples) {
  std::vector<double> result;
  for (const auto &[energy, value] : samples) {
    (void)energy;
    result.push_back(value);
  }
  return result;
}

} // namespace

TEST(TabulatedDensity, parses_only_shape_preserving_density_methods) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  EXPECT_EQ(NRG::Tools::parse_density_interpolation_method("linear"), InterpolationMethod::linear);
  EXPECT_EQ(NRG::Tools::parse_density_interpolation_method("steffen"), InterpolationMethod::steffen);
  EXPECT_THROW(NRG::Tools::parse_density_interpolation_method("akima"), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::parse_density_interpolation_method("unknown"), std::invalid_argument);
}

TEST(TabulatedDensity, linear_values_and_primitive_include_constant_tails) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  NRG::Tools::TabulatedDensity density({{0.25, 2.0}, {0.5, 1.0}, {1.0, 3.0}}, InterpolationMethod::linear);

  EXPECT_DOUBLE_EQ(density(0.0), 2.0);
  EXPECT_DOUBLE_EQ(density(0.375), 1.5);
  EXPECT_DOUBLE_EQ(density(2.0), 3.0);
  EXPECT_DOUBLE_EQ(density.cumulative(0.0), -0.5);
  EXPECT_DOUBLE_EQ(density.cumulative(0.5), 0.375);
  EXPECT_DOUBLE_EQ(density.integral(0.0, 0.25), 0.5);
  EXPECT_DOUBLE_EQ(density.integral(0.25, 0.5), 0.375);
  EXPECT_DOUBLE_EQ(density.integral(1.0, 2.0), 3.0);
  EXPECT_DOUBLE_EQ(density.integral(0.0, 2.0), 4.875);
}

TEST(TabulatedDensity, direct_local_integrals_preserve_tiny_weights) { // NOLINT
  NRG::Tools::TabulatedDensity density({{1e-99, 1.0}, {1.0, 1.0}},
                                       NRG::Tools::InterpolationMethod::linear);

  EXPECT_DOUBLE_EQ(density.integral(1e-120, 2e-120), 1e-120);
  constexpr double cancellation_bound = 1.532107773982716e-113;
  EXPECT_DOUBLE_EQ(density.integral(0.0, cancellation_bound), cancellation_bound);
  const auto lower = std::nextafter(1.0, 0.0);
  EXPECT_DOUBLE_EQ(density.integral(lower, 1.0), 1.0 - lower);

  NRG::Tools::TabulatedDensity interior(
    {{0.0, 1.0}, {1.0, 1.0}, {2.0, 0.0}, {3.0, 1e-300}, {4.0, 1e-300}, {5.0, 0.0}},
    NRG::Tools::InterpolationMethod::linear);
  EXPECT_DOUBLE_EQ(interior.integral(2.0, 5.0), 2e-300);
}

TEST(TabulatedDensity, steffen_matches_the_materialized_polynomial_and_preserves_range) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const NRG::Tools::DensityTable samples{{0.1, 0.2}, {0.25, 2.0}, {0.6, 0.4}, {0.8, 1.5}, {1.0, 0.7}};
  NRG::Tools::TabulatedDensity density(samples, InterpolationMethod::steffen);
  const auto polynomial = NRG::Tools::make_gsl_piecewise_polynomial(knots(samples), values(samples),
                                                                     InterpolationMethod::steffen);

  for (std::size_t interval = 0; interval + 1 < samples.size(); ++interval) {
    for (const double fraction : {0.0, 0.13, 0.5, 0.91}) {
      const auto x = samples[interval].first
                     + fraction * (samples[interval + 1].first - samples[interval].first);
      const auto interpolated = density(x);
      EXPECT_NEAR(interpolated, polynomial.evaluate(x), 2e-14);
      EXPECT_GE(interpolated, std::min(samples[interval].second, samples[interval + 1].second));
      EXPECT_LE(interpolated, std::max(samples[interval].second, samples[interval + 1].second));
    }
  }

  EXPECT_NEAR(density.integral(0.1, 1.0), polynomial.integral(), 2e-14);
  EXPECT_NEAR(density.integral(0.0, 1.2), 0.1 * 0.2 + polynomial.integral() + 0.2 * 0.7, 2e-14);
  EXPECT_NEAR(density.cumulative(0.8) - density.cumulative(0.25), polynomial.integral(0.25, 0.8), 2e-14);

  NRG::Tools::TabulatedDensity linear(samples, InterpolationMethod::linear);
  EXPECT_GT(std::abs(density(0.5) - linear(0.5)), 1e-3);
}

TEST(TabulatedDensity, steffen_preserves_zero_plateaus) { // NOLINT
  NRG::Tools::TabulatedDensity density({{0.0, 1.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 1.0}},
                                       NRG::Tools::InterpolationMethod::steffen);

  EXPECT_DOUBLE_EQ(density(1.5), 0.0);
  EXPECT_DOUBLE_EQ(density.integral(1.0, 2.0), 0.0);
  EXPECT_DOUBLE_EQ(density.cumulative(1.0), density.cumulative(2.0));
}

TEST(TabulatedDensity, steffen_tolerates_roundoff_near_zero_knots) { // NOLINT
  NRG::Tools::TabulatedDensity density(
    {{0.0, 0.0062727669}, {1.0, 0.2249509465}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}},
    NRG::Tools::InterpolationMethod::steffen);
  const auto x = 2.0 - std::ldexp(1.0, -30);

  EXPECT_NO_THROW((void)density(x));
  EXPECT_GE(density(x), 0.0);
}

TEST(TabulatedDensity, cumulative_evaluation_does_not_change_point_interval_flag) { // NOLINT
  NRG::Tools::TabulatedDensity density({{0.0, 1.0}, {1.0, 2.0}, {2.0, 1.0}, {3.0, 2.0}},
                                       NRG::Tools::InterpolationMethod::steffen);

  (void)density(0.25);
  EXPECT_TRUE(density.flag());
  density.clear_flag();
  (void)density(0.75);
  EXPECT_FALSE(density.flag());
  (void)density.cumulative(2.5);
  EXPECT_FALSE(density.flag());
  (void)density(2.5);
  EXPECT_TRUE(density.flag());
}

TEST(TabulatedDensity, guards_common_input_mistakes) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  using NRG::Tools::TabulatedDensity;
  const auto infinity = std::numeric_limits<double>::infinity();

  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}}, InterpolationMethod::linear), std::invalid_argument);
  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}, {1.0, 1.0}}, InterpolationMethod::steffen), std::invalid_argument);
  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}}, InterpolationMethod::steffen),
               std::invalid_argument);
  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}, {1.0, -1.0}}, InterpolationMethod::linear), std::invalid_argument);
  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}, {infinity, 1.0}}, InterpolationMethod::linear),
               std::invalid_argument);
  EXPECT_THROW(TabulatedDensity({{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}}, InterpolationMethod::cspline),
               std::invalid_argument);

  TabulatedDensity uninitialized;
  EXPECT_THROW((void)uninitialized(0.5), std::logic_error);
  TabulatedDensity density({{0.0, 1.0}, {1.0, 1.0}}, InterpolationMethod::linear);
  EXPECT_THROW(density.integral(1.0, 0.0), std::invalid_argument);
  EXPECT_THROW((void)density(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

TEST(TabulatedDensity, normalized_amplitude_preserves_normal_mass_arithmetic_exactly) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity density({{0.25, 2.0}, {0.5, 1.0}, {1.0, 3.0}}, method);
    for (const auto &[lower, upper] : {std::pair{0.0, 2.0}, std::pair{0.3, 0.8}, std::pair{1.0, 2.0}}) {
      const auto mass = density.integral(lower, upper);
      for (const auto total : {1e-200, 1.0, 4.875, 1e200})
        EXPECT_EQ(density.normalized_amplitude(lower, upper, std::sqrt(total)), std::sqrt(mass) / std::sqrt(total));
    }
    NRG::Tools::TabulatedDensity constant({{0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}}, method);
    const auto minimum = std::numeric_limits<double>::min();
    EXPECT_EQ(constant.normalized_amplitude(0.0, minimum, std::sqrt(3.0)), std::sqrt(minimum) / std::sqrt(3.0));
  }
}

TEST(TabulatedDensity, normalized_amplitude_is_invariant_under_uniform_density_scaling) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity ordinary({{0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}}, method);
    NRG::Tools::TabulatedDensity tiny({{0.0, 1e-200}, {0.5, 1e-200}, {1.0, 1e-200}}, method);
    for (const auto width : {0.25, 1e-120, 1e-150}) {
      const auto expected = ordinary.normalized_amplitude(0.0, width, 1.0);
      const auto actual = tiny.normalized_amplitude(0.0, width, std::sqrt(1e-200));
      EXPECT_GT(actual, 0.0);
      EXPECT_NEAR(actual / expected, 1.0, 2e-15);
    }
    EXPECT_EQ(tiny.integral(0.0, 1e-150), 0.0);
    EXPECT_NEAR(tiny.normalized_amplitude(0.0, 1e-150, std::sqrt(1e-200)) / 1e-75, 1.0, 2e-15);
  }
}

TEST(TabulatedDensity, normalized_amplitude_recovers_linear_pseudogaps_in_both_directions) { // NOLINT
  NRG::Tools::TabulatedDensity rising({{0.0, 0.0}, {1.0, 1.0}});
  NRG::Tools::TabulatedDensity falling({{-1.0, 1.0}, {0.0, 0.0}});
  const auto expected = std::sqrt(1.5) * 1e-200;
  EXPECT_EQ(rising.integral(1e-200, 2e-200), 0.0);
  EXPECT_EQ(falling.integral(-2e-200, -1e-200), 0.0);
  EXPECT_NEAR(rising.normalized_amplitude(1e-200, 2e-200, 1.0) / expected, 1.0, 2e-15);
  EXPECT_NEAR(falling.normalized_amplitude(-2e-200, -1e-200, 1.0) / expected, 1.0, 2e-15);
}

TEST(TabulatedDensity, normalized_amplitude_does_not_use_rounded_subnormal_mass) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity density({{0.0, 1e-200}, {0.5, 1e-200}, {1.0, 1e-200}}, method);
    const auto mass = density.integral(0.0, 1e-120);
    ASSERT_GT(mass, 0.0);
    ASSERT_LT(mass, std::numeric_limits<double>::min());
    const auto rounded = std::sqrt(mass) / std::sqrt(1e-200);
    EXPECT_GT(std::abs(rounded / 1e-60 - 1.0), 1e-7);
    EXPECT_NEAR(density.normalized_amplitude(0.0, 1e-120, std::sqrt(1e-200)) / 1e-60, 1.0, 2e-15);

    const auto minimum = std::numeric_limits<double>::denorm_min();
    NRG::Tools::TabulatedDensity subnormal({{0.0, minimum}, {1.0, minimum}, {2.0, minimum}}, method);
    EXPECT_EQ(subnormal.normalized_amplitude(0.0, 1.0, std::sqrt(minimum)), 1.0);
    EXPECT_NEAR(subnormal.normalized_amplitude(0.0, 0.5, std::sqrt(minimum)), std::sqrt(0.5), 2e-15);
  }
}

TEST(TabulatedDensity, normalized_amplitude_accumulates_multiple_pieces_and_constant_tails) { // NOLINT
  constexpr double width = 1e-150;
  constexpr double scale = 1e-200;
  const auto root_scale = std::sqrt(scale);
  NRG::Tools::TabulatedDensity density(
    {{width, scale}, {2 * width, 2 * scale}, {3 * width, 0.0}, {4 * width, 0.0}, {5 * width, 3 * scale}});
  EXPECT_EQ(density.integral(0.0, 6 * width), 0.0);
  const auto whole = density.normalized_amplitude(0.0, 6 * width, root_scale);
  EXPECT_NEAR(whole / (std::sqrt(8.0) * std::sqrt(width)), 1.0, 2e-15);
  const auto left = density.normalized_amplitude(0.0, 2.5 * width, root_scale);
  const auto right = density.normalized_amplitude(2.5 * width, 6 * width, root_scale);
  EXPECT_NEAR(std::hypot(left, right) / whole, 1.0, 2e-15);
  EXPECT_EQ(density.normalized_amplitude(3 * width, 4 * width, root_scale), 0.0);

  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity constant({{width, scale}, {2 * width, scale}, {3 * width, scale}}, method);
    EXPECT_NEAR(constant.normalized_amplitude(0.0, 4 * width, root_scale) / (2 * std::sqrt(width)), 1.0, 2e-15);
    EXPECT_NEAR(constant.normalized_amplitude(-width, 0.0, root_scale) / std::sqrt(width), 1.0, 2e-15);
    EXPECT_NEAR(constant.normalized_amplitude(4 * width, 5 * width, root_scale) / std::sqrt(width), 1.0, 2e-15);
  }
}

TEST(TabulatedDensity, normalized_amplitude_is_partition_independent_when_individual_pieces_underflow) { // NOLINT
  constexpr double width = 1e-46;
  constexpr double rho = 1e-300;
  const auto root_total_weight = std::sqrt(1e300);
  NRG::Tools::DensityTable samples;
  for (int i = 0; i <= 100; ++i) samples.emplace_back(width * (i / 100.0), rho);
  samples.emplace_back(2 * width, rho);
  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity split(samples, method);
    NRG::Tools::TabulatedDensity unsplit({{0.0, rho}, {width, rho}, {2 * width, rho}}, method);
    const auto expected = unsplit.normalized_amplitude(0.0, width, root_total_weight);
    ASSERT_EQ(expected, 2 * std::numeric_limits<double>::denorm_min());
    EXPECT_EQ(split.normalized_amplitude(0.0, width, root_total_weight), expected);
    for (int i = 0; i < 100; ++i)
      EXPECT_THROW(split.normalized_amplitude(samples[i].first, samples[i + 1].first, root_total_weight), std::runtime_error);
  }
}

TEST(TabulatedDensity, normalized_amplitude_preserves_root_normalization_without_rounding_total_mass) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const auto base_root = std::sqrt(std::numeric_limits<double>::min());
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    for (const auto width : {1e-23, 1e-24}) {
      NRG::Tools::TabulatedDensity density({{0.0, 1e-300}, {width, 1e-300}, {2 * width, 1e-300}}, method);
      const auto root_total_weight = std::sqrt(1e-300) * std::sqrt(width);
      const auto recovered_root = density.normalized_amplitude(0.0, width, base_root) * base_root;
      EXPECT_NEAR(recovered_root / root_total_weight, 1.0, 2e-15);
      for (const auto fraction : {0.25, 0.5, 1.0}) {
        EXPECT_NEAR(density.normalized_amplitude(0.0, fraction * width, root_total_weight), std::sqrt(fraction), 2e-15);
        EXPECT_NEAR(density.normalized_amplitude(0.0, fraction * width, recovered_root), std::sqrt(fraction), 2e-15);
      }
    }
  }
}

TEST(TabulatedDensity, normalized_amplitude_applies_root_normalization_only_at_the_end) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const auto minimum = std::numeric_limits<double>::denorm_min();
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity density({{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}}, method);
    // Dividing sqrt(mean) by this root first would overflow, although the complete result is finite.
    EXPECT_NEAR(density.normalized_amplitude(0.0, minimum, minimum) / (1.0 / std::sqrt(minimum)), 1.0, 2e-15);
    EXPECT_EQ(density.normalized_amplitude(0.0, 1.0, 1e200), 1.0 / 1e200);
  }
}

TEST(TabulatedDensity, normalized_amplitude_steffen_matches_analytic_polynomial) { // NOLINT
  // On the first segment Steffen gives rho/scale = 1 + u + u^2 - u^3.
  for (const auto &[width, scale] : {std::pair{1.0, 1.0}, std::pair{1e-150, 1e-200}}) {
    NRG::Tools::TabulatedDensity density({{0.0, scale}, {width, 2 * scale}, {2 * width, scale}},
                                       NRG::Tools::InterpolationMethod::steffen);
    const auto expected = std::sqrt(155.0 / 192.0) * std::sqrt(width);
    EXPECT_NEAR(density.normalized_amplitude(0.25 * width, 0.75 * width, std::sqrt(scale)) / expected, 1.0, 3e-15);
    EXPECT_NEAR(density.normalized_amplitude(0.0, 2 * width, std::sqrt(scale))
                  / (std::sqrt(19.0 / 6.0) * std::sqrt(width)), 1.0, 3e-15);
  }
}

TEST(TabulatedDensity, normalized_amplitude_returns_zero_only_for_empty_intervals_and_zero_density) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity density({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}, {3.0, 0.0}}, method);
    EXPECT_EQ(density.normalized_amplitude(1.0, 1.0, 1.0), 0.0);
    EXPECT_EQ(density.normalized_amplitude(-2.0, -1.0, 1.0), 0.0);
    EXPECT_EQ(density.normalized_amplitude(4.0, 5.0, 1.0), 0.0);
    EXPECT_EQ(density.normalized_amplitude(2.0, 3.0, 1.0), 0.0);
    EXPECT_EQ(density.normalized_amplitude(2.25, 2.75, 1.0), 0.0);
    EXPECT_EQ(density.normalized_amplitude(2.0, 5.0, 1.0), 0.0);
  }
}

TEST(TabulatedDensity, normalized_amplitude_diagnoses_underresolved_nonzero_steffen_shells) { // NOLINT
  NRG::Tools::TabulatedDensity density({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 0.0}},
                                     NRG::Tools::InterpolationMethod::steffen);
  try {
    (void)density.normalized_amplitude(-1e-100, 0.0, 1.0);
    FAIL() << "A positive overlap must not silently become a zero shell.";
  } catch (const std::runtime_error &error) {
    EXPECT_NE(std::string(error.what()).find("Underresolved nonzero density shell"), std::string::npos);
  }
}

TEST(TabulatedDensity, normalized_amplitude_diagnoses_unrepresentable_positive_output) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const auto minimum = std::numeric_limits<double>::denorm_min();
  const auto maximum = std::numeric_limits<double>::max();
  for (const auto method : {InterpolationMethod::linear, InterpolationMethod::steffen}) {
    NRG::Tools::TabulatedDensity density({{0.0, minimum}, {1.0, minimum}, {2.0, minimum}}, method);
    EXPECT_EQ(density.normalized_amplitude(-minimum, 0.0, 1.0), minimum);
    const auto expected = std::sqrt(minimum) / std::sqrt(maximum);
    ASSERT_GT(expected, 0.0);
    EXPECT_EQ(density.normalized_amplitude(-minimum, 1.0, std::sqrt(maximum)), expected);
    EXPECT_EQ(density.normalized_amplitude(-1.0, minimum, std::sqrt(maximum)), expected);
    try {
      (void)density.normalized_amplitude(-minimum, 0.0, std::sqrt(maximum));
      FAIL() << "A complete positive shell that narrows to zero must fail.";
    } catch (const std::runtime_error &error) {
      EXPECT_NE(std::string(error.what()).find("amplitude underflows double"), std::string::npos);
    }
  }
  NRG::Tools::TabulatedDensity large({{0.0, maximum}, {0.25, maximum}});
  EXPECT_THROW(large.normalized_amplitude(0.0, 0.25, std::sqrt(minimum)), std::runtime_error);
  EXPECT_THROW(large.normalized_amplitude(0.0, 2.0, 1.0), std::runtime_error);
}

TEST(TabulatedDensity, normalized_amplitude_validates_bounds_and_root_total_weight) { // NOLINT
  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  NRG::Tools::TabulatedDensity uninitialized;
  EXPECT_THROW(uninitialized.normalized_amplitude(0.0, 1.0, 1.0), std::logic_error);
  NRG::Tools::TabulatedDensity density({{0.0, 1.0}, {1.0, 1.0}});
  EXPECT_THROW(density.normalized_amplitude(1.0, 0.0, 1.0), std::invalid_argument);
  for (const auto invalid : {-infinity, infinity, nan}) {
    EXPECT_THROW(density.normalized_amplitude(invalid, 1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(density.normalized_amplitude(0.0, invalid, 1.0), std::invalid_argument);
  }
  for (const auto invalid : {-1.0, -0.0, 0.0, infinity, -infinity, nan}) {
    EXPECT_THROW(density.normalized_amplitude(0.0, 1.0, invalid), std::invalid_argument);
    EXPECT_THROW(density.normalized_amplitude(0.0, 0.0, invalid), std::invalid_argument);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
