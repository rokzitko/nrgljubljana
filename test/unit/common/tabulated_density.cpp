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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
