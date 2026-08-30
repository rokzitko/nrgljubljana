#include <gtest/gtest.h>

#include <limits>

#include <integ/integ.hpp>

namespace {

auto linear_identity() {
  return NRG::Integ::Polynomial{{-1.0, 1.0}, {{-1.0, 2.0}}};
}

} // namespace

TEST(Integ, exact_quantities_share_the_piecewise_polynomial) { // NOLINT
  const auto polynomial = linear_identity();
  NRG::Integ::Options options;

  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.0);

  options.quantity = NRG::Integ::Quantity::positive;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.5);

  options.quantity = NRG::Integ::Quantity::negative;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), -0.5);

  options.quantity = NRG::Integ::Quantity::absolute;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 1.0);

  options.quantity = NRG::Integ::Quantity::negative_absolute;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.5);

  options.quantity = NRG::Integ::Quantity::energy_moment;
  EXPECT_NEAR(NRG::Integ::calculate(polynomial, options), 2.0 / 3.0, 2e-16);

  options.quantity = NRG::Integ::Quantity::bounded;
  options.requested_lower = -0.5;
  options.requested_upper = 0.25;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), -0.09375);
}

TEST(Integ, bounded_quantities_clip_to_the_polynomial_support) { // NOLINT
  const auto polynomial = linear_identity();
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::bounded;

  options.requested_lower = -2.0;
  options.requested_upper = 0.0;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), -0.5);

  options.requested_lower = 2.0;
  options.requested_upper = 3.0;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.0);

  options.requested_lower = 0.25;
  options.requested_upper = 0.25;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.0);
}

TEST(Integ, analytic_modes_do_not_allocate_qag_workspace) { // NOLINT
  const auto polynomial = linear_identity();
  NRG::Integ::Options options;
  options.workspace_limit = 0;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), 0.0);

  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 0.0;
  EXPECT_DOUBLE_EQ(NRG::Integ::calculate(polynomial, options), -0.5);
}

TEST(Integ, fermi_factor_is_stable_at_extreme_arguments) { // NOLINT
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(1000.0, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(-1000.0, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(0.0, 1e-300), 0.5);
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(-1.0, 0.0), 1.0);
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(0.0, 0.0), 0.5);
  EXPECT_DOUBLE_EQ(NRG::Integ::fermi_factor(1.0, 0.0), 0.0);
}

TEST(Integ, positive_temperature_fermi_integration_splits_piecewise_domain) { // NOLINT
  const NRG::Integ::Polynomial constant{{-2.0, -0.25, 0.5, 3.0}, {{1.0}, {1.0}, {1.0}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 0.37;
  options.epsabs = 1e-13;
  options.epsrel = 1e-12;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  const auto expected = options.temperature
                        * (std::log1p(std::exp(2.0 / options.temperature))
                           - std::log1p(std::exp(-3.0 / options.temperature)));
  EXPECT_NEAR(NRG::Integ::calculate(constant, options), expected, 2e-13);
}

TEST(Integ, fermi_integration_resolves_a_low_temperature_endpoint_layer_without_a_knot) { // NOLINT
  const NRG::Integ::Polynomial constant{{0.0, 1.0}, {{1.0}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1e-6;
  options.epsabs = 1e-20;
  options.epsrel = 1e-10;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  EXPECT_NEAR(NRG::Integ::calculate(constant, options), options.temperature * std::log(2.0), 2e-20);
}

TEST(Integ, fermi_weighting_preserves_a_representable_product_when_the_factor_underflows) { // NOLINT
  const NRG::Integ::Polynomial constant{{750.0, 751.0}, {{1e100}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1.0;
  options.epsabs = 1e-240;
  options.epsrel = 1e-10;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;
  const auto expected = static_cast<double>(
    std::exp(std::log(1e100L) - 750.0L + std::log1p(-std::exp(-1.0L))));

  EXPECT_NEAR(NRG::Integ::calculate(constant, options), expected, expected * 2e-13);
}

TEST(Integ, fermi_integration_resolves_a_representable_far_positive_tail) { // NOLINT
  const auto maximum = std::numeric_limits<double>::max();
  const NRG::Integ::Polynomial constant{{1024.0, 1e6}, {{maximum}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1.0;
  options.epsabs = 1e-150;
  options.epsrel = 0.0;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;
  const auto expected = static_cast<double>(std::exp(std::log(static_cast<long double>(maximum)) - 1024.0L));

  EXPECT_NEAR(NRG::Integ::calculate(constant, options), expected, expected * 3e-13);
}

TEST(Integ, fermi_integration_accepts_local_roundoff_when_the_global_tolerance_is_met) { // NOLINT
  const NRG::Integ::Polynomial constant{{0.0, 1.0}, {{1.0}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1e-6;
  options.epsabs = 1e-16;
  options.epsrel = 0.0;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  EXPECT_NEAR(NRG::Integ::calculate(constant, options), options.temperature * std::log(2.0), 1e-16);
}

TEST(Integ, fermi_integration_resolves_a_subnormal_temperature_layer) { // NOLINT
  const auto maximum = std::numeric_limits<double>::max();
  const auto temperature = std::numeric_limits<double>::denorm_min();
  const NRG::Integ::Polynomial constant{{0.0, 1.0}, {{maximum}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = temperature;
  options.epsabs = 1e-16;
  options.epsrel = 0.0;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;
  const auto expected = static_cast<double>(static_cast<long double>(maximum)
                                            * static_cast<long double>(temperature) * std::log(2.0L));

  EXPECT_NEAR(NRG::Integ::calculate(constant, options), expected, 3e-17);
}

TEST(Integ, fermi_interval_accumulation_does_not_overflow_before_cancellation) { // NOLINT
  const NRG::Integ::Polynomial cancelling{{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                                           {{1.3e308}, {1.3e308}, {1.3e308},
                                            {-1.3e308}, {-1.3e308}, {-1.3e308}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1e300;
  options.epsabs = 1e300;
  options.epsrel = 0.0;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  EXPECT_TRUE(std::isfinite(NRG::Integ::calculate(cancelling, options)));
}

TEST(Integ, fermi_callback_scaling_preserves_a_finite_integral_of_huge_values) { // NOLINT
  const NRG::Integ::Polynomial cancelling{{0.0, 4.0}, {{1e308, -1.7e308}}};
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1e300;
  options.epsabs = 1e300;
  options.epsrel = 0.0;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  EXPECT_NEAR(NRG::Integ::calculate(cancelling, options), 3e307, 2e293);
}

TEST(Integ, global_fermi_tolerance_detects_cancellation_between_intervals) { // NOLINT
  const auto polynomial = linear_identity();
  NRG::Integ::Options options;
  options.quantity = NRG::Integ::Quantity::fermi;
  options.temperature = 1e15;
  options.epsabs = 0.0;
  options.epsrel = 1e-8;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;

  EXPECT_THROW(NRG::Integ::calculate(polynomial, options), std::runtime_error);
}

TEST(Integ, energy_moment_does_not_materialize_overflowing_weighted_coefficients) { // NOLINT
  const auto left = 1e10;
  const auto right = std::nextafter(left, std::numeric_limits<double>::infinity());
  const NRG::Integ::Polynomial constant{{left, right}, {{1e300}}};
  const auto width = static_cast<long double>(right - left);
  const auto expected = static_cast<double>(1e300L * width * (static_cast<long double>(left) + width / 2.0L));

  EXPECT_DOUBLE_EQ(NRG::Integ::energy_moment(constant), expected);
}

TEST(Integ, output_precision_round_trips_double) { // NOLINT
  EXPECT_EQ(NRG::Integ::OUTPUT_PRECISION, std::numeric_limits<double>::max_digits10);
}
