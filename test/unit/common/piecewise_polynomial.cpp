#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include <common/piecewise_polynomial.hpp>

namespace {

using NRG::Tools::CauchyEndpointPolicy;
using NRG::Tools::PiecewisePolynomial;

void expect_complex_near(const std::complex<double> actual, const std::complex<double> expected, const double tolerance) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

auto constant_transform(const std::complex<double> z) { return std::log(z + 1.0) - std::log(z - 1.0); }

} // namespace

TEST(PiecewisePolynomial, evaluates_and_integrates_normalized_local_coefficients) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{-2.0, 1.0, 3.0}, {{1.0, 3.0}, {4.0, -2.0, 1.0}}};

  EXPECT_DOUBLE_EQ(polynomial.evaluate(-2.0), 1.0);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(-0.5), 2.5);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(1.0), 4.0);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(3.0), 3.0);
  EXPECT_NEAR(polynomial.integral(), 3.0 * 2.5 + 2.0 * (4.0 - 1.0 + 1.0 / 3.0), 1e-15);
  EXPECT_THROW(polynomial.evaluate(-2.1), std::domain_error);
}

TEST(PiecewisePolynomial, exact_cauchy_transform_matches_polynomial_identities) { // NOLINT
  const PiecewisePolynomial<double> constant{{-1.0, 1.0}, {{1.0}}};
  const PiecewisePolynomial<double> linear{{-1.0, 1.0}, {{-1.0, 2.0}}};
  const PiecewisePolynomial<double> quadratic{{-1.0, 1.0}, {{1.0, -4.0, 4.0}}};
  const PiecewisePolynomial<double> cubic{{-1.0, 1.0}, {{-1.0, 6.0, -12.0, 8.0}}};

  for (const std::complex<double> z : {std::complex<double>{0.2, 0.3}, std::complex<double>{-0.4, -0.2},
                                       std::complex<double>{1.5, 0.7}}) {
    const auto c0 = constant_transform(z);
    const auto c1 = z * c0 - 2.0;
    const auto c2 = z * c1;
    const auto c3 = z * c2 - 2.0 / 3.0;
    expect_complex_near(NRG::Tools::cauchy_transform(constant, z), c0, 2e-14);
    expect_complex_near(NRG::Tools::cauchy_transform(linear, z), c1, 2e-14);
    expect_complex_near(NRG::Tools::cauchy_transform(quadratic, z), c2, 3e-14);
    expect_complex_near(NRG::Tools::cauchy_transform(cubic, z), c3, 5e-14);
  }
}

TEST(PiecewisePolynomial, local_far_expansion_does_not_stop_at_a_vanishing_moment) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{0.0, 1.0, 2.0}, {{1.0, -1.2}, {-0.2}}};
  const std::complex<double> z{2.5, 0.01};
  const auto first_logarithm = std::log(z) - std::log(z - 1.0);
  const auto second_logarithm = std::log(z - 1.0) - std::log(z - 2.0);
  const auto expected = (1.0 - 1.2 * z) * first_logarithm + 1.2 - 0.2 * second_logarithm;

  expect_complex_near(NRG::Tools::cauchy_transform(polynomial, z), expected, 2e-14);
}

TEST(PiecewisePolynomial, principal_value_and_endpoint_subtraction_are_explicit) { // NOLINT
  const PiecewisePolynomial<double> constant{{-1.0, 0.0, 1.0}, {{1.0}, {1.0}}};
  const PiecewisePolynomial<double> linear{{-1.0, 0.0, 1.0}, {{-1.0, 1.0}, {0.0, 1.0}}};

  for (const double x : {-0.7, 0.0, 0.4}) {
    const auto expected_constant = std::log(std::abs(x + 1.0)) - std::log(std::abs(x - 1.0));
    EXPECT_NEAR(NRG::Tools::cauchy_principal_value(constant, x).real(), expected_constant, 2e-14);
    EXPECT_NEAR(NRG::Tools::cauchy_principal_value(linear, x).real(), x * expected_constant - 2.0, 2e-14);
  }

  EXPECT_THROW(NRG::Tools::cauchy_principal_value(constant, 1.0), std::domain_error);
  EXPECT_DOUBLE_EQ(NRG::Tools::cauchy_principal_value(constant, 1.0, CauchyEndpointPolicy::subtracted).real(), 0.0);
  EXPECT_NEAR(NRG::Tools::cauchy_principal_value(linear, 1.0, CauchyEndpointPolicy::subtracted).real(), -2.0, 1e-15);
}

TEST(PiecewisePolynomial, splitting_intervals_does_not_change_the_transform) { // NOLINT
  const PiecewisePolynomial<double> one_piece{{-1.0, 1.0}, {{-1.0, 2.0}}};
  const PiecewisePolynomial<double> split{{-1.0, -0.2, 0.4, 1.0},
                                           {{-1.0, 0.8}, {-0.2, 0.6}, {0.4, 0.6}}};
  for (const std::complex<double> z : {std::complex<double>{-0.2, 1e-12}, std::complex<double>{0.3, 0.4}})
    expect_complex_near(NRG::Tools::cauchy_transform(split, z), NRG::Tools::cauchy_transform(one_piece, z), 3e-13);
  for (const double x : {-0.2, 0.0, 0.4})
    expect_complex_near(NRG::Tools::cauchy_principal_value(split, x), NRG::Tools::cauchy_principal_value(one_piece, x), 3e-13);

  const std::complex<double> far_z{1e16, 2e15};
  const auto inverse = 1.0 / far_z;
  const auto far_expected = 2.0 / 3.0 * inverse * inverse;
  expect_complex_near(NRG::Tools::cauchy_transform(one_piece, far_z), far_expected, 2e-47);
  expect_complex_near(NRG::Tools::cauchy_transform(split, far_z), far_expected, 2e-47);
  EXPECT_NEAR(NRG::Tools::cauchy_principal_value(split, 1e16).real(), 2.0 / (3.0 * 1e32), 2e-47);
}

TEST(PiecewisePolynomial, principal_value_rejects_discontinuous_internal_knots) { // NOLINT
  const PiecewisePolynomial<double> discontinuous{{-1.0, 0.0, 1.0}, {{1.0}, {2.0}}};

  EXPECT_THROW(NRG::Tools::cauchy_principal_value(discontinuous, 0.0), std::domain_error);
  EXPECT_NO_THROW(NRG::Tools::cauchy_principal_value(discontinuous, 0.2));
}

TEST(PiecewisePolynomial, energy_weighting_preserves_interpolate_then_multiply_semantics) { // NOLINT
  const PiecewisePolynomial<double> constant{{-1.0, 1.0}, {{0.5}}};
  const auto weighted = constant.multiply_by_monomial(2);
  const std::complex<double> z{0.2, 0.3};
  const auto c0 = constant_transform(z);
  const auto expected = 0.5 * (z * (z * c0 - 2.0));

  EXPECT_NEAR(weighted.integral(), 1.0 / 3.0, 2e-15);
  expect_complex_near(NRG::Tools::cauchy_transform(weighted, z), expected, 3e-14);
}

TEST(PiecewisePolynomial, supports_complex_coefficients_and_lower_half_plane_conjugation) { // NOLINT
  const std::complex<double> scale{1.0, 0.25};
  const PiecewisePolynomial<std::complex<double>> polynomial{{-1.0, 1.0}, {{scale}}};
  const std::complex<double> upper{0.3, 1e-8};
  const std::complex<double> lower{0.3, -1e-8};

  expect_complex_near(NRG::Tools::cauchy_transform(polynomial, upper), scale * constant_transform(upper), 2e-13);
  expect_complex_near(NRG::Tools::cauchy_transform(polynomial, lower), scale * constant_transform(lower), 2e-13);
}

TEST(PiecewisePolynomial, far_field_uses_stable_moment_expansion) { // NOLINT
  const PiecewisePolynomial<double> constant{{-1.0, 1.0}, {{1.0}}};
  const std::complex<double> z{1e12, 2e11};
  const auto inverse = 1.0 / z;
  const auto expected = 2.0 * inverse + 2.0 / 3.0 * inverse * inverse * inverse;
  expect_complex_near(NRG::Tools::cauchy_transform(constant, z), expected, 1e-38);
}

TEST(PiecewisePolynomial, far_field_preserves_subnormal_interval_widths) { // NOLINT
  const auto width = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> narrow{{-2.0, 0.0, width, 2.0}, {{0.0}, {1e308}, {0.0}}};
  const std::complex<double> z{10.0, 1.0};
  const auto expected = (width * 1e308) / z;

  expect_complex_near(NRG::Tools::cauchy_transform(narrow, z), expected, 2e-31);
}

TEST(PiecewisePolynomial, local_transforms_preserve_subnormal_interval_widths) { // NOLINT
  const auto width = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> spike{{0.0, width, 2.0 * width, 4.0},
                                           {{0.0, 1e308}, {1e308, -1e308}, {0.0}}};
  const auto mass = width * 1e308;
  for (const std::complex<double> z : {std::complex<double>{1.0, 0.1}, std::complex<double>{2.0, 1.0}})
    expect_complex_near(NRG::Tools::cauchy_transform(spike, z), mass / z, 2e-31);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
