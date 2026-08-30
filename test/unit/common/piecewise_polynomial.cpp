#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmake_configure.hpp>
#include <common/piecewise_polynomial.hpp>

namespace {

using NRG::Tools::CauchyEndpointPolicy;
using NRG::Tools::PiecewisePolynomial;

void expect_complex_near(const std::complex<double> actual, const std::complex<double> expected, const double tolerance) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

auto constant_transform(const std::complex<double> z) { return std::log(z + 1.0) - std::log(z - 1.0); }

auto split_tab_separated(const std::string &line) {
  std::vector<std::string> fields;
  std::istringstream input(line);
  std::string field;
  while (std::getline(input, field, '\t')) fields.push_back(field);
  return fields;
}

auto reference_file(const std::string &name) {
  return std::string(TEST_REFERENCE_DIR) + "/" + name;
}

auto load_reference_polynomials() {
  std::ifstream input(reference_file("polynomials.tsv"));
  if (!input) throw std::runtime_error("Unable to open generated polynomial references.");
  std::map<std::string, std::pair<std::vector<double>, std::vector<std::vector<double>>>> data;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') continue;
    const auto fields = split_tab_separated(line);
    if (fields.size() != 9) throw std::runtime_error("Malformed generated polynomial reference row.");
    const auto degree = static_cast<std::size_t>(std::stoul(fields[4]));
    auto &[knots, coefficients] = data[fields[0]];
    const auto left = std::stod(fields[2]);
    const auto right = std::stod(fields[3]);
    if (knots.empty())
      knots.push_back(left);
    else if (knots.back() != left)
      throw std::runtime_error("Generated polynomial intervals are not contiguous.");
    knots.push_back(right);
    std::vector<double> interval;
    for (std::size_t power = 0; power <= degree; ++power) interval.push_back(std::stod(fields[5 + power]));
    coefficients.push_back(std::move(interval));
  }

  std::map<std::string, PiecewisePolynomial<double>> result;
  for (auto &[name, polynomial] : data)
    result.emplace(name, PiecewisePolynomial<double>{std::move(polynomial.first), std::move(polynomial.second)});
  return result;
}

} // namespace

TEST(PiecewisePolynomial, evaluates_and_integrates_normalized_local_coefficients) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{-2.0, 1.0, 3.0}, {{1.0, 3.0}, {4.0, -2.0, 1.0}}};

  EXPECT_DOUBLE_EQ(polynomial.evaluate(-2.0), 1.0);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(-0.5), 2.5);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(1.0), 4.0);
  EXPECT_DOUBLE_EQ(polynomial.evaluate(3.0), 3.0);
  EXPECT_DOUBLE_EQ(polynomial.integral(), 85.0 / 6.0);
  EXPECT_NEAR(polynomial.integral(-0.5, 2.0), 203.0 / 24.0, 2e-15);
  EXPECT_DOUBLE_EQ(polynomial.integral(-2.0, 3.0), polynomial.integral());
  EXPECT_THROW(polynomial.evaluate(-2.1), std::domain_error);
}

TEST(PiecewisePolynomial, bounded_integral_supports_complex_coefficients_and_stable_primitive_differences) { // NOLINT
  using Complex = std::complex<double>;
  const PiecewisePolynomial<Complex> complex_polynomial{{10.0, 14.0},
                                                         {{{1.0, 2.0}, {3.0, -4.0}, {-2.0, 1.0}}}};
  expect_complex_near(complex_polynomial.integral(11.0, 13.0), {47.0 / 12.0, 13.0 / 24.0}, 2e-15);

  const PiecewisePolynomial<double> cubic_monomial{{0.0, 1.0}, {{0.0, 0.0, 0.0, 1.0}}};
  const auto lower = std::nextafter(1.0, 0.0);
  const auto extended_lower = static_cast<long double>(lower);
  const auto expected = static_cast<double>((1.0L - std::pow(extended_lower, 4)) / 4.0L);
  EXPECT_NEAR(cubic_monomial.integral(lower, 1.0), expected,
              2.0 * std::numeric_limits<double>::epsilon() * expected);

  const auto maximum = std::numeric_limits<double>::max();
  constexpr auto narrow_width = 1e-308;
  const PiecewisePolynomial<double> narrow{{0.0, narrow_width}, {{maximum, maximum}}};
  const auto narrow_expected = maximum * narrow_width + maximum * (narrow_width / 2.0);
  EXPECT_NEAR(narrow.integral(), narrow_expected, 2.0 * std::numeric_limits<double>::epsilon() * narrow_expected);
  EXPECT_NEAR(NRG::Tools::absolute_integral(narrow), narrow_expected,
              2.0 * std::numeric_limits<double>::epsilon() * narrow_expected);

  const auto minimum = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> wide{{0.0, 1e308}, {{maximum}}};
  EXPECT_DOUBLE_EQ(wide.integral(minimum, 2.0 * minimum), maximum * minimum);
  EXPECT_DOUBLE_EQ(NRG::Tools::absolute_integral(wide, minimum, 2.0 * minimum), maximum * minimum);

  const PiecewisePolynomial<double> scaled_linear{{0.0, 1.0}, {{0.0, 1e308}}};
  EXPECT_NEAR(scaled_linear.integral(0.0, 2e-308), 2e-308, 2.0 * minimum);

  const PiecewisePolynomial<double> huge_interval_identity{{0.0, 1e308}, {{0.0, 1e308}}};
  EXPECT_NEAR(huge_interval_identity.integral(0.0, 1e-16), 5e-33, 2e-48);

  const PiecewisePolynomial<double> cancelling{{0.0, 2.0}, {{1e308, -1.7e308}}};
  const auto cancelling_expected = static_cast<double>(
    2.0L * (static_cast<long double>(1e308) + static_cast<long double>(-1.7e308) / 2.0L));
  EXPECT_DOUBLE_EQ(cancelling.integral(), cancelling_expected);

  const PiecewisePolynomial<double> finite_triangle{{0.0, 2.0}, {{maximum, -maximum}}};
  EXPECT_DOUBLE_EQ(NRG::Tools::absolute_integral(finite_triangle), maximum);
}

TEST(PiecewisePolynomial, bounded_integral_validates_bounds) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{0.0, 2.0}, {{1.0, -1.0}}};
  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_DOUBLE_EQ(polynomial.integral(0.5, 0.5), 0.0);
  EXPECT_THROW(polynomial.integral(1.5, 0.5), std::invalid_argument);
  EXPECT_THROW(polynomial.integral(-0.1, 1.0), std::domain_error);
  EXPECT_THROW(polynomial.integral(1.0, 2.1), std::domain_error);
  EXPECT_THROW(polynomial.integral(nan, 1.0), std::invalid_argument);
  EXPECT_THROW(polynomial.integral(0.0, infinity), std::invalid_argument);
}

TEST(PiecewisePolynomial, absolute_integral_handles_linear_crossings_and_subranges) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{2.0, 6.0}, {{-0.25, 1.0}}};

  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial), 5.0 / 4.0, 2e-15);
  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial, 2.5, 4.0), 5.0 / 32.0, 2e-15);
  EXPECT_DOUBLE_EQ(NRG::Tools::absolute_integral(polynomial, 3.0, 3.0), 0.0);
}

TEST(PiecewisePolynomial, absolute_integral_handles_degenerate_quadratic_tangent) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{-1.0, 2.0}, {{0.25, -1.0, 1.0, 0.0}}};

  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial), 0.25, 2e-15);
}

TEST(PiecewisePolynomial, absolute_integral_partitions_cubic_with_three_roots) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{0.0, 1.0}, {{-3.0 / 32.0, 11.0 / 16.0, -1.5, 1.0}}};

  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial), 5.0 / 256.0, 2e-15);
}

TEST(PiecewisePolynomial, absolute_integral_handles_endpoint_and_shared_knot_roots) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{0.0, 2.0, 5.0}, {{0.0, 1.0, -1.0}, {0.0, 1.0, -1.0}}};

  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial), 5.0 / 6.0, 2e-15);
  EXPECT_NEAR(NRG::Tools::absolute_integral(polynomial, 2.0, 5.0), 0.5, 2e-15);
}

TEST(PiecewisePolynomial, absolute_integral_handles_zero_and_constant_intervals) { // NOLINT
  const PiecewisePolynomial<double> polynomial{{0.0, 1.0, 3.0}, {{0.0, 0.0, 0.0, 0.0}, {-2.0}}};

  EXPECT_DOUBLE_EQ(NRG::Tools::absolute_integral(polynomial, 0.0, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(NRG::Tools::absolute_integral(polynomial), 4.0);
}

TEST(PiecewisePolynomial, absolute_integral_validates_bounds_and_degree) { // NOLINT
  const PiecewisePolynomial<double> linear{{0.0, 1.0}, {{1.0, -2.0}}};
  const PiecewisePolynomial<double> quartic{{0.0, 1.0}, {{0.0, 0.0, 0.0, 1.0, 1.0}}};
  const auto infinity = std::numeric_limits<double>::infinity();

  EXPECT_THROW(NRG::Tools::absolute_integral(linear, 0.8, 0.2), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::absolute_integral(linear, -0.1, 0.8), std::domain_error);
  EXPECT_THROW(NRG::Tools::absolute_integral(linear, 0.2, infinity), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::absolute_integral(quartic), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::absolute_integral(quartic, 0.2, 0.8), std::invalid_argument);
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
  const std::vector<double> knots{-1.0, -0.25, 0.5, 1.0};
  const PiecewisePolynomial<double> split{knots,
                                           {{knots[0], knots[1] - knots[0]},
                                            {knots[1], knots[2] - knots[1]},
                                            {knots[2], knots[3] - knots[2]}}};
  for (const std::complex<double> z : {std::complex<double>{-0.25, 1e-12}, std::complex<double>{0.3, 0.4}})
    expect_complex_near(NRG::Tools::cauchy_transform(split, z), NRG::Tools::cauchy_transform(one_piece, z), 3e-13);
  for (const double x : {-0.25, 0.0, 0.5})
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

TEST(PiecewisePolynomial, energy_weighting_scales_before_narrowing_to_double) { // NOLINT
  const PiecewisePolynomial<double> large{{1e200, 2e200}, {{1e-200}}};
  const auto large_weighted = large.multiply_by_monomial(2);
  EXPECT_NEAR(large_weighted.coefficients()[0][0], 1e200, 2e184);
  EXPECT_NEAR(large_weighted.coefficients()[0][1], 2e200, 4e184);
  EXPECT_NEAR(large_weighted.coefficients()[0][2], 1e200, 2e184);
  EXPECT_NEAR(large_weighted.evaluate(1.5e200), 2.25e200, 5e184);

  const PiecewisePolynomial<double> small{{1e-200, 2e-200}, {{1e200}}};
  const auto small_weighted = small.multiply_by_monomial(2);
  EXPECT_NEAR(small_weighted.coefficients()[0][0], 1e-200, 2e-216);
  EXPECT_NEAR(small_weighted.coefficients()[0][1], 2e-200, 4e-216);
  EXPECT_NEAR(small_weighted.coefficients()[0][2], 1e-200, 2e-216);
  EXPECT_NEAR(small_weighted.evaluate(1.5e-200), 2.25e-200, 5e-216);

  const PiecewisePolynomial<double> zero{{1e300, 2e300}, {{0.0}}};
  const auto zero_weighted = zero.multiply_by_monomial(20);
  EXPECT_DOUBLE_EQ(zero_weighted.evaluate(1.5e300), 0.0);
  EXPECT_THROW(zero.multiply_by_monomial(std::numeric_limits<std::size_t>::max()), std::length_error);

  const auto minimum = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> subnormal{{0.4, 0.8}, {{minimum, minimum}}};
  const auto subnormal_weighted = subnormal.multiply_by_monomial(1);
  EXPECT_DOUBLE_EQ(subnormal_weighted.coefficients()[0][1], minimum);

  const PiecewisePolynomial<double> cancelling{{1.0, 2.0}, {{1e300, 0.5, -1e300}}};
  const auto cancelling_weighted = cancelling.multiply_by_monomial(2);
  EXPECT_DOUBLE_EQ(cancelling_weighted.coefficients()[0][2], 1.0);
}

TEST(PiecewisePolynomial, rejects_domains_whose_width_overflows) { // NOLINT
  const auto maximum = std::numeric_limits<double>::max();
  EXPECT_THROW((PiecewisePolynomial<double>{{-maximum, maximum}, {{1.0}}}), std::invalid_argument);
  EXPECT_THROW((PiecewisePolynomial<double>{{-maximum, 0.0, maximum}, {{1.0}, {1.0}}}), std::invalid_argument);
}

TEST(PiecewisePolynomial, matches_stored_multiprecision_transform_references) { // NOLINT
  const auto polynomials = load_reference_polynomials();
  std::ifstream input(reference_file("cauchy_transforms.tsv"));
  ASSERT_TRUE(input);

  std::string line;
  std::size_t cases = 0;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') continue;
    const auto fields = split_tab_separated(line);
    ASSERT_EQ(fields.size(), 8) << line;
    SCOPED_TRACE(fields[0]);
    const auto polynomial = polynomials.at(fields[1]).multiply_by_monomial(std::stoul(fields[2]));
    const std::complex<double> argument{std::stod(fields[4]), std::stod(fields[5])};
    const std::complex<double> expected{std::stod(fields[6]), std::stod(fields[7])};
    if (fields[0] == "legendre3_threshold_near") {
      EXPECT_GE(std::abs(2.0 / (argument + 1.0)), 0.5);
    }
    if (fields[0] == "legendre3_threshold_far") {
      EXPECT_LT(std::abs(2.0 / (argument + 1.0)), 0.5);
    }
    std::complex<double> actual;
    if (fields[3] == "complex")
      actual = NRG::Tools::cauchy_transform(polynomial, argument);
    else if (fields[3] == "real_pv")
      actual = NRG::Tools::cauchy_principal_value(polynomial, argument.real());
    else
      FAIL() << "Unknown generated argument kind: " << fields[3];
    const auto tolerance = 3e-13 * std::abs(expected) + 32.0 * std::numeric_limits<double>::denorm_min();
    expect_complex_near(actual, expected, tolerance);
    ++cases;
  }
  EXPECT_EQ(cases, 23);
}

TEST(PiecewisePolynomial, matches_stored_multiprecision_moment_references) { // NOLINT
  const auto polynomials = load_reference_polynomials();
  std::ifstream input(reference_file("moments.tsv"));
  ASSERT_TRUE(input);

  std::string line;
  std::size_t cases = 0;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') continue;
    const auto fields = split_tab_separated(line);
    ASSERT_EQ(fields.size(), 3) << line;
    SCOPED_TRACE(fields[0] + " power " + fields[1]);
    const auto actual = polynomials.at(fields[0]).multiply_by_monomial(std::stoul(fields[1])).integral();
    const auto expected = std::stod(fields[2]);
    EXPECT_NEAR(actual, expected, 2e-13 * std::max(1.0, std::abs(expected)));
    ++cases;
  }
  EXPECT_EQ(cases, 16);
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

  const auto maximum = std::numeric_limits<double>::max();
  const PiecewisePolynomial<std::complex<double>> large_complex{{0.0, 1.0}, {{{maximum, maximum}}}};
  const std::complex<double> moderate_z{4.0, 1.0};
  const auto unit_transform = std::log(moderate_z) - std::log(moderate_z - 1.0);
  const std::complex<double> large_expected{maximum * (unit_transform.real() - unit_transform.imag()),
                                             maximum * (unit_transform.real() + unit_transform.imag())};
  const auto large_actual = NRG::Tools::cauchy_transform(large_complex, moderate_z);
  EXPECT_NEAR(large_actual.real() / large_expected.real(), 1.0, 3e-14);
  EXPECT_NEAR(large_actual.imag() / large_expected.imag(), 1.0, 3e-14);

  const PiecewisePolynomial<double> large_moments{{0.0, 1.0}, {{maximum, maximum}}};
  constexpr double real_argument = 10.0;
  const auto large_factor = 11.0L * std::log1p(1.0L / 9.0L) - 1.0L;
  const auto large_moment_result = NRG::Tools::cauchy_principal_value(large_moments, real_argument).real();
  EXPECT_NEAR(large_moment_result / (maximum * static_cast<double>(large_factor)), 1.0, 3e-14);

  const PiecewisePolynomial<double> zero_transform{{0.0, 1.0}, {{-1.5, 3.5, -1.0}}};
  EXPECT_DOUBLE_EQ(NRG::Tools::cauchy_principal_value(zero_transform, 3.0).real(), 0.0);

  const auto wide_right = std::ldexp(1.0, 53);
  const PiecewisePolynomial<double> wide_domain{{-1.0, wide_right}, {{1.0}}};
  EXPECT_DOUBLE_EQ(NRG::Tools::cauchy_principal_value(wide_domain, std::ldexp(1.0, 55)).real(),
                   0.28768207245178096);

  const PiecewisePolynomial<std::complex<double>> separated_scales{
    {0.0, 1.0}, {{{-1.0, 1e300}, {12.0, 0.0}, {-30.0, 0.0}, {20.0, 0.0}}}};
  constexpr double separated_argument = 1e10;
  const auto inverse_separated = 1.0 / separated_argument;
  const auto separated_expected = 1.0 / 140.0 * std::pow(inverse_separated, 4)
                                  + 1.0 / 70.0 * std::pow(inverse_separated, 5);
  const auto separated_result = NRG::Tools::cauchy_principal_value(separated_scales, separated_argument);
  EXPECT_NEAR(separated_result.real() / separated_expected, 1.0, 3e-14);

  const auto cancellation_scale = std::ldexp(1.0, 40);
  const PiecewisePolynomial<std::complex<double>> delayed_component{
    {0.0, 1.0}, {{{9.0, 36.0 * cancellation_scale}, {-36.0, -192.0 * cancellation_scale},
                  {30.0, 180.0 * cancellation_scale}}}};
  const auto delayed_result = NRG::Tools::cauchy_transform(
    delayed_component, std::complex<double>{0.0, cancellation_scale});
  EXPECT_NEAR(delayed_result.real() / 6.5979105986e-49, 1.0, 3e-10);
  EXPECT_NEAR(delayed_result.imag() / 4.5138983072e-37, 1.0, 3e-10);

  const auto large_scale = std::ldexp(1.0, 990);
  const std::vector<double> scaled_coefficients{-large_scale, 12.0 * large_scale,
                                                 -30.0 * large_scale, 20.0 * large_scale};
  constexpr double extreme_argument = 1e100;
  auto extreme_expected = large_scale / 140.0;
  for (int power = 0; power < 4; ++power) extreme_expected /= extreme_argument;
  const PiecewisePolynomial<double> scaled_legendre{{0.0, 1.0}, {scaled_coefficients}};
  const auto scaled_result = NRG::Tools::cauchy_principal_value(scaled_legendre, extreme_argument).real();
  EXPECT_NEAR(scaled_result / extreme_expected, 1.0, 3e-14);

  const PiecewisePolynomial<double> locally_far{{0.0, 1.0, 100.0}, {scaled_coefficients, {0.0}}};
  const std::complex<double> local_argument{99.0, 1.0};
  const auto local_expected = NRG::Tools::cauchy_transform(scaled_legendre, local_argument);
  const auto local_result = NRG::Tools::cauchy_transform(locally_far, local_argument);
  EXPECT_NEAR(local_result.real() / local_expected.real(), 1.0, 3e-14);
  EXPECT_NEAR(local_result.imag() / local_expected.imag(), 1.0, 3e-14);

  std::vector<double> many_knots;
  std::vector<std::vector<double>> many_coefficients;
  for (int knot = 0; knot <= 66; ++knot) many_knots.push_back(static_cast<double>(knot));
  for (int interval = 0; interval < 65; ++interval) many_coefficients.push_back({0.0});
  many_coefficients.push_back(scaled_coefficients);
  const PiecewisePolynomial<double> many_intervals{std::move(many_knots), std::move(many_coefficients)};
  const auto many_result = NRG::Tools::cauchy_principal_value(many_intervals, extreme_argument).real();
  EXPECT_NEAR(many_result / extreme_expected, 1.0, 3e-13);
}

TEST(PiecewisePolynomial, far_field_preserves_subnormal_interval_widths) { // NOLINT
  const auto width = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> narrow{{-2.0, 0.0, width, 2.0}, {{0.0}, {1e308}, {0.0}}};
  const std::complex<double> z{10.0, 1.0};
  const auto expected = (width * 1e308) / z;

  expect_complex_near(NRG::Tools::cauchy_transform(narrow, z), expected, 2e-31);

  const auto maximum = std::numeric_limits<double>::max();
  const PiecewisePolynomial<double> extreme{{0.0, width}, {{maximum}}};
  EXPECT_DOUBLE_EQ(NRG::Tools::cauchy_principal_value(extreme, maximum).real(), width);

  const PiecewisePolynomial<double> power_underflow{{0.0, 1.0}, {{-maximum / 2.0, maximum}}};
  constexpr double huge_argument = 1e308;
  const auto power_underflow_expected = (maximum / huge_argument) / 12.0 / huge_argument;
  const auto power_underflow_result = NRG::Tools::cauchy_principal_value(power_underflow, huge_argument).real();
  EXPECT_NEAR(power_underflow_result / power_underflow_expected, 1.0, 3e-14);
}

TEST(PiecewisePolynomial, local_transforms_preserve_subnormal_interval_widths) { // NOLINT
  const auto width = std::numeric_limits<double>::denorm_min();
  const PiecewisePolynomial<double> spike{{0.0, width, 2.0 * width, 4.0},
                                           {{0.0, 1e308}, {1e308, -1e308}, {0.0}}};
  const auto mass = width * 1e308;
  for (const std::complex<double> z : {std::complex<double>{1.0, 0.1}, std::complex<double>{2.0, 1.0}})
    expect_complex_near(NRG::Tools::cauchy_transform(spike, z), mass / z, 2e-31);
}

TEST(PiecewisePolynomial, wide_near_path_preserves_binary64_endpoint_differences) { // NOLINT
  if constexpr (!NRG::Tools::detail::native_extended_precision) {
    const auto scale = std::ldexp(1.0, 1023);
    const auto complex_result = NRG::Tools::detail::wide_complex_interval_transform(
      std::vector<double>{scale}, {}, 0.0, scale,
      std::complex<double>{3.0 * std::ldexp(1.0, 1022), std::ldexp(1.0, -511)});
    EXPECT_NEAR(complex_result.imag() / -1.9888908616533885e-154, 1.0, 3e-14);
  }

  const auto real_result = NRG::Tools::detail::wide_real_interval_transform(
    std::vector<double>{1.0, -2.5, 1.0}, {}, -1.0, std::ldexp(1.0, 53), std::ldexp(1.0, 54));
  EXPECT_NEAR(real_result.real() / -4.409891434233644e-18, 1.0, 3e-14);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
