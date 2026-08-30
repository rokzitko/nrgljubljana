#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <ios>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include <hilb/hilb.hpp>

namespace {

constexpr double B = 1.0;
int gsl_handler_calls = 0;

void count_gsl_handler_calls([[maybe_unused]] const char *reason, [[maybe_unused]] const char *file,
                             [[maybe_unused]] int line, [[maybe_unused]] int error) {
  ++gsl_handler_calls;
}

auto flat_band_h0(const std::complex<double> z) { return (std::log(z + B) - std::log(z - B)) / (2.0 * B); }

auto flat_band_h0_far(const std::complex<double> z) {
  const auto inverse = B / z;
  return (inverse + inverse * inverse * inverse / 3.0
          + inverse * inverse * inverse * inverse * inverse / 5.0) / B;
}

void expect_complex_near(const std::complex<double> actual, const std::complex<double> expected, const double tolerance) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

void run_hilb(std::vector<std::string> arguments) {
  std::vector<char *> argv;
  argv.reserve(arguments.size());
  for (auto &argument : arguments) argv.push_back(argument.data());
  optind = 1;
  NRG::Hilb::Hilb(static_cast<int>(argv.size()), argv.data());
}

struct captured_output {
  std::string out;
  std::string err;
};

auto run_hilb_captured(std::vector<std::string> arguments) {
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  try {
    run_hilb(std::move(arguments));
  } catch (...) {
    testing::internal::GetCapturedStderr();
    testing::internal::GetCapturedStdout();
    throw;
  }
  const auto err = testing::internal::GetCapturedStderr();
  const auto out = testing::internal::GetCapturedStdout();
  return captured_output{out, err};
}

void write_file(const std::string &filename, const std::string &contents) {
  std::ofstream file(filename);
  file << contents;
}

auto read_file(const std::string &filename) {
  std::ifstream file(filename);
  std::ostringstream contents;
  contents << file.rdbuf();
  return contents.str();
}

auto parse_complex_output(const std::string &output) {
  std::istringstream input(output);
  char open = 0;
  char comma = 0;
  char close = 0;
  double real = 0.0;
  double imaginary = 0.0;
  if (!(input >> open >> real >> comma >> imaginary >> close) || open != '(' || comma != ',' || close != ')')
    throw std::runtime_error("Expected one complex numeric value.");
  input >> std::ws;
  if (!input.eof()) throw std::runtime_error("Unexpected trailing complex-output data.");
  return std::complex<double>{real, imaginary};
}

auto count_occurrences(const std::string &text, const std::string &needle) {
  size_t count = 0;
  size_t position = 0;
  while ((position = text.find(needle, position)) != std::string::npos) {
    ++count;
    position += needle.size();
  }
  return count;
}

} // namespace

TEST(Hilb, empty_dos_throws) { // NOLINT
  const auto filename = "hilb_empty_dos.dat";
  std::ofstream file(filename);

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_empty_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_THROW(NRG::Hilb::Hilb(5, argv), std::runtime_error);

  std::remove(filename);
}

TEST(Hilb, unsorted_dos_is_sorted_before_use) { // NOLINT
  const auto filename = "hilb_unsorted_dos.dat";
  {
    std::ofstream file(filename);
    file << "1 1\n-1 1\n0 1\n";
  }

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_unsorted_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_NO_THROW(NRG::Hilb::Hilb(5, argv));

  std::remove(filename);
}

TEST(Hilb, interpolator_move_assignment_releases_old_state) { // NOLINT
  auto a = NRG::Hilb::interpolator(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{1.0, 2.0, 3.0});
  auto b = NRG::Hilb::interpolator(std::vector<double>{-2.0, 0.0, 2.0}, std::vector<double>{4.0, 5.0, 6.0});
  a = std::move(b);
  EXPECT_DOUBLE_EQ(a(0.0), 5.0);
}

TEST(Hilb, interpolator_rejects_invalid_grids) { // NOLINT
  EXPECT_THROW(NRG::Hilb::interpolator(std::vector<double>{-1.0, 1.0}, std::vector<double>{1.0, 1.0}), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::interpolator(std::vector<double>{-1.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::interpolator(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{1.0, 1.0}), std::invalid_argument);
}

TEST(Hilb, interpolator_supports_all_interpolation_methods) { // NOLINT
  using NRG::Tools::InterpolationMethod;

  NRG::Hilb::interpolator linear({0.0, 2.0}, {1.0, 5.0}, -7.0, InterpolationMethod::linear);
  EXPECT_DOUBLE_EQ(linear(0.5), 2.0);
  EXPECT_DOUBLE_EQ(linear(-0.1), -7.0);

  const std::vector<double> cubic_x{0.0, 1.0, 2.0};
  const std::vector<double> cubic_y{0.0, 1.0, 0.0};
  NRG::Hilb::interpolator legacy_cubic(cubic_x, cubic_y, -7.0);
  NRG::Hilb::interpolator explicit_cubic(cubic_x, cubic_y, -7.0, InterpolationMethod::cspline);
  NRG::Hilb::interpolator piecewise_linear(cubic_x, cubic_y, 0.0, InterpolationMethod::linear);
  NRG::Hilb::interpolator steffen(cubic_x, cubic_y, 0.0, InterpolationMethod::steffen);
  EXPECT_DOUBLE_EQ(legacy_cubic(0.5), explicit_cubic(0.5));
  EXPECT_NEAR(explicit_cubic(0.5), 0.6875, 1e-15);
  EXPECT_DOUBLE_EQ(piecewise_linear(0.5), 0.5);
  EXPECT_NEAR(steffen(0.5), 0.625, 1e-15);
  EXPECT_NE(explicit_cubic(0.5), piecewise_linear(0.5));
  EXPECT_DOUBLE_EQ(legacy_cubic(2.1), -7.0);

  NRG::Hilb::interpolator akima({0.0, 1.0, 2.0, 3.0, 4.0}, {0.0, 1.0, 0.0, 1.0, 0.0}, 0.0,
                                 InterpolationMethod::akima);
  EXPECT_NEAR(akima(1.25), 0.84375, 1e-14);

  NRG::Hilb::interpolator legacy_brace_argument(cubic_x, cubic_y, {});
  EXPECT_DOUBLE_EQ(legacy_brace_argument(0.5), legacy_cubic(0.5));
}

TEST(Hilb, interpolator_enforces_each_method_minimum_size) { // NOLINT
  using NRG::Tools::InterpolationMethod;

  EXPECT_THROW(NRG::Hilb::interpolator({0.0}, {1.0}, 0.0, InterpolationMethod::linear), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::interpolator({0.0, 1.0}, {0.0, 1.0}, 0.0, InterpolationMethod::cspline), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::interpolator({0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 0.0, 1.0}, 0.0, InterpolationMethod::akima),
               std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::interpolator({0.0, 1.0}, {0.0, 1.0}, 0.0, InterpolationMethod::steffen), std::invalid_argument);

  EXPECT_NO_THROW(NRG::Hilb::interpolator({0.0, 1.0}, {0.0, 1.0}, 0.0, InterpolationMethod::linear));
  EXPECT_NO_THROW(NRG::Hilb::interpolator({0.0, 1.0, 2.0}, {0.0, 1.0, 0.0}, 0.0, InterpolationMethod::cspline));
  EXPECT_NO_THROW(NRG::Hilb::interpolator({0.0, 1.0, 2.0, 3.0, 4.0}, {0.0, 1.0, 0.0, 1.0, 0.0}, 0.0,
                                           InterpolationMethod::akima));
  EXPECT_NO_THROW(NRG::Hilb::interpolator({0.0, 1.0, 2.0}, {0.0, 1.0, 0.0}, 0.0, InterpolationMethod::steffen));
}

TEST(Hilb, energy_power_direct_quadrature) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{0.2, 0.1};
  const auto expected = z * flat_band_h0(z) - 1.0;

  const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, z, 1e-3, 1);

  expect_complex_near(result, expected, 1e-10);
}

TEST(Hilb, energy_power_singularity_subtraction) { // NOLINT
  auto rhor = [](const double) { return 0.5; };
  auto rhoi = [](const double) { return 0.125; };
  const std::complex<double> z{0.2, 1e-6};
  const auto real_density_result = z * z * flat_band_h0(z) - z;
  const auto expected = std::complex<double>{1.0, 0.25} * real_density_result;

  const auto result = NRG::Hilb::hilbert_transform(rhor, rhoi, B, z, 1e-3, 2);

  expect_complex_near(result, expected, 1e-9);
}

TEST(Hilb, singularity_subtraction_integrates_both_band_edges) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  for (const double x : {-B, B}) {
    const std::complex<double> positive_z{x, 1e-6};
    const auto positive_expected = positive_z * flat_band_h0(positive_z) - 1.0;
    const auto positive_result = NRG::Hilb::hilbert_transform(rho, zero, B, positive_z, 1e-3, 1);
    expect_complex_near(positive_result, positive_expected, 1e-9);

    const std::complex<double> negative_z{x, -1e-6};
    const auto negative_result = NRG::Hilb::hilbert_transform(rho, zero, B, negative_z, 1e-3, 1);
    expect_complex_near(negative_result, std::conj(positive_expected), 1e-9);
  }
}

TEST(Hilb, subtraction_does_not_evaluate_density_outside_support) { // NOLINT
  auto rho = [](const double energy) {
    if (std::abs(energy) > B) throw std::runtime_error("outside support");
    return 0.5;
  };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{2.0, 1e-6};
  EXPECT_NO_THROW({
    const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, z);
    expect_complex_near(result, flat_band_h0(z), 1e-9);
  });
}

TEST(Hilb, far_outside_callable_transform_does_not_collapse_or_overflow) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  const auto minimum_imaginary = NRG::Hilb::minimum_safe_imaginary_part();

  for (const double x : {9007199254740992.0, -9007199254740992.0, 1e200, -1e200}) {
    const std::complex<double> z{x, minimum_imaginary};
    const auto expected = flat_band_h0_far(z);
    const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, z);

    EXPECT_NE(result.real(), 0.0);
    EXPECT_NEAR(result.real() / expected.real(), 1.0, 2e-14);
    if (expected.imag() == 0.0)
      EXPECT_EQ(result.imag(), 0.0);
    else
      EXPECT_NEAR(result.imag() / expected.imag(), 1.0, 2e-14);
  }
}

TEST(Hilb, direct_callable_transform_avoids_finite_complex_division_overflow) { // NOLINT
  const auto maximum = std::numeric_limits<double>::max();
  const auto minimum = std::numeric_limits<double>::denorm_min();
  const auto quotient = NRG::Hilb::robust_complex_divide({maximum, minimum}, 0.0, 1.0);
  EXPECT_DOUBLE_EQ(quotient.real(), minimum);
  EXPECT_DOUBLE_EQ(quotient.imag(), -maximum);

  auto density = [maximum](const double) { return maximum; };
  constexpr double narrow_bandwidth = 1e-100;
  const auto result = NRG::Hilb::hilbert_transform(density, density, narrow_bandwidth,
                                                    std::complex<double>{1.0, 1.0});

  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));
  EXPECT_NEAR(result.real() / (2.0 * narrow_bandwidth * maximum), 1.0, 2e-14);
  EXPECT_NEAR(result.imag() / result.real(), 0.0, 2e-14);

  auto large_density = [](const double) { return 1e200; };
  auto small_density = [](const double) { return 1e-200; };
  NRG::Hilb::integrator integration(NRG::Hilb::integrator::configured, 64, NRG::Tools::QagRule::gauss61);
  constexpr double mixed_bandwidth = 0.3;
  const auto mixed = NRG::Hilb::hilbert_transform(integration, large_density, small_density, mixed_bandwidth,
                                                  std::complex<double>{0.0, 1.0}, 0.1, 0, 1e-220, 1e-12);
  const auto factor = 2.0 * std::atan(mixed_bandwidth);
  EXPECT_NEAR(mixed.real() / (1e-200 * factor), 1.0, 3e-14);
  EXPECT_NEAR(mixed.imag() / (-1e200 * factor), 1.0, 3e-14);

  constexpr auto coupled_bandwidth = std::ldexp(1.0, -346);
  auto coupled_real = [](const double) { return std::ldexp(1.0, 991); };
  auto coupled_imaginary = [](const double) { return std::ldexp(1.0, 660); };
  NRG::Hilb::integrator coupled_integration(NRG::Hilb::integrator::configured, 64,
                                             NRG::Tools::QagRule::gauss61);
  const auto coupled = NRG::Hilb::hilbert_transform(
    coupled_integration, coupled_real, coupled_imaginary, coupled_bandwidth,
    std::complex<double>{std::ldexp(1.0, -180), std::ldexp(1.0, -511)}, 0.0, 0, 0.0, 1e-12);
  EXPECT_NEAR(coupled.imag() / -7.794675399098149e48, 1.0, 5e-13);

  auto cancelling_real = [](const double energy) { return 1.0 + energy; };
  auto cancelling_imaginary = [](const double energy) { return -(1.0 + energy); };
  NRG::Hilb::integrator cancelling_integration(NRG::Hilb::integrator::configured, 64,
                                               NRG::Tools::QagRule::gauss61);
  const auto cancelling_result = NRG::Hilb::hilbert_transform(
    cancelling_integration, cancelling_real, cancelling_imaginary, 1.0,
    std::complex<double>{2e15, 2e15}, 0.1, 0, 1e-40, 1e-12);
  EXPECT_NEAR(cancelling_result.real() / (-8.333333333333333e-32), 1.0, 5e-13);

  auto extreme_cancelling_real = [](const double energy) { return 5e99 * (1.0 + energy); };
  auto extreme_cancelling_imaginary = [](const double energy) { return -5e99 * (1.0 + energy); };
  NRG::Hilb::integrator extreme_cancelling_integration(NRG::Hilb::integrator::configured, 64,
                                                       NRG::Tools::QagRule::gauss61);
  const auto extreme_cancelling_result = NRG::Hilb::hilbert_transform(
    extreme_cancelling_integration, extreme_cancelling_real, extreme_cancelling_imaginary, 1.0,
    std::complex<double>{1e100, 1e100}, 0.1, 0, 1e-115, 1e-12);
  EXPECT_NEAR(extreme_cancelling_result.real() / (-1.0 / 6e100), 1.0, 5e-13);
  EXPECT_NEAR(extreme_cancelling_result.imag(), -1.0, 5e-13);

  constexpr double tiny_bandwidth = 1e-200;
  constexpr double huge_bandwidth = 1e308;
  const auto minimum_imaginary = NRG::Hilb::minimum_safe_imaginary_part();
  auto zero_density = [](const double) { return 0.0; };
  auto tiny_flat = [](const double) { return 0.5 / tiny_bandwidth; };
  auto huge_flat = [](const double) { return 0.5 / huge_bandwidth; };
  NRG::Hilb::integrator extreme_scale_integration(NRG::Hilb::integrator::configured, 64,
                                                   NRG::Tools::QagRule::gauss61);

  const auto tiny_direct = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, tiny_flat, zero_density, tiny_bandwidth,
    std::complex<double>{0.0, minimum_imaginary}, 0.0, 0, 0.0, 1e-12);
  EXPECT_DOUBLE_EQ(tiny_direct.real(), 0.0);
  EXPECT_NEAR(tiny_direct.imag() / (-std::atan(tiny_bandwidth / minimum_imaginary) / tiny_bandwidth), 1.0,
              3e-14);

  const auto huge_direct = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, huge_flat, zero_density, huge_bandwidth,
    std::complex<double>{0.0, huge_bandwidth}, 0.0, 0, 0.0, 1e-12);
  EXPECT_DOUBLE_EQ(huge_direct.real(), 0.0);
  EXPECT_NEAR(huge_direct.imag() / (-std::atan(1.0) / huge_bandwidth), 1.0, 3e-14);

  const auto huge_subtracted = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, huge_flat, zero_density, huge_bandwidth,
    std::complex<double>{0.0, minimum_imaginary}, 1e-3, 0, 0.0, 1e-12);
  EXPECT_DOUBLE_EQ(huge_subtracted.real(), 0.0);
  EXPECT_NEAR(huge_subtracted.imag() / (-std::atan(huge_bandwidth / minimum_imaginary) / huge_bandwidth), 1.0,
              3e-14);

  auto huge_constant = [](const double) { return 1e307; };
  const auto asymmetric_huge = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, huge_constant, zero_density, std::numeric_limits<double>::max(),
    std::complex<double>{1.0, minimum_imaginary}, 1e-3, 0, 0.0, 1e-12);
  EXPECT_NEAR(asymmetric_huge.real(), 0.11125369292536008, 3e-15);

  auto half_density = [](const double) { return 0.5; };
  const auto lower_edge = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, half_density, zero_density, 1.0,
    std::complex<double>{-1.0, minimum_imaginary}, 1e-3, 0, 0.0, 1e-12);
  EXPECT_NEAR(lower_edge.real(), -177.445678223346, 3e-13);
  EXPECT_NEAR(lower_edge.imag(), -0.7853981633974483, 3e-15);

  auto sampled_zero_imaginary = [](const double energy) {
    const auto square = energy * energy;
    return square * (square - 0.25) * (square - 1.0);
  };
  const auto sampled_zero = NRG::Hilb::hilbert_transform(
    extreme_scale_integration, large_density, sampled_zero_imaginary, 1.0,
    std::complex<double>{0.0, 1.0}, 0.0, 0, 0.0, 1e-12);
  EXPECT_NEAR(sampled_zero.real(), -0.026990816987241543, 3e-15);
}

TEST(Hilb, tabulated_density_is_weighted_after_interpolation) { // NOLINT
  const std::vector<double> energies{-1.0, -0.5, 0.0, 0.5, 1.0};
  const std::vector<double> density(energies.size(), 0.5);
  const std::vector<double> zero(energies.size(), 0.0);
  const std::complex<double> z{-0.3, 0.1};
  const auto expected = z * z * flat_band_h0(z) - z;

  const auto result = NRG::Hilb::hilbert_transform(energies, density, zero, z, 1e-3, 2);

  expect_complex_near(result, expected, 1e-10);
}

TEST(Hilb, piecewise_polynomial_support_must_fit_bandwidth) { // NOLINT
  const NRG::Tools::PiecewisePolynomial<double> inside{{-1.0, 1.0}, {{0.5}}};
  const NRG::Tools::PiecewisePolynomial<double> below{{-1.1, 0.5}, {{0.5}}};
  const NRG::Tools::PiecewisePolynomial<double> above{{-0.5, 1.1}, {{0.5}}};
  const std::complex<double> z{0.2, 0.4};

  EXPECT_NO_THROW(NRG::Hilb::hilbert_transform(inside, B, z));
  EXPECT_THROW(NRG::Hilb::hilbert_transform(below, B, z), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(above, B, z), std::invalid_argument);
}

TEST(Hilb, tabulated_transform_accepts_explicit_interpolation_and_defaults_to_cubic) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const std::vector<double> energies{-1.0, -0.5, 0.0, 0.5, 1.0};
  const std::vector<double> real_density{0.0, 0.2, 1.0, 0.1, 0.0};
  const std::vector<double> imaginary_density{0.0, -0.1, 0.3, -0.2, 0.0};
  const std::complex<double> z{0.2, 0.4};

  constexpr std::array methods{InterpolationMethod::linear, InterpolationMethod::cspline, InterpolationMethod::akima,
                               InterpolationMethod::steffen};
  std::array<std::complex<double>, methods.size()> results;
  for (std::size_t index = 0; index < methods.size(); ++index) {
    const auto method = methods[index];
    NRG::Hilb::interpolator rhor(energies, real_density, 0.0, method);
    NRG::Hilb::interpolator rhoi(energies, imaginary_density, 0.0, method);
    NRG::Hilb::integrator integration;
    results[index] = NRG::Hilb::hilbert_transform_with_interpolation(energies, real_density, imaginary_density, z, method);
    EXPECT_EQ(results[index], NRG::Hilb::hilbert_transform(integration, rhor, rhoi, B, z));
    EXPECT_EQ(results[index], NRG::Hilb::hilbert_transform(integration,
                                                           static_cast<const NRG::Hilb::interpolator &>(rhor),
                                                           static_cast<const NRG::Hilb::interpolator &>(rhoi), B, z));
    EXPECT_EQ(results[index], NRG::Hilb::hilbert_transform(integration, rhor,
                                                           static_cast<const NRG::Hilb::interpolator &>(rhoi), B, z));
    EXPECT_EQ(results[index], NRG::Hilb::hilbert_transform(integration,
                                                           static_cast<const NRG::Hilb::interpolator &>(rhor), rhoi, B, z));
  }

  const auto legacy = NRG::Hilb::hilbert_transform(energies, real_density, imaginary_density, z);
  EXPECT_EQ(legacy, results[1]);
  EXPECT_GT(std::abs(results[0] - results[1]), 1e-3);
  EXPECT_GT(std::abs(results[2] - results[1]), 1e-3);
  EXPECT_GT(std::abs(results[3] - results[1]), 1e-3);
}

TEST(Hilb, zero_energy_power_preserves_transform) { // NOLINT
  auto rho = [](const double energy) { return 0.5 + 0.1 * energy; };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{0.2, 0.1};

  const auto implicit_zero = NRG::Hilb::hilbert_transform(rho, zero, B, z);
  const auto explicit_zero = NRG::Hilb::hilbert_transform(rho, zero, B, z, 1e-3, 0);

  EXPECT_EQ(implicit_zero, explicit_zero);
}

TEST(Hilb, energy_power_parser_requires_nonnegative_integer) { // NOLINT
  auto rho = [](double) { return 1.0; };
  auto zero = [](double) { return 0.0; };
  EXPECT_EQ(NRG::Hilb::parse_energy_power("0"), 0);
  EXPECT_EQ(NRG::Hilb::parse_energy_power("12"), 12);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("-1"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("1.5"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("2junk"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("999999999999999999999999"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, std::complex<double>{0.0, 0.1}, 1e-3, -1), std::invalid_argument);
}

TEST(Hilb, minimum_safe_imaginary_part_avoids_underflow) { // NOLINT
  auto rho = [](const double energy) { return 0.5 + 0.1 * energy; };
  auto zero = [](const double) { return 0.0; };
  const auto minimum = NRG::Hilb::minimum_safe_imaginary_part();

  EXPECT_DOUBLE_EQ(minimum, std::sqrt(std::numeric_limits<double>::min()));
  const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, std::complex<double>{0.2, minimum});
  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, std::complex<double>{0.2, minimum / 2.0}), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, std::complex<double>{0.2, 0.0}), std::invalid_argument);
}

TEST(Hilb, integration_settings_are_configurable_and_validated) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{0.2, 0.1};
  const auto expected = flat_band_h0(z);

  const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, z, 0.2, 0, 1e-12, 1e-8);
  expect_complex_near(result, expected, 1e-8);
  EXPECT_NO_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, z, 0.2, 0, 0.0, 1e-10));
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, z, -1.0, 0, 1e-12, 1e-8), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, z, 0.2, 0, -1.0, 1e-8), std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, z, 0.2, 0, 0.0, 0.0), std::invalid_argument);
}

TEST(Hilb, integration_workspace_is_reusable_and_copyable) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  NRG::Hilb::integrator integration;
  const std::complex<double> first_z{0.2, 0.1};
  const std::complex<double> second_z{-0.3, 0.2};
  expect_complex_near(NRG::Hilb::hilbert_transform(integration, rho, zero, B, first_z), flat_band_h0(first_z), 1e-10);
  expect_complex_near(NRG::Hilb::hilbert_transform(integration, rho, zero, B, second_z), flat_band_h0(second_z), 1e-10);

  auto copied = integration;
  expect_complex_near(NRG::Hilb::hilbert_transform(copied, rho, zero, B, first_z), flat_band_h0(first_z), 1e-10);
  NRG::Hilb::integrator assigned(10);
  assigned = integration;
  expect_complex_near(NRG::Hilb::hilbert_transform(assigned, rho, zero, B, second_z), flat_band_h0(second_z), 1e-10);
}

TEST(Hilb, integrator_accepts_configured_workspace_rule_and_policy) { // NOLINT
  using NRG::Tools::GslErrorPolicy;
  using NRG::Tools::QagRule;

  EXPECT_THROW(NRG::Hilb::integrator(NRG::Hilb::integrator::configured, 0, QagRule::gauss61, GslErrorPolicy::fail),
               std::invalid_argument);
  EXPECT_THROW(NRG::Hilb::integrator(NRG::Hilb::integrator::configured,
                                     NRG::Tools::qag_workspace_limit_maximum() + 1,
                                     QagRule::gauss61, GslErrorPolicy::fail),
               std::invalid_argument);

  NRG::Hilb::integrator integration(NRG::Hilb::integrator::configured, 32, QagRule::gauss61, GslErrorPolicy::fail);
  std::size_t gauss61_evaluations = 0;
  const auto result = integration([&gauss61_evaluations](const double x) {
    ++gauss61_evaluations;
    return x * x * x * x;
  }, -1.0, 1.0, 1e-12, 1e-12);
  EXPECT_NEAR(result, 0.4, 1e-14);
  EXPECT_EQ(gauss61_evaluations, 61);

  NRG::Hilb::integrator gauss15(NRG::Hilb::integrator::configured, 32, QagRule::gauss15, GslErrorPolicy::fail);
  std::size_t gauss15_evaluations = 0;
  EXPECT_NEAR(gauss15([&gauss15_evaluations](const double x) {
    ++gauss15_evaluations;
    return x * x * x * x;
  }, -1.0, 1.0, 1e-12, 1e-12), 0.4, 1e-14);
  EXPECT_EQ(gauss15_evaluations, 15);

  auto copied = integration;
  EXPECT_NEAR(copied([](const double x) { return x * x; }, 0.0, 1.0, 1e-12, 1e-12), 1.0 / 3.0, 1e-14);
}

TEST(Hilb, configured_integrator_policies_classify_nonfinite_results) { // NOLINT
  using NRG::Tools::GslErrorPolicy;
  using NRG::Tools::QagRule;
  gsl_set_error_handler_off();
  const auto nonfinite = [](const double) { return std::numeric_limits<double>::quiet_NaN(); };

  NRG::Hilb::gsl_failure_summary ignored_failures;
  NRG::Hilb::integrator ignored(NRG::Hilb::integrator::configured, 32, QagRule::gauss21, GslErrorPolicy::ignore,
                                &ignored_failures);
  EXPECT_FALSE(std::isfinite(ignored(nonfinite, -1.0, 1.0)));
  EXPECT_EQ(ignored_failures.count(), 0);

  NRG::Hilb::gsl_failure_summary warned_failures;
  NRG::Hilb::integrator warned(NRG::Hilb::integrator::configured, 32, QagRule::gauss31, GslErrorPolicy::warn,
                               &warned_failures);
  EXPECT_FALSE(std::isfinite(warned(nonfinite, -1.0, 1.0)));
  EXPECT_EQ(warned_failures.count(), 1);

  NRG::Hilb::integrator failing(NRG::Hilb::integrator::configured, 32, QagRule::gauss41, GslErrorPolicy::fail);
  auto *previous_handler = gsl_set_error_handler(&count_gsl_handler_calls);
  gsl_handler_calls = 0;
  EXPECT_THROW(failing(nonfinite, -1.0, 1.0), std::runtime_error);
  if (previous_handler)
    gsl_set_error_handler(previous_handler);
  else
    gsl_set_error_handler_off();
  EXPECT_EQ(gsl_handler_calls, 0);
}

TEST(Hilb, gsl_failures_are_aggregated_and_reported_once) { // NOLINT
  gsl_set_error_handler_off();
  NRG::Hilb::gsl_failure_summary failures;
  NRG::Hilb::integrator integration(1000, false, &failures);
  const auto result = integration([](const double) { return std::numeric_limits<double>::quiet_NaN(); }, -1.0, 1.0);
  EXPECT_FALSE(std::isfinite(result));
  EXPECT_EQ(failures.count(), 1);

  std::ostringstream report;
  failures.report(report);
  failures.report(report);
  const auto report_output = report.str();
  EXPECT_EQ(count_occurrences(report_output, "hilb: warning:"), 1);
  EXPECT_EQ(std::count(report_output.begin(), report_output.end(), '\n'), 1);
}

TEST(Hilb, numeric_row_reader_enforces_physical_records) { // NOLINT
  std::istringstream valid("\n  1 2 3  \n\t\n4 5 6\r\n");
  NRG::Hilb::numeric_row_reader<3> rows(valid, "points.dat");
  const auto first = rows.next();
  ASSERT_TRUE(first);
  EXPECT_EQ(first->line_number, 2);
  EXPECT_EQ(first->values, (std::array<double, 3>{1.0, 2.0, 3.0}));
  const auto second = rows.next();
  ASSERT_TRUE(second);
  EXPECT_EQ(second->line_number, 4);
  EXPECT_FALSE(rows.next());

  std::istringstream extra("1 2 3 4\n");
  NRG::Hilb::numeric_row_reader<3> extra_rows(extra, "points.dat");
  EXPECT_THROW(extra_rows.next(), std::runtime_error);
  std::istringstream split("1 2\n3\n");
  NRG::Hilb::numeric_row_reader<3> split_rows(split, "points.dat");
  EXPECT_THROW(split_rows.next(), std::runtime_error);
  std::istringstream comment("# comment\n");
  NRG::Hilb::numeric_row_reader<2> comment_rows(comment, "DOS.dat");
  EXPECT_THROW(comment_rows.next(), std::runtime_error);
  std::istringstream nonfinite("1 nan\n");
  NRG::Hilb::numeric_row_reader<2> nonfinite_rows(nonfinite, "DOS.dat");
  EXPECT_THROW(nonfinite_rows.next(), std::runtime_error);
}

TEST(Hilb, dmft_clipping_is_reported_once) { // NOLINT
  const std::string resigma = "hilb_clipping_resigma.dat";
  const std::string imsigma = "hilb_clipping_imsigma.dat";
  const std::string reoutput = "hilb_clipping_reoutput.dat";
  const std::string imoutput = "hilb_clipping_imoutput.dat";
  write_file(resigma, "0 0\n0.25 0\n0.5 0\n");
  write_file(imsigma, "0 0\n0.25 -0.05\n0.5 -0.2\n");

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(run_hilb({"hilb", "-G", "-c", "0.1", resigma, imsigma, reoutput, imoutput}));
  const auto stderr_output = testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();

  EXPECT_NE(stderr_output.find("hilb: clipped 2 ImSigma value(s)"), std::string::npos);
  EXPECT_NE(stderr_output.find("first at omega=0: 0 -> -0.10000000000000001"), std::string::npos);
  EXPECT_EQ(count_occurrences(stderr_output, "hilb: clipped"), 1);

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, default_dmft_clipping_produces_a_finite_result) { // NOLINT
  const std::string resigma = "hilb_safe_clipping_resigma.dat";
  const std::string imsigma = "hilb_safe_clipping_imsigma.dat";
  const std::string reoutput = "hilb_safe_clipping_reoutput.dat";
  const std::string imoutput = "hilb_safe_clipping_imoutput.dat";
  write_file(resigma, "0 0\n");
  write_file(imsigma, "0 0\n");

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}));
  const auto stderr_output = testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_NE(stderr_output.find("hilb: clipped 1 ImSigma value(s)"), std::string::npos);
  EXPECT_NE(stderr_output.find("-1.4916681462400413e-154"), std::string::npos);

  double label = 0.0;
  double value = 0.0;
  std::ifstream refile(reoutput);
  EXPECT_TRUE(static_cast<bool>(refile >> label >> value));
  EXPECT_TRUE(std::isfinite(value));
  std::ifstream imfile(imoutput);
  EXPECT_TRUE(static_cast<bool>(imfile >> label >> value));
  EXPECT_TRUE(std::isfinite(value));

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, dmft_without_clipping_has_no_clipping_report) { // NOLINT
  const std::string resigma = "hilb_default_clipping_resigma.dat";
  const std::string imsigma = "hilb_default_clipping_imsigma.dat";
  const std::string reoutput = "hilb_default_clipping_reoutput.dat";
  const std::string imoutput = "hilb_default_clipping_imoutput.dat";
  write_file(resigma, "0 0\n");
  write_file(imsigma, "0 -0.1\n");

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}));
  const auto stderr_output = testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_EQ(stderr_output.find("hilb: clipped"), std::string::npos);

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, frequency_tolerance_is_configurable) { // NOLINT
  const std::string resigma = "hilb_frequency_resigma.dat";
  const std::string imsigma = "hilb_frequency_imsigma.dat";
  const std::string reoutput = "hilb_frequency_reoutput.dat";
  const std::string imoutput = "hilb_frequency_imoutput.dat";
  write_file(resigma, "0 0\n");
  write_file(imsigma, "0.00005 -0.1\n");

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(run_hilb({"hilb", "-G", "-f", "0.0001", resigma, imsigma, reoutput, imoutput}));
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, numerical_cli_options_reject_invalid_values) { // NOLINT
  EXPECT_THROW(run_hilb({"hilb", "-c", "1e-200", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-c", "invalid", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-t", "-1", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-a", "-1", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-r", "-1", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-f", "-1", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-a", "0", "-r", "0", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-t", "nan", "0", "0.1"}), std::runtime_error);
}

TEST(Hilb, numerical_cli_options_are_accepted) { // NOLINT
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(run_hilb({"hilb", "-t", "0.2", "-a", "1e-12", "-r", "1e-8", "0.2", "0.1"}));
  testing::internal::GetCapturedStdout();
}

TEST(Hilb, shared_gsl_cli_options_are_accepted) { // NOLINT
  const auto output = run_hilb_captured({"hilb", "--workspace-limit", "64", "--quadrature-rule", "61",
                                          "--gsl-error-policy", "fail", "--epsabs", "1e-10", "--epsrel", "1e-8",
                                          "-G", "0.2", "0.4"});
  EXPECT_TRUE(output.err.empty());
  const auto result = parse_complex_output(output.out);
  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));

  const auto oversized_workspace = std::to_string(NRG::Tools::qag_workspace_limit_maximum() + 1);
  EXPECT_THROW(run_hilb({"hilb", "--workspace-limit", oversized_workspace, "0.2", "0.4"}),
               std::runtime_error);
}

TEST(Hilb, built_in_dos_cli_qag_warning_and_failure_policies_are_wired) { // NOLINT
  const std::vector<std::string> difficult_options{"--workspace-limit", "1", "--quadrature-rule", "15",
                                                   "--epsabs", "0", "--epsrel", "2e-14"};
  auto warning_arguments = std::vector<std::string>{"hilb"};
  warning_arguments.insert(warning_arguments.end(), difficult_options.begin(), difficult_options.end());
  warning_arguments.insert(warning_arguments.end(), {"--gsl-error-policy", "warn", "0.2", "0.4"});
  const auto warning = run_hilb_captured(std::move(warning_arguments));
  EXPECT_TRUE(std::isfinite(std::stod(warning.out)));
  EXPECT_NE(warning.err.find("hilb: warning: 2 GSL integration call(s) reported failure"), std::string::npos);
  EXPECT_NE(warning.err.find("exceeded max number of iterations"), std::string::npos);

  auto failure_arguments = std::vector<std::string>{"hilb"};
  failure_arguments.insert(failure_arguments.end(), difficult_options.begin(), difficult_options.end());
  failure_arguments.insert(failure_arguments.end(), {"--gsl-error-policy", "fail", "0.2", "0.4"});
  EXPECT_THROW(run_hilb(std::move(failure_arguments)), std::runtime_error);
}

TEST(Hilb, shared_gsl_cli_options_reject_invalid_values) { // NOLINT
  const auto overflow = std::to_string(std::numeric_limits<std::size_t>::max()) + "0";
  EXPECT_ANY_THROW(run_hilb({"hilb", "--interpolation", "cubic", "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--workspace-limit", "0", "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--workspace-limit", overflow, "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--quadrature-rule", "17", "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--gsl-error-policy", "error", "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--epsabs", "-1", "0", "0.1"}));
  EXPECT_ANY_THROW(run_hilb({"hilb", "--epsrel", "-1", "0", "0.1"}));
}

TEST(Hilb, cli_legacy_interpolation_default_is_cubic) { // NOLINT
  const std::string dos = "hilb_default_interpolation_dos.dat";
  write_file(dos, "-1 0\n-0.5 0.2\n0 1\n0.5 0.1\n1 0\n");

  const auto legacy = run_hilb_captured({"hilb", "-d", dos, "0.2", "0.4"});
  const auto cubic = run_hilb_captured({"hilb", "-i", "cspline", "-d", dos, "0.2", "0.4"});
  const auto linear = run_hilb_captured({"hilb", "-i", "linear", "-d", dos, "0.2", "0.4"});

  EXPECT_TRUE(legacy.err.empty());
  EXPECT_TRUE(cubic.err.empty());
  EXPECT_TRUE(linear.err.empty());
  EXPECT_EQ(legacy.out, cubic.out);
  EXPECT_GT(std::abs(std::stod(legacy.out) - std::stod(linear.out)), 1e-3);

  std::remove(dos.c_str());
}

TEST(Hilb, cli_argument_shifts_and_bandwidth_scaling_are_effective) { // NOLINT
  const std::string dos = "hilb_shift_dos.dat";
  write_file(dos, "-1 0.25\n0 1\n1 0.25\n");

  const auto shifted = run_hilb_captured({"hilb", "-d", dos, "-G", "-x", "0.125", "-y", "-0.125",
                                          "0.25", "0.5"});
  const auto translated = run_hilb_captured({"hilb", "-d", dos, "-G", "0.375", "0.375"});
  EXPECT_TRUE(shifted.err.empty());
  EXPECT_EQ(shifted.out, translated.out);

  const auto scaled = run_hilb_captured({"hilb", "-G", "-s", "2", "0.2", "0.4"});
  const auto bandwidth = run_hilb_captured({"hilb", "-G", "-B", "0.5", "0.2", "0.4"});
  const auto unscaled = run_hilb_captured({"hilb", "-G", "0.2", "0.4"});
  EXPECT_TRUE(scaled.err.empty());
  EXPECT_EQ(scaled.out, bandwidth.out);
  EXPECT_GT(std::abs(parse_complex_output(scaled.out) - parse_complex_output(unscaled.out)), 0.1);

  std::remove(dos.c_str());
}

TEST(Hilb, dmft_shifts_and_clipping_match_equivalent_arguments) { // NOLINT
  const std::string dos = "hilb_dmft_options_dos.dat";
  const std::string shifted_real = "hilb_dmft_shifted_real.dat";
  const std::string shifted_imaginary = "hilb_dmft_shifted_imaginary.dat";
  const std::string translated_real = "hilb_dmft_translated_real.dat";
  const std::string translated_imaginary = "hilb_dmft_translated_imaginary.dat";
  const std::string clipped_imaginary = "hilb_dmft_clipped_imaginary.dat";
  const std::string shifted_reoutput = "hilb_dmft_shifted_reoutput.dat";
  const std::string shifted_imoutput = "hilb_dmft_shifted_imoutput.dat";
  const std::string translated_reoutput = "hilb_dmft_translated_reoutput.dat";
  const std::string translated_imoutput = "hilb_dmft_translated_imoutput.dat";
  const std::string clipped_reoutput = "hilb_dmft_clipped_reoutput.dat";
  const std::string clipped_imoutput = "hilb_dmft_clipped_imoutput.dat";
  const std::array files{dos, shifted_real, shifted_imaginary, translated_real, translated_imaginary,
                         clipped_imaginary, shifted_reoutput, shifted_imoutput, translated_reoutput,
                         translated_imoutput, clipped_reoutput, clipped_imoutput};
  for (const auto &file : files) std::remove(file.c_str());
  write_file(dos, "-1 0.25\n0 1\n1 0.25\n");
  write_file(shifted_real, "0.375 -0.125\n");
  write_file(shifted_imaginary, "0.375 -0.25\n");
  write_file(translated_real, "0.375 -0.25\n");
  write_file(translated_imaginary, "0.375 -0.375\n");
  write_file(clipped_imaginary, "0.375 0\n");

  const auto shifted = run_hilb_captured({"hilb", "-G", "-d", dos, "-x", "0.125", "-y", "0.125",
                                          shifted_real, shifted_imaginary, shifted_reoutput, shifted_imoutput});
  const auto translated = run_hilb_captured({"hilb", "-G", "-d", dos, translated_real, translated_imaginary,
                                             translated_reoutput, translated_imoutput});
  EXPECT_TRUE(shifted.err.empty());
  EXPECT_TRUE(translated.err.empty());
  EXPECT_EQ(read_file(shifted_reoutput), read_file(translated_reoutput));
  EXPECT_EQ(read_file(shifted_imoutput), read_file(translated_imoutput));

  const auto clipped = run_hilb_captured({"hilb", "-G", "-d", dos, "-c", "0.1", shifted_real,
                                          clipped_imaginary, clipped_reoutput, clipped_imoutput});
  EXPECT_NE(clipped.err.find("hilb: clipped 1 ImSigma value(s)"), std::string::npos);
  write_file(translated_imaginary, "0.375 -0.1\n");
  const auto unclipped = run_hilb_captured({"hilb", "-G", "-d", dos, shifted_real, translated_imaginary,
                                            translated_reoutput, translated_imoutput});
  EXPECT_TRUE(unclipped.err.empty());
  EXPECT_EQ(read_file(clipped_reoutput), read_file(translated_reoutput));
  EXPECT_EQ(read_file(clipped_imoutput), read_file(translated_imoutput));

  for (const auto &file : files) std::remove(file.c_str());
}

TEST(Hilb, invalid_arity_is_rejected) { // NOLINT
  EXPECT_THROW(run_hilb({"hilb"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "1", "2", "3"}), std::runtime_error);
}

TEST(Hilb, legacy_numeric_options_are_strict) { // NOLINT
  EXPECT_THROW(run_hilb({"hilb", "junk", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-B", "0", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-s", "-1", "0", "0.1"}), std::runtime_error);
  EXPECT_THROW(run_hilb({"hilb", "-x", "nan", "0", "0.1"}), std::runtime_error);
  EXPECT_EQ(NRG::Hilb::OUTPUT_PRECISION, std::numeric_limits<double>::max_digits10);
}

TEST(Hilb, dmft_inputs_must_have_equal_record_counts) { // NOLINT
  const std::string resigma = "hilb_unequal_resigma.dat";
  const std::string imsigma = "hilb_unequal_imsigma.dat";
  const std::string reoutput = "hilb_unequal_reoutput.dat";
  const std::string imoutput = "hilb_unequal_imoutput.dat";
  const std::string preserved = "preserve existing output\n";
  write_file(resigma, "0 0\n1 0\n");
  write_file(imsigma, "0 -0.1\n");
  write_file(reoutput, preserved);
  write_file(imoutput, preserved);

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_EQ(read_file(reoutput), preserved);
  EXPECT_EQ(read_file(imoutput), preserved);

  write_file(resigma, "0 0\n");
  write_file(imsigma, "0 -0.1\n1 -0.1\n");
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_EQ(read_file(reoutput), preserved);
  EXPECT_EQ(read_file(imoutput), preserved);

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, existing_outputs_survive_parse_and_calculation_failures) { // NOLINT
  const std::string input = "hilb_preserved_stream_input.dat";
  const std::string stream_output = "hilb_preserved_stream_output.dat";
  const std::string single_output = "hilb_preserved_single_output.dat";
  const std::string preserved = "preserve existing output\n";
  write_file(input, "0 0 0.1\nmalformed row\n");
  write_file(stream_output, preserved);
  write_file(single_output, preserved);

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-o", stream_output, input}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_EQ(read_file(stream_output), preserved);

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-o", single_output, "0", "0"}), std::invalid_argument);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();
  EXPECT_EQ(read_file(single_output), preserved);

  std::remove(input.c_str());
  std::remove(stream_output.c_str());
  std::remove(single_output.c_str());
}

TEST(Hilb, output_file_must_not_alias_stream_input) { // NOLINT
  const std::string input = "hilb_collision_input.dat";
  const std::string symlink_alias = "hilb_collision_input_symlink.dat";
  const std::string hardlink_alias = "hilb_collision_input_hardlink.dat";
  const std::string contents = "1 0 0.1\n";
  std::remove(symlink_alias.c_str());
  std::remove(hardlink_alias.c_str());
  write_file(input, contents);
  ASSERT_EQ(::symlink(input.c_str(), symlink_alias.c_str()), 0);
  ASSERT_EQ(::link(input.c_str(), hardlink_alias.c_str()), 0);

  for (const auto &output : {input, "./" + input, symlink_alias, hardlink_alias}) {
    EXPECT_THROW(run_hilb({"hilb", "-o", output, input}), std::runtime_error);
    EXPECT_EQ(read_file(input), contents);
  }

  std::remove(symlink_alias.c_str());
  std::remove(hardlink_alias.c_str());
  std::remove(input.c_str());
}

TEST(Hilb, output_write_failures_are_reported) { // NOLINT
  if (!std::filesystem::exists("/dev/full")) GTEST_SKIP() << "/dev/full is unavailable";
  EXPECT_THROW(run_hilb({"hilb", "-o", "/dev/full", "0", "0.1"}), std::runtime_error);
}

TEST(Hilb, dmft_output_collisions_are_filesystem_aware) { // NOLINT
  const std::string resigma = "hilb_collision_resigma.dat";
  const std::string imsigma = "hilb_collision_imsigma.dat";
  const std::string reoutput = "hilb_collision_reoutput.dat";
  const std::string imoutput = "hilb_collision_imoutput.dat";
  const std::string dos = "hilb_collision_dos.dat";
  const std::string real_contents = "0 0\n";
  const std::string output_contents = "must remain intact\n";
  const std::string dos_contents = "-1 0\n0 1\n1 0\n";
  for (const auto &file : {resigma, imsigma, reoutput, imoutput, dos}) std::remove(file.c_str());
  write_file(resigma, real_contents);
  write_file(imsigma, "0 -0.1\n");

  ASSERT_EQ(::link(resigma.c_str(), reoutput.c_str()), 0);
  EXPECT_THROW(run_hilb({"hilb", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  EXPECT_EQ(read_file(resigma), real_contents);
  std::remove(reoutput.c_str());

  write_file(reoutput, output_contents);
  ASSERT_EQ(::symlink(reoutput.c_str(), imoutput.c_str()), 0);
  EXPECT_THROW(run_hilb({"hilb", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  EXPECT_EQ(read_file(reoutput), output_contents);
  std::remove(imoutput.c_str());
  std::remove(reoutput.c_str());

  EXPECT_THROW(run_hilb({"hilb", resigma, imsigma, reoutput, "./" + reoutput}), std::runtime_error);
  EXPECT_FALSE(std::filesystem::exists(reoutput));

  write_file(dos, dos_contents);
  ASSERT_EQ(::link(dos.c_str(), reoutput.c_str()), 0);
  EXPECT_THROW(run_hilb({"hilb", "-d", dos, resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  EXPECT_EQ(read_file(dos), dos_contents);

  for (const auto &file : {resigma, imsigma, reoutput, imoutput, dos}) std::remove(file.c_str());
}

TEST(Hilb, output_file_must_not_alias_tabulated_dos) { // NOLINT
  const std::string dos = "hilb_output_alias_dos.dat";
  const std::string symlink_alias = "hilb_output_alias_dos_symlink.dat";
  const std::string hardlink_alias = "hilb_output_alias_dos_hardlink.dat";
  const std::string contents = "-1 0\n0 1\n1 0\n";
  std::remove(symlink_alias.c_str());
  std::remove(hardlink_alias.c_str());
  write_file(dos, contents);
  ASSERT_EQ(::symlink(dos.c_str(), symlink_alias.c_str()), 0);
  ASSERT_EQ(::link(dos.c_str(), hardlink_alias.c_str()), 0);

  for (const auto &output : {dos, "./" + dos, symlink_alias, hardlink_alias}) {
    EXPECT_THROW(run_hilb({"hilb", "-d", dos, "-o", output, "0", "0.1"}), std::runtime_error);
    EXPECT_EQ(read_file(dos), contents);
  }

  std::remove(symlink_alias.c_str());
  std::remove(hardlink_alias.c_str());
  std::remove(dos.c_str());
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
