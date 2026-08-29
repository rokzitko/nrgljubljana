#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
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

TEST(Hilb, tabulated_density_is_weighted_after_interpolation) { // NOLINT
  const std::vector<double> energies{-1.0, -0.5, 0.0, 0.5, 1.0};
  const std::vector<double> density(energies.size(), 0.5);
  const std::vector<double> zero(energies.size(), 0.0);
  const std::complex<double> z{-0.3, 0.1};
  const auto expected = z * z * flat_band_h0(z) - z;

  const auto result = NRG::Hilb::hilbert_transform(energies, density, zero, z, 1e-3, 2);

  expect_complex_near(result, expected, 1e-10);
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
  const std::string dos = "hilb_gsl_options_dos.dat";
  write_file(dos, "-1 0\n-0.5 0.5\n0 1\n0.5 0.5\n1 0\n");

  const auto output = run_hilb_captured({"hilb", "--interpolation", "steffen", "--workspace-limit", "64",
                                         "--quadrature-rule", "61", "--gsl-error-policy", "fail", "--epsabs", "1e-10",
                                         "--epsrel", "1e-8", "-d", dos, "0.2", "0.4"});
  EXPECT_TRUE(output.err.empty());
  EXPECT_TRUE(std::isfinite(std::stod(output.out)));

  const auto oversized_workspace = std::to_string(NRG::Tools::qag_workspace_limit_maximum() + 1);
  EXPECT_THROW(run_hilb({"hilb", "-d", dos, "--workspace-limit", oversized_workspace, "0.2", "0.4"}),
               std::runtime_error);

  std::remove(dos.c_str());
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

TEST(Hilb, version_and_invalid_arity_have_distinct_behavior) { // NOLINT
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(run_hilb({"hilb", "-V"}));
  EXPECT_EQ(testing::internal::GetCapturedStdout(), "hilb 2026.09\n");
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
  write_file(resigma, "0 0\n1 0\n");
  write_file(imsigma, "0 -0.1\n");

  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();

  write_file(resigma, "0 0\n");
  write_file(imsigma, "0 -0.1\n1 -0.1\n");
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  EXPECT_THROW(run_hilb({"hilb", "-G", resigma, imsigma, reoutput, imoutput}), std::runtime_error);
  testing::internal::GetCapturedStderr();
  testing::internal::GetCapturedStdout();

  std::remove(resigma.c_str());
  std::remove(imsigma.c_str());
  std::remove(reoutput.c_str());
  std::remove(imoutput.c_str());
}

TEST(Hilb, output_file_must_not_alias_tabulated_dos) { // NOLINT
  const std::string dos = "hilb_output_alias_dos.dat";
  const std::string contents = "-1 0\n0 1\n1 0\n";
  write_file(dos, contents);

  EXPECT_THROW(run_hilb({"hilb", "-d", dos, "-o", dos, "0", "0.1"}), std::runtime_error);

  std::ifstream file(dos);
  std::ostringstream preserved;
  preserved << file.rdbuf();
  EXPECT_EQ(preserved.str(), contents);
  std::remove(dos.c_str());
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
