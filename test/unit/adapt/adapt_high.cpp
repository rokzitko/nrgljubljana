#include <gtest/gtest.h>

#include <cstddef>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include <adapt/adapt.hpp>

using namespace NRG::Adapt;

namespace {

void write_file(const std::string &filename, const std::string &contents) {
  std::ofstream file(filename);
  file << contents;
}

} // namespace

TEST(Adapt, parser_skips_blank_lines) { // NOLINT
  const auto filename = "adapt_blank_lines.param";
  write_file(filename, "[param]\n\nLambda=3\n\n# comment\n\nxmax=7\n");

  Params P(filename);
  EXPECT_EQ(P.P("Lambda", 0.0), 3.0);
  EXPECT_EQ(P.P("xmax", 0.0), 7.0);

  std::remove(filename);
}

TEST(Adapt, parser_accepts_extended_bool_values_and_rejects_invalid_values) { // NOLINT
  const auto filename = "adapt_bool.param";
  write_file(filename,
             "[param]\n"
             "yes = YeS\n"
             "true = tRuE\n"
             "one = 1\n"
             "no = nO\n"
             "false = FaLsE\n"
             "zero = 0\n"
             "invalid = maybe\n");

  Params P(filename);
  EXPECT_TRUE(P.Pbool("yes", false));
  EXPECT_TRUE(P.Pbool("true", false));
  EXPECT_TRUE(P.Pbool("one", false));
  EXPECT_FALSE(P.Pbool("no", true));
  EXPECT_FALSE(P.Pbool("false", true));
  EXPECT_FALSE(P.Pbool("zero", true));
  EXPECT_TRUE(P.Pbool("missing", true));
  EXPECT_THROW(P.Pbool("invalid", false), std::runtime_error);

  std::remove(filename);
}

TEST(Adapt, linint_requires_two_points) { // NOLINT
  EXPECT_THROW(LinInt(Vec{}), std::runtime_error);
  EXPECT_THROW(LinInt(Vec{{1.0, 2.0}}), std::runtime_error);
}

TEST(Adapt, int_with_to_throws_when_step_limit_is_exceeded) { // NOLINT
  const auto filename = "adapt_int_with_to.param";
  write_file(filename, "[param]\nLambda=2\nadapt=false\nxmax=2\nxfine=2\n");

  Params P(filename);
  Adapt calc(P, Sign::POS);
  calc.x = 0.0;
  calc.y = 0.0;

  EXPECT_THROW(calc.int_with_to(0.0, 1.0, []([[maybe_unused]] const auto x, [[maybe_unused]] const auto y) { return 0.0; }, false, 1e-10, 0), std::runtime_error);

  std::remove(filename);
}

TEST(Adapt, split_string_keeps_high_bit_bytes) { // NOLINT
  const std::string input = std::string(1, static_cast<char>(0x80)) + "1 2";
  const auto columns = split_string(input, 2);

  ASSERT_EQ(columns.size(), 2U);
  EXPECT_EQ(columns[0], std::string(1, static_cast<char>(0x80)) + "1");
  EXPECT_EQ(columns[1], "2");
}

TEST(Adapt, integral_method_matches_flat_band) { // NOLINT
  const auto filename = "adapt_integral_flat.param";
  write_file(filename,
             "[param]\n"
             "Lambda=10\n"
             "adapt=false\n"
             "f_method=integral\n"
             "xmax=30\n"
             "outputstep=0.25\n"
             "allowed_error=1e-10\n");

  Params P(filename);
  Adapt calc(P, Sign::POS, 0.01);
  calc.run();

  const auto values = load_g("FSOL.dat");
  const double lambda = 10.0;
  const double log_lambda = std::log(lambda);
  for (const auto &[x, f] : values) {
    const double expected = x <= 2.0
                              ? 2.0 - x + (1.0 - std::pow(lambda, 1.0 - x)) / log_lambda
                              : std::pow(lambda, 2.0 - x) * (1.0 - 1.0 / lambda) / log_lambda;
    const double expected_f = expected / std::pow(lambda, 2.0 - x);
    EXPECT_NEAR(f, expected_f, 1e-10) << "x=" << x;
  }

  std::remove("FSOL.dat");
  std::remove(filename);
}

TEST(Adapt, cumulative_inverse_uses_upper_plateau_edge) { // NOLINT
  const auto param_filename = "adapt_integral_plateau.param";
  const auto dos_filename = "adapt_integral_plateau.dat";
  write_file(param_filename,
             "[param]\n"
             "dos=adapt_integral_plateau.dat\n"
             "f_method=integral\n"
             "xmax=2\n");
  write_file(dos_filename,
             "0.1 0\n"
             "0.25 0\n"
             "0.4 1\n"
             "0.5 0\n"
             "0.7 0\n"
             "0.8 1\n"
             "1.0 1\n");

  Params P(param_filename);
  Adapt calc(P, Sign::POS);
  calc.load_init_rho();
  calc.init_cumulative();

  EXPECT_NEAR(calc.inverse_normalized_cumulative(0.0), 0.25, 1e-14);
  const auto plateau_weight = calc.normalized_cumulative(0.6);
  EXPECT_NEAR(calc.inverse_normalized_cumulative(plateau_weight), 0.7, 1e-14);
  EXPECT_LT(calc.inverse_normalized_cumulative(std::nextafter(plateau_weight, 0.0)), 0.5);
  EXPECT_GT(calc.inverse_normalized_cumulative(std::nextafter(plateau_weight, 1.0)), 0.7);

  std::remove(dos_filename);
  std::remove(param_filename);
}

TEST(Adapt, cumulative_inverse_extends_terminal_plateau_to_band_edge) { // NOLINT
  const auto param_filename = "adapt_integral_terminal_plateau.param";
  const auto dos_filename = "adapt_integral_terminal_plateau.dat";
  write_file(param_filename,
             "[param]\n"
             "dos=adapt_integral_terminal_plateau.dat\n"
             "f_method=integral\n"
             "xmax=2\n");
  write_file(dos_filename,
             "0.1 1\n"
             "0.5 1\n"
             "0.8 0\n"
             "0.9 0\n");

  Params P(param_filename);
  Adapt calc(P, Sign::POS);
  calc.load_init_rho();
  calc.init_cumulative();

  EXPECT_DOUBLE_EQ(calc.inverse_normalized_cumulative(1.0), 1.0);

  std::remove(dos_filename);
  std::remove(param_filename);
}

TEST(Adapt, integral_method_matches_linear_density_oracle) { // NOLINT
  const auto filename = "adapt_integral_linear.param";
  write_file(filename,
             "[param]\n"
             "Lambda=2\n"
             "adapt=false\n"
             "f_method=integral\n"
             "xmax=40\n"
             "allowed_error=1e-10\n");

  Params P(filename);
  Adapt calc(P, Sign::POS);
  calc.vecrho = {{0.0, 0.0}, {1.0, 2.0}};
  calc.rho = LinInt(calc.vecrho);
  auto integrated_rho = calc.vecrho;
  integrate(integrated_rho);
  calc.intrho1 = IntLinInt(calc.vecrho, integrated_rho);
  calc.intrho2 = calc.intrho1;
  calc.init_cumulative();
  calc.max_error = 0.0;

  constexpr std::size_t limit = 1000;
  std::unique_ptr<gsl_integration_cquad_workspace, GslWorkspaceDeleter> workspace(
    gsl_integration_cquad_workspace_alloc(limit));
  ASSERT_TRUE(workspace);
  const double log_lambda = std::log(2.0);
  for (const double x : {1.25, 2.0, 10.0, 20.0, 30.0, 40.0}) {
    const double expected_weight = x < 2.0
                                     ? 2.0 - x + (1.0 - std::pow(2.0, 2.0 - 2.0 * x)) / (2.0 * log_lambda)
                                     : std::pow(2.0, 4.0 - 2.0 * x) * (1.0 - 0.25) / (2.0 * log_lambda);
    const double expected = std::sqrt(expected_weight);
    EXPECT_NEAR(calc.Eps_integral(x, workspace.get()) / expected, 1.0, 1e-9) << "x=" << x;
  }

  std::remove(filename);
}

TEST(Adapt, rejects_unknown_f_method) { // NOLINT
  const auto filename = "adapt_invalid_f_method.param";
  write_file(filename, "[param]\nf_method=unknown\n");

  Params P(filename);
  EXPECT_THROW(Adapt(P, Sign::POS), std::invalid_argument);

  std::remove(filename);
}

TEST(Adapt, cquad_error_policies_handle_failed_results) { // NOLINT
  const auto filename = "adapt_cquad_policy.param";
  write_file(filename, "[param]\nf_method=integral\n");
  Params params(filename);

  Adapt default_fail(params, Sign::POS);
  EXPECT_THROW(default_fail.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0), std::runtime_error);

  CquadOptions ignore_options;
  ignore_options.gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;
  Adapt ignored(params, Sign::POS, std::nullopt, false, ignore_options);
  EXPECT_NO_THROW(ignored.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0));

  CquadOptions warn_options;
  warn_options.gsl_error_policy = NRG::Tools::GslErrorPolicy::warn;
  Adapt warned(params, Sign::POS, std::nullopt, false, warn_options);
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(warned.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0));
  const auto warning = testing::internal::GetCapturedStderr();
  EXPECT_NE(warning.find("adapt: warning: Integral method failed at x=2.000000"), std::string::npos);

  std::remove(filename);
}

TEST(Adapt, cquad_override_does_not_bypass_allowed_error_validation) { // NOLINT
  const auto filename = "adapt_invalid_allowed_error.param";
  write_file(filename, "[param]\nadapt=true\nf_method=integral\nallowed_error=-1\n");
  Params params(filename);
  CquadOptions options;
  options.epsrel = 1e-8;

  EXPECT_THROW(Adapt(params, Sign::POS, std::nullopt, false, options), std::invalid_argument);
  std::remove(filename);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
