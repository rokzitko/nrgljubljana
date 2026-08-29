#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include <kk/kk.hpp>

namespace {

auto nonlinear_symmetric_input() {
  return NRG::KK::XYFUNC{{-5.0, 0.5}, {-3.0, -2.0}, {-1.25, 1.5}, {-0.4, 0.25},
                         {0.4, 3.0},  {1.25, -1.0}, {3.0, 2.5},   {5.0, -0.75}};
}

} // namespace

TEST(KK, empty_input_throws) { // NOLINT
  auto construct = [] { return NRG::KK::KK(NRG::KK::XYFUNC{}); };
  EXPECT_THROW(construct(), std::runtime_error);
}

TEST(KK, invalid_akima_input_throws) { // NOLINT
  auto construct = [] {
    return NRG::KK::KK(NRG::KK::XYFUNC{{-1.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}, {1.0, 0.0}});
  };
  EXPECT_THROW(construct(), std::runtime_error);
}

TEST(KK, interpolation_method_minimum_sizes) { // NOLINT
  using NRG::Tools::InterpolationMethod;

  EXPECT_EQ(NRG::Tools::interpolation_minimum_size(InterpolationMethod::linear), 2);
  EXPECT_EQ(NRG::Tools::interpolation_minimum_size(InterpolationMethod::cspline), 3);
  EXPECT_EQ(NRG::Tools::interpolation_minimum_size(InterpolationMethod::akima), 5);
  EXPECT_EQ(NRG::Tools::interpolation_minimum_size(InterpolationMethod::steffen), 3);

  NRG::KK::NumericalOptions options;
  options.interpolation = InterpolationMethod::linear;
  EXPECT_NO_THROW(NRG::KK::KK(NRG::KK::XYFUNC{{-1.0, 0.0}, {1.0, 1.0}}, options));

  options.interpolation = InterpolationMethod::cspline;
  EXPECT_NO_THROW(NRG::KK::KK(NRG::KK::XYFUNC{{-2.0, 0.0}, {-1.0, 1.0}, {1.0, 4.0}, {2.0, 2.0}}, options));
  try {
    NRG::KK::KK too_small(NRG::KK::XYFUNC{{-1.0, 0.0}, {1.0, 1.0}}, options);
    FAIL() << "cspline accepted fewer than three input points";
  } catch (const std::runtime_error &error) {
    EXPECT_NE(std::string(error.what()).find("cspline requires at least 3"), std::string::npos);
  }

  options.interpolation = InterpolationMethod::steffen;
  EXPECT_NO_THROW(NRG::KK::KK(NRG::KK::XYFUNC{{-2.0, 0.0}, {-1.0, 1.0}, {1.0, 4.0}, {2.0, 2.0}}, options));
  try {
    NRG::KK::KK too_small(NRG::KK::XYFUNC{{-1.0, 0.0}, {1.0, 1.0}}, options);
    FAIL() << "steffen accepted fewer than three input points";
  } catch (const std::runtime_error &error) {
    EXPECT_NE(std::string(error.what()).find("steffen requires at least 3"), std::string::npos);
  }

  options.interpolation = InterpolationMethod::akima;
  EXPECT_NO_THROW(NRG::KK::KK(NRG::KK::XYFUNC{{-3.0, 0.0}, {-1.5, 1.0}, {-0.5, 4.0},
                                               {0.5, 2.0},  {1.5, -1.0}, {3.0, 3.0}},
                                  options));
  try {
    NRG::KK::KK too_small(NRG::KK::XYFUNC{{-2.0, 0.0}, {-1.0, 1.0}, {1.0, 4.0}, {2.0, 2.0}}, options);
    FAIL() << "akima accepted fewer than five input points";
  } catch (const std::runtime_error &error) {
    EXPECT_NE(std::string(error.what()).find("akima requires at least 5"), std::string::npos);
  }
}

TEST(KK, all_interpolation_methods_handle_symmetric_nonlinear_input) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  constexpr std::array methods{InterpolationMethod::linear, InterpolationMethod::cspline, InterpolationMethod::akima,
                               InterpolationMethod::steffen};
  const NRG::KK::DVEC grid{-2.2, 0.7, 2.4};
  std::array<double, methods.size()> probe_values{};

  for (std::size_t method_index = 0; method_index < methods.size(); ++method_index) {
    const auto method = methods[method_index];
    SCOPED_TRACE(NRG::Tools::interpolation_method_name(method));
    NRG::KK::NumericalOptions options;
    options.interpolation = method;
    NRG::KK::KK transformer(nonlinear_symmetric_input(), options);
    const auto result = transformer.calc(grid);

    ASSERT_EQ(result.size(), grid.size());
    for (std::size_t point_index = 0; point_index < result.size(); ++point_index) {
      EXPECT_DOUBLE_EQ(result[point_index].first, grid[point_index]);
      EXPECT_TRUE(std::isfinite(result[point_index].second));
    }
    probe_values[method_index] = result[1].second;
  }

  EXPECT_GT(std::abs(probe_values[0] - probe_values[1]), 1e-4);
  EXPECT_GT(std::abs(probe_values[0] - probe_values[2]), 1e-4);
  EXPECT_GT(std::abs(probe_values[1] - probe_values[2]), 1e-4);
  EXPECT_GT(std::abs(probe_values[0] - probe_values[3]), 1e-4);
}

TEST(KK, numerical_options_defaults) { // NOLINT
  const NRG::KK::NumericalOptions options;
  EXPECT_EQ(options.interpolation, NRG::Tools::InterpolationMethod::akima);
  EXPECT_DOUBLE_EQ(options.epsabs, 1e-12);
  EXPECT_DOUBLE_EQ(options.epsrel, 1e-8);
  EXPECT_EQ(options.workspace_limit, 1000);
  EXPECT_EQ(options.quadrature_rule, NRG::Tools::QagRule::gauss15);
  EXPECT_EQ(options.gsl_error_policy, NRG::Tools::GslErrorPolicy::ignore);
}

TEST(KK, old_constructor_matches_explicit_defaults) { // NOLINT
  const auto input = nonlinear_symmetric_input();
  const NRG::KK::DVEC grid{-2.6, -0.7, 0.6, 2.2};
  NRG::KK::KK old_constructor(input);
  NRG::KK::KK explicit_defaults(input, NRG::KK::NumericalOptions{});

  const auto old_result = old_constructor.calc(grid);
  const auto default_result = explicit_defaults.calc(grid);
  ASSERT_EQ(old_result.size(), default_result.size());
  for (std::size_t index = 0; index < old_result.size(); ++index) {
    EXPECT_DOUBLE_EQ(old_result[index].first, default_result[index].first);
    EXPECT_NEAR(old_result[index].second, default_result[index].second, 1e-13);
  }
}

TEST(KK, numerical_options_validate_tolerances_and_workspace) { // NOLINT
  const auto input = nonlinear_symmetric_input();
  NRG::KK::NumericalOptions options;

  options.workspace_limit = 0;
  EXPECT_THROW(NRG::KK::KK(input, options), std::invalid_argument);

  options.workspace_limit = NRG::Tools::qag_workspace_limit_maximum() + 1;
  EXPECT_THROW(NRG::KK::KK(input, options), std::invalid_argument);

  options = {};
  options.epsabs = -1.0;
  EXPECT_THROW(NRG::KK::KK(input, options), std::invalid_argument);

  options = {};
  options.epsrel = -1.0;
  EXPECT_THROW(NRG::KK::KK(input, options), std::invalid_argument);

  options = {};
  options.epsabs = 0.0;
  options.epsrel = 0.0;
  EXPECT_THROW(NRG::KK::KK(input, options), std::invalid_argument);
}

TEST(KK, configured_rule_workspace_and_tolerances_are_used) { // NOLINT
  NRG::KK::NumericalOptions options;
  options.interpolation = NRG::Tools::InterpolationMethod::cspline;
  options.epsabs = 2e-9;
  options.epsrel = 3e-7;
  options.workspace_limit = 64;
  options.quadrature_rule = NRG::Tools::QagRule::gauss61;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;
  NRG::KK::KK transformer(nonlinear_symmetric_input(), options);

  const NRG::KK::DVEC grid{-1.7, 0.35, 2.1};
  const auto result = transformer.calc(grid);
  ASSERT_EQ(result.size(), grid.size());
  for (std::size_t index = 0; index < result.size(); ++index) {
    EXPECT_DOUBLE_EQ(result[index].first, grid[index]);
    EXPECT_TRUE(std::isfinite(result[index].second));
  }
  EXPECT_DOUBLE_EQ(transformer.calc(grid[1]), result[1].second);
}

TEST(KK, fail_policy_reports_deterministic_qag_failure) { // NOLINT
  NRG::KK::NumericalOptions options;
  options.epsabs = 1e-300;
  options.epsrel = 0.0;
  options.workspace_limit = 1;
  options.quadrature_rule = NRG::Tools::QagRule::gauss61;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::fail;
  NRG::KK::KK transformer(nonlinear_symmetric_input(), options);

  std::string message;
  try {
    (void)transformer.calc(0.73);
  } catch (const std::runtime_error &error) {
    message = error.what();
  }

  ASSERT_FALSE(message.empty()) << "the constrained QAG calculation unexpectedly converged";
  EXPECT_NE(message.find("GSL QAG failed for z="), std::string::npos);
  EXPECT_NE(message.find("epsabs=1e-300"), std::string::npos);
  EXPECT_NE(message.find("epsrel=0"), std::string::npos);
  EXPECT_NE(message.find("workspace_limit=1"), std::string::npos);
  EXPECT_NE(message.find("quadrature_rule=61"), std::string::npos);
}

TEST(KK, configured_qag_rule_changes_the_constrained_estimate) { // NOLINT
  NRG::KK::NumericalOptions options;
  options.epsabs = 1e-300;
  options.epsrel = 0.0;
  options.workspace_limit = 1;
  options.gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;

  options.quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::KK::KK gauss15(nonlinear_symmetric_input(), options);
  const auto result15 = gauss15.calc(0.73);

  options.quadrature_rule = NRG::Tools::QagRule::gauss61;
  NRG::KK::KK gauss61(nonlinear_symmetric_input(), options);
  const auto result61 = gauss61.calc(0.73);

  EXPECT_GT(std::abs(result15 - result61), 1e-8);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
