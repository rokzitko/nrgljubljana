#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include <common/gsl_config.hpp>

namespace {

using NRG::Tools::GslErrorPolicy;
using NRG::Tools::InterpolationMethod;
using NRG::Tools::QagRule;

} // namespace

TEST(GslConfig, interpolation_methods_parse_map_and_report_minimum_sizes) { // NOLINT
  struct interpolation_case {
    std::string_view name;
    InterpolationMethod method;
    const gsl_interp_type *gsl_type;
    std::size_t minimum_size;
  };
  const std::array cases{
    interpolation_case{"linear", InterpolationMethod::linear, gsl_interp_linear, 2},
    interpolation_case{"cspline", InterpolationMethod::cspline, gsl_interp_cspline, 3},
    interpolation_case{"akima", InterpolationMethod::akima, gsl_interp_akima, 5},
  };

  for (const auto &[name, method, gsl_type, minimum_size] : cases) {
    EXPECT_EQ(NRG::Tools::parse_interpolation_method(name), method);
    EXPECT_EQ(NRG::Tools::interpolation_method_name(method), name);
    EXPECT_EQ(NRG::Tools::gsl_interpolation_type(method), gsl_type);
    EXPECT_EQ(NRG::Tools::interpolation_minimum_size(method), minimum_size);
  }
}

TEST(GslConfig, invalid_interpolation_method_is_rejected) { // NOLINT
  EXPECT_THROW(NRG::Tools::parse_interpolation_method("cubic"), std::invalid_argument);

  const auto invalid = static_cast<InterpolationMethod>(99);
  EXPECT_THROW(NRG::Tools::interpolation_method_name(invalid), std::logic_error);
  EXPECT_THROW(NRG::Tools::gsl_interpolation_type(invalid), std::logic_error);
  EXPECT_THROW(NRG::Tools::interpolation_minimum_size(invalid), std::logic_error);
}

TEST(GslConfig, qag_rules_parse_and_map_to_gsl) { // NOLINT
  struct qag_case {
    std::string_view name;
    QagRule rule;
    int gsl_rule;
  };
  constexpr std::array cases{
    qag_case{"15", QagRule::gauss15, GSL_INTEG_GAUSS15},
    qag_case{"21", QagRule::gauss21, GSL_INTEG_GAUSS21},
    qag_case{"31", QagRule::gauss31, GSL_INTEG_GAUSS31},
    qag_case{"41", QagRule::gauss41, GSL_INTEG_GAUSS41},
    qag_case{"51", QagRule::gauss51, GSL_INTEG_GAUSS51},
    qag_case{"61", QagRule::gauss61, GSL_INTEG_GAUSS61},
  };

  for (const auto &[name, rule, gsl_rule] : cases) {
    EXPECT_EQ(NRG::Tools::parse_qag_rule(name), rule);
    EXPECT_EQ(NRG::Tools::gsl_qag_rule(rule), gsl_rule);
  }

  EXPECT_THROW(NRG::Tools::parse_qag_rule("17"), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::gsl_qag_rule(static_cast<QagRule>(17)), std::logic_error);
}

TEST(GslConfig, gsl_error_policies_parse_and_report_names) { // NOLINT
  struct policy_case {
    std::string_view name;
    GslErrorPolicy policy;
  };
  constexpr std::array cases{
    policy_case{"ignore", GslErrorPolicy::ignore},
    policy_case{"warn", GslErrorPolicy::warn},
    policy_case{"fail", GslErrorPolicy::fail},
  };

  for (const auto &[name, policy] : cases) {
    EXPECT_EQ(NRG::Tools::parse_gsl_error_policy(name), policy);
    EXPECT_EQ(NRG::Tools::gsl_error_policy_name(policy), name);
  }

  EXPECT_THROW(NRG::Tools::parse_gsl_error_policy("error"), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::gsl_error_policy_name(static_cast<GslErrorPolicy>(99)), std::logic_error);
}

TEST(GslConfig, positive_size_parser_rejects_invalid_and_overflowing_values) { // NOLINT
  EXPECT_EQ(NRG::Tools::parse_positive_size("1", "Workspace limit"), std::size_t{1});
  EXPECT_EQ(NRG::Tools::parse_positive_size("4096", "Workspace limit"), std::size_t{4096});

  const auto maximum = std::numeric_limits<std::size_t>::max();
  const auto maximum_text = std::to_string(maximum);
  EXPECT_EQ(NRG::Tools::parse_positive_size(maximum_text, "Workspace limit"), maximum);

  for (const std::string_view invalid : {"", "0", "-1", "+1", "1.0", " 1", "1 ", "1x"})
    EXPECT_THROW(NRG::Tools::parse_positive_size(invalid, "Workspace limit"), std::invalid_argument) << invalid;

  const auto overflow = maximum_text + "0";
  EXPECT_THROW(NRG::Tools::parse_positive_size(overflow, "Workspace limit"), std::invalid_argument);
}

TEST(GslConfig, workspace_limits_prevent_invalid_or_overflowing_gsl_allocations) { // NOLINT
  EXPECT_NO_THROW(NRG::Tools::validate_qag_workspace_limit(1));
  EXPECT_NO_THROW(NRG::Tools::validate_qag_workspace_limit(NRG::Tools::qag_workspace_limit_maximum()));
  EXPECT_THROW(NRG::Tools::validate_qag_workspace_limit(0), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_qag_workspace_limit(NRG::Tools::qag_workspace_limit_maximum() + 1), std::invalid_argument);

  EXPECT_NO_THROW(NRG::Tools::validate_cquad_workspace_limit(3));
  EXPECT_NO_THROW(NRG::Tools::validate_cquad_workspace_limit(NRG::Tools::cquad_workspace_limit_maximum()));
  EXPECT_THROW(NRG::Tools::validate_cquad_workspace_limit(2), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_cquad_workspace_limit(NRG::Tools::cquad_workspace_limit_maximum() + 1), std::invalid_argument);
}

TEST(GslConfig, integration_tolerances_are_validated) { // NOLINT
  const auto minimum_relative = std::max(50.0 * std::numeric_limits<double>::epsilon(), 0.5e-28);
  EXPECT_NO_THROW(NRG::Tools::validate_tolerances(1e-12, 0.0));
  EXPECT_NO_THROW(NRG::Tools::validate_tolerances(0.0, minimum_relative));
  EXPECT_NO_THROW(NRG::Tools::validate_tolerances(0.0, 1e-8));

  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(NRG::Tools::validate_tolerances(-1.0, 1e-8), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(infinity, 1e-8), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(nan, 1e-8), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(1e-12, -1.0), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(1e-12, infinity), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(1e-12, nan), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_tolerances(0.0, minimum_relative / 2.0), std::invalid_argument);
}

TEST(GslConfig, cquad_uses_its_own_relative_tolerance_floor) { // NOLINT
  const auto epsilon = std::numeric_limits<double>::epsilon();
  EXPECT_NO_THROW(NRG::Tools::validate_cquad_tolerances(0.0, epsilon));
  EXPECT_NO_THROW(NRG::Tools::validate_cquad_tolerances(1e-12, 0.0));
  EXPECT_THROW(NRG::Tools::validate_cquad_tolerances(0.0, std::nextafter(epsilon, 0.0)), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_cquad_tolerances(-1.0, epsilon), std::invalid_argument);
  EXPECT_THROW(NRG::Tools::validate_cquad_tolerances(0.0, -1.0), std::invalid_argument);
}

TEST(GslConfig, integration_failure_classification_includes_status_and_nonfinite_outputs) { // NOLINT
  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(NRG::Tools::gsl_integration_failed(GSL_SUCCESS, 1.0, 0.0));
  EXPECT_TRUE(NRG::Tools::gsl_integration_failed(GSL_EFAILED, 1.0, 0.0));
  EXPECT_TRUE(NRG::Tools::gsl_integration_failed(GSL_SUCCESS, nan, 0.0));
  EXPECT_TRUE(NRG::Tools::gsl_integration_failed(GSL_SUCCESS, 1.0, infinity));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
