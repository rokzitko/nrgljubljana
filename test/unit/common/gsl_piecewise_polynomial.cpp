#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <common/gsl_piecewise_polynomial.hpp>

namespace {

struct AccelDeleter {
  void operator()(gsl_interp_accel *acceleration) const {
    if (acceleration) gsl_interp_accel_free(acceleration);
  }
};

struct SplineDeleter {
  void operator()(gsl_spline *spline) const {
    if (spline) gsl_spline_free(spline);
  }
};

} // namespace

TEST(GslPiecewisePolynomial, materialized_intervals_match_all_selected_gsl_interpolants) { // NOLINT
  using NRG::Tools::InterpolationMethod;
  const std::vector<double> knots{-3.0, -1.0, 0.5, 2.0, 5.0, 9.0};
  const std::vector<double> values{4.0, -2.0, 1.0, 8.0, -1.0, 6.0};
  constexpr std::array methods{InterpolationMethod::linear, InterpolationMethod::cspline, InterpolationMethod::akima,
                               InterpolationMethod::steffen};

  for (const auto method : methods) {
    SCOPED_TRACE(NRG::Tools::interpolation_method_name(method));
    const auto polynomial = NRG::Tools::make_gsl_piecewise_polynomial(knots, values, method);
    const NRG::Tools::GslErrorHandlerGuard error_handler;
    std::unique_ptr<gsl_interp_accel, AccelDeleter> acceleration{gsl_interp_accel_alloc()};
    std::unique_ptr<gsl_spline, SplineDeleter> spline{gsl_spline_alloc(NRG::Tools::gsl_interpolation_type(method), knots.size())};
    ASSERT_TRUE(acceleration);
    ASSERT_TRUE(spline);
    ASSERT_EQ(gsl_spline_init(spline.get(), knots.data(), values.data(), knots.size()), GSL_SUCCESS);

    for (std::size_t interval = 0; interval + 1 < knots.size(); ++interval) {
      for (const double fraction : {0.0, 0.1, 0.37, 0.8}) {
        const auto x = knots[interval] + fraction * (knots[interval + 1] - knots[interval]);
        EXPECT_NEAR(polynomial.evaluate(x), gsl_spline_eval(spline.get(), x, acceleration.get()), 2e-13);
      }
    }
    EXPECT_NEAR(polynomial.evaluate(knots.back()), values.back(), 1e-14);
    EXPECT_NEAR(polynomial.integral(), gsl_spline_eval_integ(spline.get(), knots.front(), knots.back(), acceleration.get()), 5e-13);
  }
}

TEST(GslPiecewisePolynomial, preserves_nonrounded_akima_corner_intervals) { // NOLINT
  const std::vector<double> knots{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<double> values{0.0, 0.0, 0.0, 1.0, 2.0, 3.0};
  const auto polynomial = NRG::Tools::make_gsl_piecewise_polynomial(knots, values, NRG::Tools::InterpolationMethod::akima);
  const NRG::Tools::GslErrorHandlerGuard error_handler;
  std::unique_ptr<gsl_interp_accel, AccelDeleter> acceleration{gsl_interp_accel_alloc()};
  std::unique_ptr<gsl_spline, SplineDeleter> spline{gsl_spline_alloc(gsl_interp_akima, knots.size())};
  ASSERT_TRUE(acceleration);
  ASSERT_TRUE(spline);
  ASSERT_EQ(gsl_spline_init(spline.get(), knots.data(), values.data(), knots.size()), GSL_SUCCESS);

  for (const double x : {1.75, 1.99, 2.0, 2.01, 2.25})
    EXPECT_NEAR(polynomial.evaluate(x), gsl_spline_eval(spline.get(), x, acceleration.get()), 2e-14);
}

TEST(GslPiecewisePolynomial, validates_grid_and_method_minimum) { // NOLINT
  EXPECT_THROW(NRG::Tools::make_gsl_piecewise_polynomial({0.0, 1.0}, {0.0}, NRG::Tools::InterpolationMethod::linear),
               std::invalid_argument);
  EXPECT_THROW(NRG::Tools::make_gsl_piecewise_polynomial({0.0, 1.0}, {0.0, 1.0}, NRG::Tools::InterpolationMethod::steffen),
               std::invalid_argument);
  EXPECT_THROW(NRG::Tools::make_gsl_piecewise_polynomial({0.0, 1.0, 1.0}, {0.0, 1.0, 2.0}, NRG::Tools::InterpolationMethod::cspline),
               std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
