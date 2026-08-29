#ifndef _tools_common_gsl_piecewise_polynomial_hpp_
#define _tools_common_gsl_piecewise_polynomial_hpp_

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "gsl_config.hpp"
#include "piecewise_polynomial.hpp"

namespace NRG::Tools {

struct GslSplineDeleter {
  void operator()(gsl_spline *spline) const {
    if (spline) gsl_spline_free(spline);
  }
};

inline auto make_gsl_piecewise_polynomial(const std::vector<double> &knots, const std::vector<double> &values,
                                          const InterpolationMethod method) {
  if (knots.size() != values.size()) throw std::invalid_argument("Interpolation grids must have equal sizes.");
  const auto minimum_size = interpolation_minimum_size(method);
  if (knots.size() < minimum_size)
    throw std::invalid_argument(std::string(interpolation_method_name(method)) + " interpolation requires at least "
                                + std::to_string(minimum_size) + " points.");
  for (std::size_t index = 0; index < knots.size(); ++index) {
    if (!std::isfinite(knots[index]) || !std::isfinite(values[index]))
      throw std::invalid_argument("Interpolation data must be finite.");
    if (index != 0 && !(knots[index - 1] < knots[index]))
      throw std::invalid_argument("Interpolation energies must be strictly increasing.");
  }
  if (!std::isfinite(knots.back() - knots.front()))
    throw std::invalid_argument("Interpolation domain width must be finite.");

  const GslErrorHandlerGuard error_handler;
  std::unique_ptr<gsl_spline, GslSplineDeleter> spline{gsl_spline_alloc(gsl_interpolation_type(method), knots.size())};
  if (!spline) throw std::runtime_error("Failed to allocate GSL interpolation spline.");
  if (const auto status = gsl_spline_init(spline.get(), knots.data(), values.data(), knots.size()); status != GSL_SUCCESS)
    throw std::runtime_error("Failed to initialize GSL interpolation spline: " + std::string(gsl_strerror(status)));

  std::vector<std::vector<double>> coefficients;
  coefficients.reserve(knots.size() - 1);
  for (std::size_t interval = 0; interval + 1 < knots.size(); ++interval) {
    const auto delta = values[interval + 1] - values[interval];
    if (method == InterpolationMethod::linear) {
      coefficients.push_back({values[interval], delta});
      continue;
    }

    double first_derivative = 0.0;
    double second_derivative = 0.0;
    if (const auto status = gsl_spline_eval_deriv_e(spline.get(), knots[interval], nullptr, &first_derivative); status != GSL_SUCCESS)
      throw std::runtime_error("Failed to evaluate GSL interpolation derivative: " + std::string(gsl_strerror(status)));
    if (const auto status = gsl_spline_eval_deriv2_e(spline.get(), knots[interval], nullptr, &second_derivative); status != GSL_SUCCESS)
      throw std::runtime_error("Failed to evaluate GSL interpolation second derivative: " + std::string(gsl_strerror(status)));
    const auto width = knots[interval + 1] - knots[interval];
    const auto linear = width * first_derivative;
    const auto quadratic = 0.5 * width * (width * second_derivative);
    coefficients.push_back({values[interval], linear, quadratic, delta - linear - quadratic});
  }
  return PiecewisePolynomial<double>{knots, std::move(coefficients)};
}

} // namespace NRG::Tools

#endif
