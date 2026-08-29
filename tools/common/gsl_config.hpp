#ifndef _tools_common_gsl_config_hpp_
#define _tools_common_gsl_config_hpp_

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace NRG::Tools {

class GslErrorHandlerGuard {
  gsl_error_handler_t *previous;

  public:
  GslErrorHandlerGuard() : previous{gsl_set_error_handler_off()} {}
  GslErrorHandlerGuard(const GslErrorHandlerGuard &) = delete;
  GslErrorHandlerGuard &operator=(const GslErrorHandlerGuard &) = delete;
  ~GslErrorHandlerGuard() {
    if (previous)
      gsl_set_error_handler(previous);
    else
      gsl_set_error_handler_off();
  }
};

enum class InterpolationMethod { linear, cspline, akima, steffen };

inline auto interpolation_method_name(const InterpolationMethod method) -> std::string_view {
  switch (method) {
    case InterpolationMethod::linear: return "linear";
    case InterpolationMethod::cspline: return "cspline";
    case InterpolationMethod::akima: return "akima";
    case InterpolationMethod::steffen: return "steffen";
  }
  throw std::logic_error("Unknown interpolation method.");
}

inline auto parse_interpolation_method(const std::string_view value) {
  if (value == "linear") return InterpolationMethod::linear;
  if (value == "cspline") return InterpolationMethod::cspline;
  if (value == "akima") return InterpolationMethod::akima;
  if (value == "steffen") return InterpolationMethod::steffen;
  throw std::invalid_argument("Interpolation method must be one of: linear, cspline, akima, steffen.");
}

inline auto gsl_interpolation_type(const InterpolationMethod method) -> const gsl_interp_type * {
  switch (method) {
    case InterpolationMethod::linear: return gsl_interp_linear;
    case InterpolationMethod::cspline: return gsl_interp_cspline;
    case InterpolationMethod::akima: return gsl_interp_akima;
    case InterpolationMethod::steffen: return gsl_interp_steffen;
  }
  throw std::logic_error("Unknown interpolation method.");
}

inline auto interpolation_minimum_size(const InterpolationMethod method) {
  return static_cast<std::size_t>(gsl_interp_type_min_size(gsl_interpolation_type(method)));
}

enum class QagRule { gauss15 = 15, gauss21 = 21, gauss31 = 31, gauss41 = 41, gauss51 = 51, gauss61 = 61 };

inline auto parse_qag_rule(const std::string_view value) {
  if (value == "15") return QagRule::gauss15;
  if (value == "21") return QagRule::gauss21;
  if (value == "31") return QagRule::gauss31;
  if (value == "41") return QagRule::gauss41;
  if (value == "51") return QagRule::gauss51;
  if (value == "61") return QagRule::gauss61;
  throw std::invalid_argument("Quadrature rule must be one of: 15, 21, 31, 41, 51, 61.");
}

inline auto gsl_qag_rule(const QagRule rule) {
  switch (rule) {
    case QagRule::gauss15: return GSL_INTEG_GAUSS15;
    case QagRule::gauss21: return GSL_INTEG_GAUSS21;
    case QagRule::gauss31: return GSL_INTEG_GAUSS31;
    case QagRule::gauss41: return GSL_INTEG_GAUSS41;
    case QagRule::gauss51: return GSL_INTEG_GAUSS51;
    case QagRule::gauss61: return GSL_INTEG_GAUSS61;
  }
  throw std::logic_error("Unknown quadrature rule.");
}

enum class GslErrorPolicy { ignore, warn, fail };

inline auto gsl_error_policy_name(const GslErrorPolicy policy) -> std::string_view {
  switch (policy) {
    case GslErrorPolicy::ignore: return "ignore";
    case GslErrorPolicy::warn: return "warn";
    case GslErrorPolicy::fail: return "fail";
  }
  throw std::logic_error("Unknown GSL error policy.");
}

inline auto parse_gsl_error_policy(const std::string_view value) {
  if (value == "ignore") return GslErrorPolicy::ignore;
  if (value == "warn") return GslErrorPolicy::warn;
  if (value == "fail") return GslErrorPolicy::fail;
  throw std::invalid_argument("GSL error policy must be one of: ignore, warn, fail.");
}

inline auto parse_finite_double(const std::string_view value, const std::string_view name) {
  const std::string text(value);
  char *end = nullptr;
  errno = 0;
  const auto result = std::strtod(text.c_str(), &end);
  const bool underflowed_to_zero = errno == ERANGE && result == 0.0;
  if (underflowed_to_zero || end == text.c_str() || end != text.c_str() + text.size() || !std::isfinite(result))
    throw std::invalid_argument(std::string(name) + " must be a finite representable number: " + text);
  return result;
}

inline auto parse_positive_size(const std::string_view value, const std::string_view name) {
  std::size_t result = 0;
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
  if (error != std::errc{} || end != value.data() + value.size() || result == 0)
    throw std::invalid_argument(std::string(name) + " must be a positive integer: " + std::string(value));
  return result;
}

inline constexpr auto qag_workspace_limit_maximum() {
  constexpr auto element_size = std::max(sizeof(double), sizeof(std::size_t));
  return std::numeric_limits<std::size_t>::max() / element_size;
}

inline constexpr auto cquad_workspace_limit_maximum() {
  constexpr auto element_size = std::max(sizeof(gsl_integration_cquad_ival), sizeof(std::size_t));
  return std::numeric_limits<std::size_t>::max() / element_size;
}

inline void validate_qag_workspace_limit(const std::size_t limit) {
  if (limit == 0) throw std::invalid_argument("QAG workspace limit must be positive.");
  if (limit > qag_workspace_limit_maximum()) throw std::invalid_argument("QAG workspace limit is too large.");
}

inline void validate_cquad_workspace_limit(const std::size_t limit) {
  if (limit < 3) throw std::invalid_argument("CQUAD workspace limit must be at least 3.");
  if (limit > cquad_workspace_limit_maximum()) throw std::invalid_argument("CQUAD workspace limit is too large.");
}

inline void validate_tolerances(const double epsabs, const double epsrel) {
  if (!std::isfinite(epsabs) || epsabs < 0.0) throw std::invalid_argument("Absolute integration tolerance must be finite and nonnegative.");
  if (!std::isfinite(epsrel) || epsrel < 0.0) throw std::invalid_argument("Relative integration tolerance must be finite and nonnegative.");
  const auto minimum_relative_tolerance = std::max(50.0 * std::numeric_limits<double>::epsilon(), 0.5e-28);
  if (epsabs == 0.0 && epsrel < minimum_relative_tolerance)
    throw std::invalid_argument("Relative integration tolerance is too small when absolute tolerance is zero.");
}

inline void validate_cquad_tolerances(const double epsabs, const double epsrel) {
  if (!std::isfinite(epsabs) || epsabs < 0.0) throw std::invalid_argument("Absolute integration tolerance must be finite and nonnegative.");
  if (!std::isfinite(epsrel) || epsrel < 0.0) throw std::invalid_argument("Relative integration tolerance must be finite and nonnegative.");
  if (epsabs == 0.0 && epsrel < std::numeric_limits<double>::epsilon())
    throw std::invalid_argument("Relative CQUAD tolerance is too small when absolute tolerance is zero.");
}

inline auto gsl_integration_failed(const int status, const double result, const double error) {
  return status != GSL_SUCCESS || !std::isfinite(result) || !std::isfinite(error);
}

} // namespace NRG::Tools

#endif
