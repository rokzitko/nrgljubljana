#ifndef NRG_TOOLS_COMMON_ADAPT_PARAMETERS_HPP
#define NRG_TOOLS_COMMON_ADAPT_PARAMETERS_HPP

#include <cfloat>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

#include <parse_bool.hpp>
#include "parser.hpp"
#include "tabulated_density.hpp"

namespace NRG::Tools {

struct AdaptParameters {
  double Lambda;
  bool adapt, hardgap;
  double boundary, bandrescale;
  InterpolationMethod density_interpolation;
  double xmax, xfine, output_step, dx_fine, dx_fast, allowed_error;
  int max_subdiv;
  double max_abs, convergence_eps, factor0;
  int max_iter;
};

// The parameter-only part of Adapt::set_parameters. Keep its defaults, parsing,
// checks and scalar arithmetic here; do not construct an Adapt or a GSL guard.
inline AdaptParameters resolve_adapt_parameters(const std::map<std::string, std::string> &parameters) {
  const auto number = [&](const std::string &key, double fallback) {
    return get_or_default(parameters, key, fallback, [&](const auto &value) { return parse_parameter_double(value, key); });
  };
  const auto integer = [&](const std::string &key, int fallback) {
    return get_or_default(parameters, key, fallback, [&](const auto &value) { return parse_parameter_int(value, key); });
  };
  const auto boolean = [&](const std::string &key, bool fallback) {
    return get_or_default(parameters, key, fallback, [&](const auto &value) {
      try { return NRG::parse_bool(value); }
      catch (const std::runtime_error &error) { throw std::runtime_error(key + ": " + error.what()); }
    });
  };

  AdaptParameters result;
  result.Lambda = number("Lambda", 2.0);
  if (!(result.Lambda > 1.0)) throw std::invalid_argument("Lambda must be greater than 1.");
  result.adapt = boolean("adapt", false);
  result.hardgap = boolean("hardgap", false);
  result.boundary = number("boundary", 0.0);
  if (result.hardgap && !(result.boundary >= 0.0 && result.boundary < 1.0))
    throw std::invalid_argument("boundary must be in [0,1) when hardgap=true.");
  result.bandrescale = number("bandrescale", 1.0);
  if (!(std::isfinite(result.bandrescale) && result.bandrescale > 0.0))
    throw std::invalid_argument("bandrescale must be positive and finite.");
  result.density_interpolation = parse_density_interpolation_method(
    get_or_default(parameters, "density_interpolation", std::string("linear"), [](const auto &value) { return value; }));
  result.xmax = number("xmax", 30);
  if (!(result.xmax > 1)) throw std::invalid_argument("xmax must be greater than 1.");
  result.xfine = number("xfine", 5);
  if (!(result.xfine > 1)) throw std::invalid_argument("xfine must be greater than 1.");
  result.output_step = number("outputstep", 1.0 / 64.0);
  if (!std::isfinite(result.output_step) || result.output_step <= 0.0 || result.output_step > 1.0)
    throw std::invalid_argument("outputstep must be finite and in the range (0, 1].");
  result.dx_fine = number("dx_fine", 1e-5);
  result.dx_fast = number("dx_fast", 1e-4);
  if (!(result.dx_fine > 0.0 && result.dx_fast > 0.0))
    throw std::invalid_argument("dx_fine and dx_fast must be positive and finite.");
  result.allowed_error = number("allowed_error", 1e-10);
  if (!(std::isfinite(result.allowed_error) && result.allowed_error > 0.0))
    throw std::invalid_argument("allowed_error must be positive and finite.");
  result.max_subdiv = integer("max_subdiv", 10);
  if (result.max_subdiv < 0) throw std::invalid_argument("max_subdiv must be nonnegative.");
  if (!(result.dx_fine * std::pow(0.5, result.max_subdiv) > DBL_EPSILON))
    throw std::invalid_argument("dx_fine and max_subdiv imply a sub-step below machine precision.");
  result.max_abs = number("max_abs", 100.0);
  if (!(result.max_abs > 0.0)) throw std::invalid_argument("max_abs must be greater than 0.");
  result.convergence_eps = number("secant_eps", 1e-4);
  if (!(result.convergence_eps > 0.0)) throw std::invalid_argument("secant_eps must be positive and finite.");
  const auto secant_factor = number("secant_factor", 1e-7);
  if (!(secant_factor > 0.0)) throw std::invalid_argument("secant_factor must be positive and finite.");
  result.factor0 = 1.0 + secant_factor;
  result.max_iter = integer("secant_max_iter", 10);
  if (result.max_iter <= 0) throw std::invalid_argument("secant_max_iter must be positive.");
  return result;
}

inline void validate_adapt_method(bool adapt, bool hardgap, double boundary, bool integral) {
  if (hardgap && boundary > 0.0 && (adapt || !integral))
    throw std::invalid_argument("Nonzero hardgap supports only a fixed mesh and the integral method; "
                                "use adapt=false and --integral (or f_method=integral).");
}

} // namespace NRG::Tools

#endif
