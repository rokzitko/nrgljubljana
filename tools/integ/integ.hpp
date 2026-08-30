#ifndef NRG_TOOLS_INTEG_INTEG_HPP
#define NRG_TOOLS_INTEG_INTEG_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "../common/diagnostics.hpp"
#include "../common/gsl_config.hpp"
#include "../common/gsl_piecewise_polynomial.hpp"
#include "../common/output_file.hpp"
#include "../common/tabulated.hpp"

namespace NRG::Integ {

using Polynomial = NRG::Tools::PiecewisePolynomial<double>;
using XYPoint = std::pair<double, double>;
using XYFunction = std::vector<XYPoint>;

inline constexpr auto OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;

enum class Quantity { total, bounded, positive, negative, absolute, negative_absolute, energy_moment, fermi };

inline auto quantity_name(const Quantity quantity) -> std::string_view {
  switch (quantity) {
    case Quantity::total: return "total";
    case Quantity::bounded: return "bounded";
    case Quantity::positive: return "positive";
    case Quantity::negative: return "negative";
    case Quantity::absolute: return "absolute";
    case Quantity::negative_absolute: return "negative-absolute";
    case Quantity::energy_moment: return "energy-moment";
    case Quantity::fermi: return "fermi";
  }
  throw std::logic_error("Unknown integration quantity.");
}

inline auto format_double(const double value) {
  std::ostringstream output;
  output << std::setprecision(OUTPUT_PRECISION) << value;
  return output.str();
}

struct Options {
  int verbosity = 0;
  Quantity quantity = Quantity::total;
  bool quantity_selected = false;
  std::optional<double> requested_lower;
  std::optional<double> requested_upper;
  double temperature = 1e-99;
  NRG::Tools::InterpolationMethod interpolation = NRG::Tools::InterpolationMethod::akima;
  double epsabs = 1e-12;
  double epsrel = 1e-8;
  std::size_t workspace_limit = 1000;
  NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::warn;
  std::vector<std::string> inputs;
  bool help = false;
};

inline void usage(std::ostream &output) {
  output << "Usage: integ [options] [input ...]\n"
            "\n"
            "With no input, or with input '-', read a strict two-column table from standard input.\n"
            "Exactly one quantity may be selected; the default is the total integral.\n"
            "\n"
            "Quantity selectors:\n"
            "      --total                         integral over the full input domain\n"
            "      --lower-bound A --upper-bound B integral over [A,B], clipped to the input domain\n"
            "  -p, --positive                      integral over positive x\n"
            "  -n, --negative                      integral over negative x\n"
            "  -a, --absolute                      integral of the absolute interpolant\n"
            "      --negative-absolute             absolute integral over negative x\n"
            "  -e, --energy-moment                 first moment: integral of x*f(x)\n"
            "  -f, --fermi                         Fermi-weighted integral\n"
            "\n"
            "Options:\n"
            "  -h, --help                          show this help\n"
            "  -v                                  increase diagnostic verbosity\n"
            "  -V, --version                       show project version\n"
            "  -i, --interpolation METHOD          linear, cspline, akima, or steffen\n"
            "  -t, -T TEMPERATURE                  nonnegative Fermi temperature\n"
            "  -w                                  ignore QAG errors\n"
            "      --epsabs VALUE                  absolute QAG tolerance\n"
            "      --epsrel VALUE                  relative QAG tolerance\n"
            "      --workspace-limit N             QAG workspace subdivision limit\n"
            "      --quadrature-rule RULE          15, 21, 31, 41, 51, or 61\n"
            "      --gsl-error-policy POLICY       ignore, warn, or fail\n";
}

inline void select_quantity(Options &options, const Quantity quantity) {
  if (options.quantity_selected)
    throw std::invalid_argument("Exactly one integration quantity may be requested per invocation.");
  options.quantity = quantity;
  options.quantity_selected = true;
}

inline auto parse_options(const int argc, char *const argv[]) {
  enum LongOption {
    OPT_EPSABS = 256,
    OPT_EPSREL,
    OPT_WORKSPACE_LIMIT,
    OPT_QUADRATURE_RULE,
    OPT_GSL_ERROR_POLICY,
    OPT_TOTAL,
    OPT_LOWER_BOUND,
    OPT_UPPER_BOUND,
    OPT_NEGATIVE_ABSOLUTE,
    OPT_ENERGY_MOMENT,
  };
  static const option long_options[] = {
    {"help", no_argument, nullptr, 'h'},
    {"interpolation", required_argument, nullptr, 'i'},
    {"positive", no_argument, nullptr, 'p'},
    {"negative", no_argument, nullptr, 'n'},
    {"absolute", no_argument, nullptr, 'a'},
    {"fermi", no_argument, nullptr, 'f'},
    {"energy-moment", no_argument, nullptr, OPT_ENERGY_MOMENT},
    {"negative-absolute", no_argument, nullptr, OPT_NEGATIVE_ABSOLUTE},
    {"total", no_argument, nullptr, OPT_TOTAL},
    {"lower-bound", required_argument, nullptr, OPT_LOWER_BOUND},
    {"upper-bound", required_argument, nullptr, OPT_UPPER_BOUND},
    {"epsabs", required_argument, nullptr, OPT_EPSABS},
    {"epsrel", required_argument, nullptr, OPT_EPSREL},
    {"workspace-limit", required_argument, nullptr, OPT_WORKSPACE_LIMIT},
    {"quadrature-rule", required_argument, nullptr, OPT_QUADRATURE_RULE},
    {"gsl-error-policy", required_argument, nullptr, OPT_GSL_ERROR_POLICY},
    {nullptr, 0, nullptr, 0},
  };

  Options options;
  opterr = 0;
  optind = 1;
  while (true) {
    const auto selected = getopt_long(argc, argv, "aehi:npfvwt:T:", long_options, nullptr);
    if (selected == -1) break;
    switch (selected) {
      case 'h':
        options.help = true;
        return options;
      case 'v': ++options.verbosity; break;
      case 'w': options.gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore; break;
      case 't':
      case 'T': options.temperature = NRG::Tools::parse_finite_double(optarg, "Temperature"); break;
      case 'i': options.interpolation = NRG::Tools::parse_interpolation_method(optarg); break;
      case 'p': select_quantity(options, Quantity::positive); break;
      case 'n': select_quantity(options, Quantity::negative); break;
      case 'a': select_quantity(options, Quantity::absolute); break;
      case 'e': select_quantity(options, Quantity::energy_moment); break;
      case 'f': select_quantity(options, Quantity::fermi); break;
      case OPT_TOTAL: select_quantity(options, Quantity::total); break;
      case OPT_NEGATIVE_ABSOLUTE: select_quantity(options, Quantity::negative_absolute); break;
      case OPT_ENERGY_MOMENT: select_quantity(options, Quantity::energy_moment); break;
      case OPT_LOWER_BOUND:
        if (options.requested_lower)
          throw std::invalid_argument("The lower integration bound may be specified only once.");
        options.requested_lower = NRG::Tools::parse_finite_double(optarg, "Lower integration bound");
        break;
      case OPT_UPPER_BOUND:
        if (options.requested_upper)
          throw std::invalid_argument("The upper integration bound may be specified only once.");
        options.requested_upper = NRG::Tools::parse_finite_double(optarg, "Upper integration bound");
        break;
      case OPT_EPSABS:
        options.epsabs = NRG::Tools::parse_finite_double(optarg, "Absolute integration tolerance");
        break;
      case OPT_EPSREL:
        options.epsrel = NRG::Tools::parse_finite_double(optarg, "Relative integration tolerance");
        break;
      case OPT_WORKSPACE_LIMIT:
        options.workspace_limit = NRG::Tools::parse_positive_size(optarg, "Workspace limit");
        break;
      case OPT_QUADRATURE_RULE: options.quadrature_rule = NRG::Tools::parse_qag_rule(optarg); break;
      case OPT_GSL_ERROR_POLICY:
        options.gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg);
        break;
      case '?':
      default: throw std::invalid_argument("Unknown or incomplete command-line option.");
    }
  }

  if (options.temperature < 0.0) throw std::invalid_argument("Temperature must be nonnegative.");
  if (options.requested_lower || options.requested_upper) {
    if (!options.requested_lower || !options.requested_upper)
      throw std::invalid_argument("Bounded integration requires both --lower-bound and --upper-bound.");
    if (options.quantity_selected)
      throw std::invalid_argument("Bounds cannot be combined with another integration quantity.");
    if (*options.requested_lower > *options.requested_upper)
      throw std::invalid_argument("The lower integration bound must not exceed the upper bound.");
    options.quantity = Quantity::bounded;
    options.quantity_selected = true;
  }

  NRG::Tools::validate_tolerances(options.epsabs, options.epsrel);
  NRG::Tools::validate_qag_workspace_limit(options.workspace_limit);
  for (auto index = optind; index < argc; ++index) options.inputs.emplace_back(argv[index]);
  return options;
}

inline auto read_input(const std::vector<std::string> &sources) {
  XYFunction data;
  bool standard_input_read = false;
  const auto append = [&](XYFunction records) {
    data.insert(data.end(), std::make_move_iterator(records.begin()), std::make_move_iterator(records.end()));
  };

  if (sources.empty()) {
    append(NRG::Tools::read_strict_pairs(std::cin, "<stdin>"));
  } else {
    for (const auto &source : sources) {
      if (source == "-") {
        if (standard_input_read) throw std::invalid_argument("Standard input may be specified only once.");
        standard_input_read = true;
        append(NRG::Tools::read_strict_pairs(std::cin, "<stdin>"));
        continue;
      }
      std::ifstream input(source);
      if (!input) throw std::runtime_error("Can't open " + source + " for reading.");
      append(NRG::Tools::read_strict_pairs(input, source));
    }
  }
  if (data.empty()) throw std::runtime_error("No input data points provided.");
  std::sort(data.begin(), data.end());
  return data;
}

inline auto make_polynomial(const XYFunction &data, const NRG::Tools::InterpolationMethod method) {
  std::vector<double> knots;
  std::vector<double> values;
  knots.reserve(data.size());
  values.reserve(data.size());
  for (const auto &[knot, value] : data) {
    knots.push_back(knot);
    values.push_back(value);
  }
  return NRG::Tools::make_gsl_piecewise_polynomial(knots, values, method);
}

inline auto clipped_bounds(const Polynomial &polynomial, const double lower, const double upper)
  -> std::optional<std::pair<double, double>> {
  const auto clipped_lower = std::max(lower, polynomial.lower_bound());
  const auto clipped_upper = std::min(upper, polynomial.upper_bound());
  if (!(clipped_lower < clipped_upper)) return std::nullopt;
  return std::pair{clipped_lower, clipped_upper};
}

inline auto integrate_signed_clipped(const Polynomial &polynomial, const double lower, const double upper) {
  const auto bounds = clipped_bounds(polynomial, lower, upper);
  return bounds ? polynomial.integral(bounds->first, bounds->second) : 0.0;
}

inline auto integrate_absolute_clipped(const Polynomial &polynomial, const double lower, const double upper) {
  const auto bounds = clipped_bounds(polynomial, lower, upper);
  return bounds ? NRG::Tools::absolute_integral(polynomial, bounds->first, bounds->second) : 0.0;
}

inline auto fermi_factor(const double energy, const double temperature) {
  if (temperature == 0.0) return energy < 0.0 ? 1.0 : energy > 0.0 ? 0.0 : 0.5;
  const auto argument = energy / temperature;
  if (argument >= 0.0) {
    const auto exponential = std::exp(-argument);
    return exponential / (1.0 + exponential);
  }
  return 1.0 / (1.0 + std::exp(argument));
}

inline auto fermi_weighted_value(const double value, const double energy, const double temperature) {
  if (temperature == 0.0) return value * fermi_factor(energy, temperature);
  const auto argument = static_cast<long double>(energy) / static_cast<long double>(temperature);
  if (argument < 0.0L)
    return static_cast<double>(static_cast<long double>(value) / (1.0L + std::exp(argument)));
  if (value == 0.0) return value;
  if (!std::isfinite(value)) return value * fermi_factor(energy, temperature);
  const auto logarithm = std::log(std::abs(static_cast<long double>(value))) - argument
                         - std::log1p(std::exp(-argument));
  return std::copysign(static_cast<double>(std::exp(logarithm)), value);
}

inline auto energy_moment(const Polynomial &polynomial) {
  if constexpr (NRG::Tools::detail::native_extended_precision) {
    NRG::Tools::CompensatedLongComplexSum sum;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left = static_cast<long double>(polynomial.knots()[interval]);
      const auto width = static_cast<long double>(polynomial.knots()[interval + 1]) - left;
      NRG::Tools::CompensatedLongComplexSum local_sum;
      for (std::size_t power = 0; power < polynomial.coefficients()[interval].size(); ++power) {
        const auto coefficient = static_cast<long double>(polynomial.coefficients()[interval][power]);
        const auto integrated_weight = left / static_cast<long double>(power + 1)
                                       + width / static_cast<long double>(power + 2);
        local_sum.add({coefficient * (width * integrated_weight), 0.0L});
      }
      sum.add(local_sum.value());
    }
    return static_cast<double>(sum.value().real());
  } else {
    NRG::Tools::PiecewiseWideFloat sum{};
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const NRG::Tools::PiecewiseWideFloat left{polynomial.knots()[interval]};
      const auto width = NRG::Tools::PiecewiseWideFloat{polynomial.knots()[interval + 1]} - left;
      for (std::size_t power = 0; power < polynomial.coefficients()[interval].size(); ++power) {
        const NRG::Tools::PiecewiseWideFloat coefficient{polynomial.coefficients()[interval][power]};
        const auto integrated_weight = left / NRG::Tools::PiecewiseWideFloat{power + 1}
                                       + width / NRG::Tools::PiecewiseWideFloat{power + 2};
        sum += coefficient * width * integrated_weight;
      }
    }
    return static_cast<double>(sum);
  }
}

struct GslWorkspaceDeleter {
  void operator()(gsl_integration_workspace *workspace) const {
    if (workspace) gsl_integration_workspace_free(workspace);
  }
};

struct LocalFermiIntegrand {
  const std::vector<double> *coefficients;
  double coefficient_scale;
  double interval_left;
  double interval_right;
  double segment_left;
  double segment_right;
  double temperature;

  auto operator()(const double unit_coordinate) const {
    if constexpr (NRG::Tools::detail::native_extended_precision) {
      const auto normalized_left = NRG::Tools::detail::normalized_interval_coordinate(
        segment_left, interval_left, interval_right);
      const auto normalized_right = NRG::Tools::detail::normalized_interval_coordinate(
        segment_right, interval_left, interval_right);
      const auto coordinate = static_cast<long double>(unit_coordinate);
      const auto normalized = normalized_left + coordinate * (normalized_right - normalized_left);
      auto value = 0.0L;
      for (auto coefficient = coefficients->rbegin(); coefficient != coefficients->rend(); ++coefficient)
        value = value * normalized
                + static_cast<long double>(*coefficient) / static_cast<long double>(coefficient_scale);

      const auto scaled_left = static_cast<long double>(segment_left) / static_cast<long double>(temperature);
      const auto scaled_right = static_cast<long double>(segment_right) / static_cast<long double>(temperature);
      const auto argument = scaled_left + coordinate * (scaled_right - scaled_left);
      const auto log_fermi = [](const long double scaled_argument) {
        return scaled_argument >= 0.0L
                 ? -scaled_argument - std::log1p(std::exp(-scaled_argument))
                 : -std::log1p(std::exp(scaled_argument));
      };
      const auto factor_ratio = std::exp(log_fermi(argument) - log_fermi(scaled_left));
      return static_cast<double>(value * factor_ratio);
    } else {
      const auto normalized_left = NRG::Tools::detail::wide_normalized_interval_coordinate(
        segment_left, interval_left, interval_right);
      const auto normalized_right = NRG::Tools::detail::wide_normalized_interval_coordinate(
        segment_right, interval_left, interval_right);
      const NRG::Tools::PiecewiseWideFloat coordinate{unit_coordinate};
      const auto normalized = normalized_left + coordinate * (normalized_right - normalized_left);
      NRG::Tools::PiecewiseWideFloat value{};
      for (auto coefficient = coefficients->rbegin(); coefficient != coefficients->rend(); ++coefficient)
        value = value * normalized
                + NRG::Tools::PiecewiseWideFloat{*coefficient}
                    / NRG::Tools::PiecewiseWideFloat{coefficient_scale};

      const auto scaled_left = NRG::Tools::PiecewiseWideFloat{segment_left}
                               / NRG::Tools::PiecewiseWideFloat{temperature};
      const auto scaled_right = NRG::Tools::PiecewiseWideFloat{segment_right}
                                / NRG::Tools::PiecewiseWideFloat{temperature};
      const auto argument = scaled_left + coordinate * (scaled_right - scaled_left);
      const auto log_fermi = [](const NRG::Tools::PiecewiseWideFloat &input) {
        return input >= 0 ? -input - boost::multiprecision::log1p(boost::multiprecision::exp(-input))
                          : -boost::multiprecision::log1p(boost::multiprecision::exp(input));
      };
      const auto factor_ratio = boost::multiprecision::exp(log_fermi(argument) - log_fermi(scaled_left));
      return static_cast<double>(value * factor_ratio);
    }
  }
};

template <typename Function>
inline double invoke_gsl_function(const double argument, void *function) {
  return (*static_cast<Function *>(function))(argument);
}

inline auto integrate_fermi(const Polynomial &polynomial, const Options &options) {
  if (options.temperature == 0.0)
    return integrate_signed_clipped(polynomial, polynomial.lower_bound(), 0.0);

  const NRG::Tools::GslErrorHandlerGuard error_handler;
  std::unique_ptr<gsl_integration_workspace, GslWorkspaceDeleter> workspace{
    gsl_integration_workspace_alloc(options.workspace_limit)};
  if (!workspace) throw std::runtime_error("Failed to allocate GSL integration workspace.");

  NRG::Tools::CompensatedLongComplexSum extended_result_sum;
  NRG::Tools::CompensatedLongComplexSum extended_error_sum;
  NRG::Tools::PiecewiseWideFloat wide_result_sum{};
  NRG::Tools::PiecewiseWideFloat wide_error_sum{};
  std::size_t failed_calls = 0;
  int first_status = GSL_SUCCESS;
  auto first_lower = 0.0;
  auto first_upper = 0.0;

  // Fixed-width x/T segments resolve the Fermi layer. Beyond 4096, even the
  // largest finite binary64 cubic and domain contribute less than denorm_min.
  constexpr auto thermal_step = 16;
  constexpr auto thermal_cutoff = 4096;
  std::vector<double> thermal_boundaries{0.0};
  const NRG::Tools::PiecewiseWideFloat wide_temperature{options.temperature};
  const NRG::Tools::PiecewiseWideFloat maximum_double{std::numeric_limits<double>::max()};
  for (auto multiple = thermal_step; multiple <= thermal_cutoff; multiple += thermal_step) {
    const auto scale = wide_temperature * multiple;
    if (scale > maximum_double) break;
    const auto boundary = static_cast<double>(scale);
    if (boundary == 0.0) continue;
    thermal_boundaries.push_back(-boundary);
    thermal_boundaries.push_back(boundary);
  }
  std::sort(thermal_boundaries.begin(), thermal_boundaries.end());
  thermal_boundaries.erase(std::unique(thermal_boundaries.begin(), thermal_boundaries.end()),
                           thermal_boundaries.end());

  for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
    const auto interval_left = polynomial.knots()[interval];
    const auto interval_right = polynomial.knots()[interval + 1];
    std::vector<double> boundaries{interval_left, interval_right};
    const auto first_thermal = std::upper_bound(thermal_boundaries.begin(), thermal_boundaries.end(), interval_left);
    const auto last_thermal = std::lower_bound(first_thermal, thermal_boundaries.end(), interval_right);
    boundaries.insert(boundaries.end(), first_thermal, last_thermal);
    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());

    auto coefficient_scale = 0.0;
    for (const auto coefficient : polynomial.coefficients()[interval])
      coefficient_scale = std::max(coefficient_scale, std::abs(coefficient));
    if (coefficient_scale == 0.0) continue;

    for (std::size_t segment = 0; segment + 1 < boundaries.size(); ++segment) {
      const auto lower = boundaries[segment];
      const auto upper = boundaries[segment + 1];
      const auto wide_scaled_lower = NRG::Tools::PiecewiseWideFloat{lower} / wide_temperature;
      if (lower > 0.0 && wide_scaled_lower >= thermal_cutoff) continue;

      // Keep QAG's callback near unit scale, then restore amplitude, width,
      // and the leading Fermi factor in extended or multiprecision arithmetic.
      long double extended_scale = 0.0L;
      NRG::Tools::PiecewiseWideFloat wide_scale{};
      double local_epsabs = 0.0;
      if constexpr (NRG::Tools::detail::native_extended_precision) {
        const auto scaled_lower = static_cast<long double>(lower) / static_cast<long double>(options.temperature);
        const auto log_fermi = scaled_lower >= 0.0L
                                 ? -scaled_lower - std::log1p(std::exp(-scaled_lower))
                                 : -std::log1p(std::exp(scaled_lower));
        const auto span = static_cast<long double>(upper) - static_cast<long double>(lower);
        extended_scale = static_cast<long double>(coefficient_scale) * span * std::exp(log_fermi);
        if (extended_scale == 0.0L) continue;
        if (options.epsabs != 0.0) {
          const auto scaled_tolerance = static_cast<long double>(options.epsabs) / std::abs(extended_scale);
          local_epsabs = static_cast<double>(std::clamp(
            scaled_tolerance, static_cast<long double>(std::numeric_limits<double>::denorm_min()),
            static_cast<long double>(std::numeric_limits<double>::max())));
        }
      } else {
        const auto log_fermi = wide_scaled_lower >= 0
                                 ? -wide_scaled_lower
                                     - boost::multiprecision::log1p(boost::multiprecision::exp(-wide_scaled_lower))
                                 : -boost::multiprecision::log1p(boost::multiprecision::exp(wide_scaled_lower));
        const auto span = NRG::Tools::PiecewiseWideFloat{upper} - NRG::Tools::PiecewiseWideFloat{lower};
        wide_scale = NRG::Tools::PiecewiseWideFloat{coefficient_scale} * span
                     * boost::multiprecision::exp(log_fermi);
        if (wide_scale == 0) continue;
        if (options.epsabs != 0.0) {
          const auto scaled_tolerance = NRG::Tools::PiecewiseWideFloat{options.epsabs}
                                        / boost::multiprecision::abs(wide_scale);
          local_epsabs = static_cast<double>(std::max(
            NRG::Tools::PiecewiseWideFloat{std::numeric_limits<double>::denorm_min()},
            std::min(scaled_tolerance, maximum_double)));
        }
      }

      LocalFermiIntegrand integrand{&polynomial.coefficients()[interval], coefficient_scale, interval_left,
                                    interval_right, lower, upper, options.temperature};
      gsl_function function{&invoke_gsl_function<LocalFermiIntegrand>, &integrand};
      auto normalized_result = std::numeric_limits<double>::quiet_NaN();
      auto normalized_error = std::numeric_limits<double>::quiet_NaN();
      const auto status = gsl_integration_qag(&function, 0.0, 1.0, local_epsabs, options.epsrel,
                                              options.workspace_limit,
                                              NRG::Tools::gsl_qag_rule(options.quadrature_rule), workspace.get(),
                                              &normalized_result, &normalized_error);
      const auto failed = NRG::Tools::gsl_integration_failed(status, normalized_result, normalized_error);
      if (failed) {
        if (failed_calls == 0) {
          first_status = status;
          first_lower = lower;
          first_upper = upper;
        }
        ++failed_calls;
      }
      if constexpr (NRG::Tools::detail::native_extended_precision) {
        const auto scaled_result = static_cast<long double>(normalized_result) * extended_scale;
        const auto scaled_error = std::abs(static_cast<long double>(normalized_error) * extended_scale);
        extended_result_sum.add({scaled_result, 0.0L});
        extended_error_sum.add({scaled_error, 0.0L});
        if (options.verbosity >= 2)
          std::cerr << std::setprecision(OUTPUT_PRECISION) << "integ: qag interval=[" << lower << ',' << upper
                    << "] result=" << scaled_result << " estimated_error=" << scaled_error << '\n';
      } else {
        const auto scaled_result = NRG::Tools::PiecewiseWideFloat{normalized_result} * wide_scale;
        const auto scaled_error = boost::multiprecision::abs(
          NRG::Tools::PiecewiseWideFloat{normalized_error} * wide_scale);
        wide_result_sum += scaled_result;
        wide_error_sum += scaled_error;
        if (options.verbosity >= 2)
          std::cerr << std::setprecision(OUTPUT_PRECISION) << "integ: qag interval=[" << lower << ',' << upper
                    << "] result=" << scaled_result << " estimated_error=" << scaled_error << '\n';
      }
    }
  }
  const auto result = NRG::Tools::detail::native_extended_precision
                        ? static_cast<double>(extended_result_sum.value().real())
                        : static_cast<double>(wide_result_sum);
  const auto estimated_error = NRG::Tools::detail::native_extended_precision
                                 ? static_cast<double>(std::abs(extended_error_sum.value().real()))
                                 : static_cast<double>(boost::multiprecision::abs(wide_error_sum));
  const auto requested_error = std::max(options.epsabs, options.epsrel * std::abs(result));
  const auto global_failure = !std::isfinite(result) || !std::isfinite(estimated_error)
                              || estimated_error > requested_error;
  if (global_failure && options.gsl_error_policy != NRG::Tools::GslErrorPolicy::ignore) {
    auto message = "qag error: summed error estimate " + format_double(estimated_error)
                   + " exceeds the requested global tolerance " + format_double(requested_error);
    if (failed_calls != 0)
      message += "; " + std::to_string(failed_calls) + " local call(s) failed; first interval ["
                 + format_double(first_lower) + ", " + format_double(first_upper) + "]: "
                 + std::to_string(first_status) + " -- " + gsl_strerror(first_status);
    if (options.gsl_error_policy == NRG::Tools::GslErrorPolicy::fail) throw std::runtime_error(message);
    std::cerr << "integ: warning: " << message << '\n';
  }
  if (options.verbosity >= 1)
    std::cerr << std::setprecision(OUTPUT_PRECISION) << "integ: qag estimated_error=" << estimated_error << '\n';
  return result;
}

inline auto calculate(const Polynomial &polynomial, const Options &options) {
  double result = 0.0;
  switch (options.quantity) {
    case Quantity::total: result = polynomial.integral(); break;
    case Quantity::bounded:
      result = integrate_signed_clipped(polynomial, *options.requested_lower, *options.requested_upper);
      break;
    case Quantity::positive:
      result = integrate_signed_clipped(polynomial, 0.0, polynomial.upper_bound());
      break;
    case Quantity::negative:
      result = integrate_signed_clipped(polynomial, polynomial.lower_bound(), 0.0);
      break;
    case Quantity::absolute: result = NRG::Tools::absolute_integral(polynomial); break;
    case Quantity::negative_absolute:
      result = integrate_absolute_clipped(polynomial, polynomial.lower_bound(), 0.0);
      break;
    case Quantity::energy_moment: result = energy_moment(polynomial); break;
    case Quantity::fermi: result = integrate_fermi(polynomial, options); break;
  }
  if (!std::isfinite(result)) throw std::runtime_error("Integration produced a nonfinite result.");
  return result;
}

inline auto input_description(const std::vector<std::string> &inputs) {
  if (inputs.empty()) return std::string{"<stdin>"};
  std::string result;
  for (const auto &input : inputs) {
    if (!result.empty()) result += ',';
    result += input == "-" ? "<stdin>" : input;
  }
  return result;
}

inline void report_configuration(const Options &options, const Polynomial &polynomial, const std::size_t points) {
  if (options.verbosity == 0) return;
  NRG::Tools::ConfigurationReport report("integ");
  report.value("verbosity", options.verbosity);
  report.value("input", input_description(options.inputs));
  report.value("input.points", points);
  report.resolved("input.lower_bound", polynomial.lower_bound(), "sorted input data");
  report.resolved("input.upper_bound", polynomial.upper_bound(), "sorted input data");
  report.value("output.quantity", quantity_name(options.quantity));
  report.value("output.precision", OUTPUT_PRECISION);
  report.value("interpolation", NRG::Tools::interpolation_method_name(options.interpolation));
  report.value("algorithm", options.quantity == Quantity::fermi && options.temperature > 0.0
                              ? "adaptive-qag-per-polynomial-interval"
                              : "analytic-piecewise-polynomial");
  if (options.quantity == Quantity::bounded) {
    report.value("integration.requested_lower_bound", *options.requested_lower);
    report.value("integration.requested_upper_bound", *options.requested_upper);
  }
  report.value("temperature", options.temperature);
  report.value("integration.epsabs", options.epsabs);
  report.value("integration.epsrel", options.epsrel);
  report.value("integration.workspace_limit", options.workspace_limit);
  report.value("integration.quadrature_rule", static_cast<int>(options.quadrature_rule));
  report.value("integration.gsl_error_policy", NRG::Tools::gsl_error_policy_name(options.gsl_error_policy));
  report.write(std::cerr);
}

inline void run(const int argc, char *const argv[]) {
  const auto options = parse_options(argc, argv);
  if (options.help) {
    usage(std::cout);
    NRG::Tools::finish_output(std::cout, "<stdout>");
    return;
  }
  const auto data = read_input(options.inputs);
  const auto polynomial = make_polynomial(data, options.interpolation);
  report_configuration(options, polynomial, data.size());
  const auto result = calculate(polynomial, options);
  std::cout << std::setprecision(OUTPUT_PRECISION) << result << '\n';
  NRG::Tools::finish_output(std::cout, "<stdout>");
}

} // namespace NRG::Integ

#endif
