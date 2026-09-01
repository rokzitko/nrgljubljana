// Computes \int E^n rho(E)/(z-E) dE for given z and tabulated density of states rho(E).
//
// Legacy modes:
// Mode 1: <Re z> <Im z> as input. Returns Im part by default, or both Re and Im parts if the -G switch is used.
// Mode 2: read x and y from a file.
// Mode 3: convert imsigma/resigma.dat to imaw/reaw.dat files in the DMFT loop.
//
// Rok Zitko, rok.zitko@ijs.si

#ifndef _hilb_hilb_hpp_
#define _hilb_hilb_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <exception>
#include <ios>
#include <istream>
#include <ostream>
#include <array>
#include <complex>
#include <utility>
#include <functional>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <optional>
#include <memory>
#include <mutex>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <cassert>
#include <charconv>
#include <cerrno>
#include <cfloat>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <thread>
#include <type_traits>
#include <getopt.h>
#include <unistd.h>

#include <gsl/gsl_errno.h> // GNU scientific library
#include <gsl/gsl_integration.h>

#include "../common/gsl_config.hpp"
#include "../common/diagnostics.hpp"
#include "../common/gsl_piecewise_polynomial.hpp"
#include "../common/output_file.hpp"
#include "../common/parallel.hpp"

namespace NRG::Hilb {

using std::size_t;

inline constexpr int OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;

enum class Algorithm { analytic, qag };

inline auto algorithm_name(const Algorithm algorithm) -> std::string_view {
  switch (algorithm) {
    case Algorithm::analytic: return "analytic";
    case Algorithm::qag: return "qag";
  }
  throw std::logic_error("Unknown Hilbert-transform algorithm.");
}

inline auto parse_algorithm(const std::string_view value) {
  if (value == "analytic") return Algorithm::analytic;
  if (value == "qag") return Algorithm::qag;
  throw std::invalid_argument("Algorithm must be one of: analytic, qag.");
}

inline auto format_double(const double value) {
  std::ostringstream output;
  output << std::setprecision(OUTPUT_PRECISION) << value;
  return output.str();
}

using NRG::Tools::files_refer_to_same_location;
using NRG::Tools::finish_output;

inline auto parse_finite_double(const std::string_view value, const std::string_view name) {
  const std::string text(value);
  char *end = nullptr;
  errno = 0;
  const auto result = std::strtod(text.c_str(), &end);
  const bool underflowed_to_zero = errno == ERANGE && result == 0.0;
  if (underflowed_to_zero || end == text.c_str() || end != text.c_str() + text.size() || !std::isfinite(result))
    throw std::runtime_error(std::string(name) + " must be a finite representable number: " + text);
  return result;
}

inline auto minimum_safe_imaginary_part() { return std::sqrt(std::numeric_limits<double>::min()); }

inline void validate_integration_settings(const double lim_direct, const double epsabs, const double epsrel) {
  if (!std::isfinite(lim_direct) || lim_direct < 0.0) throw std::invalid_argument("Direct-integration threshold must be finite and nonnegative.");
  if (!std::isfinite(epsabs) || epsabs < 0.0) throw std::invalid_argument("Absolute integration tolerance must be finite and nonnegative.");
  if (!std::isfinite(epsrel) || epsrel < 0.0) throw std::invalid_argument("Relative integration tolerance must be finite and nonnegative.");
  const auto minimum_relative_tolerance = std::max(50.0 * std::numeric_limits<double>::epsilon(), 0.5e-28);
  if (epsabs == 0.0 && epsrel < minimum_relative_tolerance)
    throw std::invalid_argument("Relative integration tolerance is too small when absolute tolerance is zero.");
}

inline auto parse_energy_power(const std::string_view value) {
  int result = 0;
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
  if (error != std::errc{} || end != value.data() + value.size() || result < 0)
    throw std::runtime_error("Energy power must be a nonnegative integer: " + std::string(value));
  return result;
}

template <size_t N>
struct numeric_row {
  std::array<double, N> values;
  size_t line_number;
};

template <size_t N>
class numeric_row_reader {
  private:
  std::istream &input;
  std::string source;
  size_t line_number = 0;
  size_t record_count = 0;

  public:
  numeric_row_reader(std::istream &input_, std::string source_) : input{input_}, source{std::move(source_)} {}

  auto next() -> std::optional<numeric_row<N>> {
    std::string line;
    while (std::getline(input, line)) {
      ++line_number;
      std::istringstream fields(line);
      std::vector<std::string> tokens;
      std::string token;
      while (fields >> token) tokens.push_back(token);
      if (tokens.empty()) continue;
      if (tokens.size() != N)
        throw std::runtime_error(source + ":" + std::to_string(line_number) + ": expected exactly " + std::to_string(N) + " numeric fields; found "
                                 + std::to_string(tokens.size()) + ".");

      numeric_row<N> row{{}, line_number};
      for (size_t i = 0; i < N; ++i) {
        const auto context = source + ":" + std::to_string(line_number) + ": field " + std::to_string(i + 1);
        row.values[i] = parse_finite_double(tokens[i], context);
      }
      ++record_count;
      return row;
    }
    if (input.bad()) throw std::runtime_error(source + ": I/O error after line " + std::to_string(line_number) + ".");
    return std::nullopt;
  }

  auto records_read() const noexcept { return record_count; }
  const auto &source_name() const noexcept { return source; }
};

template <typename F>
struct gsl_callback_context {
  F *function;
  std::exception_ptr failure;
};

// Exceptions must not unwind through GSL's C frames.
template <typename F>
inline double unwrap(const double x, void *raw_context) noexcept {
  auto &context = *static_cast<gsl_callback_context<F> *>(raw_context);
  if (context.failure) return std::numeric_limits<double>::quiet_NaN();
  try {
    return (*context.function)(x);
  } catch (...) {
    context.failure = std::current_exception();
    return std::numeric_limits<double>::quiet_NaN();
  }
}

struct gsl_failure {
  int status;
  double lower;
  double upper;
  double result;
  double estimated_error;
  double epsabs;
  double epsrel;
  bool nonfinite;
};

struct integration_attempt {
  int status;
  double result;
  double estimated_error;
};

class gsl_failure_summary {
  private:
  size_t total = 0;
  size_t nonfinite = 0;
  std::map<int, size_t> by_status;
  std::optional<gsl_failure> first;
  size_t global_targets_missed = 0;
  std::optional<std::string> first_global_failure;
  bool reported = false;

  public:
  void record(const gsl_failure &failure) {
    ++total;
    if (failure.status != GSL_SUCCESS) ++by_status[failure.status];
    if (failure.nonfinite) ++nonfinite;
    if (!first) first = failure;
  }

  auto count() const noexcept { return total; }

  void record_global_failure(std::string message) {
    ++global_targets_missed;
    if (!first_global_failure) first_global_failure = std::move(message);
  }

  void merge(const gsl_failure_summary &other) {
    total += other.total;
    nonfinite += other.nonfinite;
    for (const auto &[status, count] : other.by_status) by_status[status] += count;
    if (!first && other.first) first = other.first;
    global_targets_missed += other.global_targets_missed;
    if (!first_global_failure && other.first_global_failure) first_global_failure = other.first_global_failure;
  }

  void report(std::ostream &out) {
    if (reported || (total == 0 && global_targets_missed == 0)) return;
    reported = true;
    const auto old_precision = out.precision();
    out << std::setprecision(OUTPUT_PRECISION);
    if (total != 0) {
      out << "hilb: warning: " << total << " GSL integration call(s) reported failure";
      if (!by_status.empty()) {
        out << "; statuses: ";
        bool separator = false;
        for (const auto &[status, count] : by_status) {
          if (separator) out << ", ";
          out << status << " (" << gsl_strerror(status) << ") x" << count;
          separator = true;
        }
      }
      if (nonfinite != 0) out << "; nonfinite result/error x" << nonfinite;
      if (first) {
        out << "; first: status=" << first->status << " (" << gsl_strerror(first->status) << "), interval=[" << first->lower << ',' << first->upper
            << "], result=" << first->result << ", estimated_error=" << first->estimated_error << ", epsabs=" << first->epsabs << ", epsrel=" << first->epsrel;
      }
    }
    if (global_targets_missed != 0) {
      if (total != 0) out << "; ";
      else out << "hilb: warning: ";
      out << global_targets_missed << " tabulated QAG point(s) missed the global error target";
      if (first_global_failure) out << "; first: " << *first_global_failure;
    }
    out << '\n';
    out.precision(old_precision);
  }
};

struct gsl_workspace_deleter {
  void operator()(gsl_integration_workspace *workspace) const {
    if (workspace) gsl_integration_workspace_free(workspace);
  }
};

// Wrap around GSL integration routines
class integrator {
  private:
  size_t limit;
  NRG::Tools::QagRule rule;
  NRG::Tools::GslErrorPolicy error_policy;
  bool record_failures;
  gsl_failure_summary *failures;
  std::unique_ptr<gsl_integration_workspace, gsl_workspace_deleter> work;

  template <typename F>
  auto attempt_without_handler_guard(F f, const double a, const double b, const double epsabs,
                                     const double epsrel) {
    if (!work) throw std::runtime_error("Cannot use a moved-from GSL integrator.");
    gsl_callback_context<F> context{&f, {}};
    gsl_function function{&unwrap<F>, &context};
    double result = std::numeric_limits<double>::quiet_NaN();
    double error = std::numeric_limits<double>::quiet_NaN();
    const auto status = gsl_integration_qag(&function, a, b, epsabs, epsrel, limit,
                                             NRG::Tools::gsl_qag_rule(rule), work.get(), &result, &error);
    if (context.failure) std::rethrow_exception(context.failure);
    return integration_attempt{status, result, error};
  }

  template <typename F>
  auto integrate_without_handler_guard(F f, const double a, const double b, const double epsabs,
                                       const double epsrel) {
    const auto attempt = attempt_without_handler_guard(std::move(f), a, b, epsabs, epsrel);
    const auto status = attempt.status;
    const auto result = attempt.result;
    const auto error = attempt.estimated_error;
    const bool nonfinite = !std::isfinite(result) || !std::isfinite(error);
    const bool failed = NRG::Tools::gsl_integration_failed(status, result, error);
    if (failed && record_failures && failures)
      failures->record({status, a, b, result, error, epsabs, epsrel, nonfinite});
    if (failed && error_policy == NRG::Tools::GslErrorPolicy::warn && !record_failures)
      std::cerr << "hilb: warning: qag error: " << status << " -- " << gsl_strerror(status)
                << (nonfinite ? "; nonfinite result or error estimate" : "") << '\n';
    if (failed && error_policy == NRG::Tools::GslErrorPolicy::fail)
      throw std::runtime_error(status == GSL_SUCCESS
                                 ? "qag returned a nonfinite result or error estimate"
                                 : "qag error: " + std::to_string(status) + " -- " + gsl_strerror(status));
    return result;
  }

  static auto allocate_workspace(const size_t limit) {
    NRG::Tools::validate_qag_workspace_limit(limit);
    const NRG::Tools::GslErrorHandlerGuard error_handler;
    auto *workspace = gsl_integration_workspace_alloc(limit);
    if (!workspace) throw std::runtime_error("Failed to allocate GSL integration workspace.");
    return std::unique_ptr<gsl_integration_workspace, gsl_workspace_deleter>{workspace};
  }

  public:
  struct configured_t {};
  inline static constexpr configured_t configured{};
  struct handler_already_disabled_t {};
  inline static constexpr handler_already_disabled_t handler_already_disabled{};

  integrator(size_t limit_ = 1000, bool throw_on_error_ = false, gsl_failure_summary *failures_ = nullptr)
      : limit{limit_}, rule{NRG::Tools::QagRule::gauss15},
        error_policy{throw_on_error_ ? NRG::Tools::GslErrorPolicy::fail
                                     : (failures_ ? NRG::Tools::GslErrorPolicy::warn : NRG::Tools::GslErrorPolicy::ignore)},
        record_failures{failures_ != nullptr}, failures{failures_}, work{allocate_workspace(limit)} {}

  integrator(configured_t, const size_t limit_, const NRG::Tools::QagRule rule_,
             const NRG::Tools::GslErrorPolicy error_policy_ = NRG::Tools::GslErrorPolicy::ignore,
             gsl_failure_summary *failures_ = nullptr)
      : limit{limit_}, rule{rule_}, error_policy{error_policy_},
        record_failures{error_policy == NRG::Tools::GslErrorPolicy::warn && failures_ != nullptr}, failures{failures_},
        work{allocate_workspace(limit)} {}

  integrator(const integrator &other)
      : limit{other.limit}, rule{other.rule}, error_policy{other.error_policy}, record_failures{other.record_failures}, failures{other.failures},
        work{allocate_workspace(limit)} {}

  integrator(integrator &&other) noexcept = default;

  integrator &operator=(const integrator &other) {
    if (this == &other) return *this;
    integrator copy(other);
    swap(copy);
    return *this;
  }

  integrator &operator=(integrator &&other) noexcept = default;
  ~integrator() = default;

  void set_failure_summary(gsl_failure_summary *failures_) {
    failures = failures_;
    record_failures = error_policy == NRG::Tools::GslErrorPolicy::warn && failures != nullptr;
  }

  auto policy() const noexcept { return error_policy; }

  void record_warning(const gsl_failure &failure) {
    if (error_policy != NRG::Tools::GslErrorPolicy::warn) return;
    if (failures) {
      failures->record(failure);
      return;
    }
    std::cerr << "hilb: warning: qag error: " << failure.status << " -- " << gsl_strerror(failure.status)
              << (failure.nonfinite ? "; nonfinite result or error estimate" : "") << '\n';
  }

  void record_global_warning(std::string message) {
    if (error_policy != NRG::Tools::GslErrorPolicy::warn) return;
    if (failures) {
      failures->record_global_failure(std::move(message));
      return;
    }
    std::cerr << "hilb: warning: " << message << '\n';
  }

  void swap(integrator &other) noexcept {
    using std::swap;
    swap(limit, other.limit);
    swap(rule, other.rule);
    swap(error_policy, other.error_policy);
    swap(record_failures, other.record_failures);
    swap(failures, other.failures);
    swap(work, other.work);
  }

  /**
     * Integrate function f on [a:b].
     *
     * @param f Function to be integrated
     * @param a Lower integration range boundary
     * @param b Upper integration range boundary
     * @param epsabs numeric integration epsilon (absolute)
     * @param epsrel numeric integration epsilon (relative)
     */
  template <typename F>
  auto operator()(F f, const double a, const double b, const double epsabs = 1e-14, const double epsrel = 1e-10) {
    const NRG::Tools::GslErrorHandlerGuard error_handler;
    return integrate_without_handler_guard(std::move(f), a, b, epsabs, epsrel);
  }

  template <typename F>
  auto operator()([[maybe_unused]] handler_already_disabled_t tag, F f, const double a, const double b,
                  const double epsabs = 1e-14, const double epsrel = 1e-10) {
    return integrate_without_handler_guard(std::move(f), a, b, epsabs, epsrel);
  }

  template <typename F>
  auto attempt([[maybe_unused]] handler_already_disabled_t tag, F f, const double a, const double b,
               const double epsabs = 1e-14, const double epsrel = 1e-10) {
    return attempt_without_handler_guard(std::move(f), a, b, epsabs, epsrel);
  }
};

// Wrap around GSL interpolation routines
class interpolator {
  private:
  NRG::Tools::PiecewisePolynomial<double> polynomial;
  double oob_value;

  public:
  interpolator(const std::vector<double> &_X, const std::vector<double> &_Y, const double _oob_value = 0.0)
      : polynomial{NRG::Tools::make_gsl_piecewise_polynomial(_X, _Y, NRG::Tools::InterpolationMethod::cspline)}, oob_value{_oob_value} {}
  interpolator(const std::vector<double> &_X, const std::vector<double> &_Y, const double _oob_value,
                const NRG::Tools::InterpolationMethod method_)
      : polynomial{NRG::Tools::make_gsl_piecewise_polynomial(_X, _Y, method_)}, oob_value{_oob_value} {}
  auto operator()(const double x) const {
    return polynomial.lower_bound() <= x && x <= polynomial.upper_bound() ? polynomial.evaluate(x) : oob_value;
  }
  const auto &piecewise_polynomial() const noexcept { return polynomial; }
  auto out_of_bounds_value() const noexcept { return oob_value; }
};

// Square of x
inline auto sqr(const double x) { return x * x; }

struct scaled_value {
  double fraction = 0.0;
  int exponent = 0;
};

struct scaled_complex_value {
  scaled_value real;
  scaled_value imaginary;
};

inline auto make_scaled_product_ratio(const double value, const double numerator, const double denominator) {
  if (value == 0.0 || numerator == 0.0) return scaled_value{std::copysign(0.0, value * numerator), 0};
  int value_exponent = 0;
  int numerator_exponent = 0;
  int denominator_exponent = 0;
  const auto value_fraction = std::frexp(value, &value_exponent);
  const auto numerator_fraction = std::frexp(numerator, &numerator_exponent);
  const auto denominator_fraction = std::frexp(denominator, &denominator_exponent);
  int normalization_exponent = 0;
  const auto fraction = std::frexp(value_fraction * numerator_fraction / denominator_fraction,
                                   &normalization_exponent);
  return scaled_value{fraction, value_exponent + numerator_exponent - denominator_exponent + normalization_exponent};
}

inline auto add_scaled_values(const scaled_value first, const scaled_value second) {
  if (first.fraction == 0.0) return second;
  if (second.fraction == 0.0) return first;
  const auto exponent = std::max(first.exponent, second.exponent);
  const auto sum = std::scalbn(first.fraction, first.exponent - exponent)
                   + std::scalbn(second.fraction, second.exponent - exponent);
  if (sum == 0.0) return scaled_value{std::copysign(0.0, sum), 0};
  int normalization_exponent = 0;
  const auto fraction = std::frexp(sum, &normalization_exponent);
  return scaled_value{static_cast<double>(fraction), exponent + normalization_exponent};
}

inline auto scale_value(const scaled_value value, const double numerator, const double denominator = 1.0) {
  const auto scaled = make_scaled_product_ratio(value.fraction, numerator, denominator);
  return scaled_value{scaled.fraction, scaled.exponent + value.exponent};
}

inline auto multiply_scaled_values(const scaled_value first, const scaled_value second) {
  if (first.fraction == 0.0 || second.fraction == 0.0) return scaled_value{};
  int normalization_exponent = 0;
  const auto fraction = std::frexp(first.fraction * second.fraction, &normalization_exponent);
  return scaled_value{fraction, first.exponent + second.exponent + normalization_exponent};
}

inline auto scaled_exponential(const double logarithm) {
  constexpr auto logarithm_of_two = 0.693147180559945309417232121458176568L;
  const auto exponent = static_cast<int>(std::floor(static_cast<long double>(logarithm) / logarithm_of_two));
  int normalization_exponent = 0;
  const auto fraction = std::frexp(std::exp(static_cast<long double>(logarithm)
                                            - static_cast<long double>(exponent) * logarithm_of_two),
                                   &normalization_exponent);
  return scaled_value{static_cast<double>(fraction), exponent + normalization_exponent};
}

inline auto scaled_to_double(const scaled_value value) { return std::scalbn(value.fraction, value.exponent); }

inline auto robust_complex_divide_scaled(const std::complex<double> numerator, const double denominator_real,
                                         const double denominator_imaginary) {
  const auto numerator_scale = std::max(std::abs(numerator.real()), std::abs(numerator.imag()));
  if (numerator_scale == 0.0) return scaled_complex_value{};
  const auto denominator_scale = std::max(std::abs(denominator_real), std::abs(denominator_imaginary));
  const auto denominator_real_scaled = denominator_real / denominator_scale;
  const auto denominator_imaginary_scaled = denominator_imaginary / denominator_scale;
  const auto norm = sqr(denominator_real_scaled) + sqr(denominator_imaginary_scaled);
  const auto real_first = make_scaled_product_ratio(numerator.real(), denominator_real_scaled / norm,
                                                     denominator_scale);
  const auto real_second = make_scaled_product_ratio(numerator.imag(), denominator_imaginary_scaled / norm,
                                                      denominator_scale);
  const auto imaginary_first = make_scaled_product_ratio(numerator.imag(), denominator_real_scaled / norm,
                                                          denominator_scale);
  const auto imaginary_second = make_scaled_product_ratio(numerator.real(), -denominator_imaginary_scaled / norm,
                                                           denominator_scale);
  const auto real = add_scaled_values(real_first, real_second);
  const auto imaginary = add_scaled_values(imaginary_first, imaginary_second);
  return scaled_complex_value{real, imaginary};
}

inline auto robust_complex_divide(const std::complex<double> numerator, const double denominator_real,
                                  const double denominator_imaginary) {
  const auto result = robust_complex_divide_scaled(numerator, denominator_real, denominator_imaginary);
  return std::complex<double>{scaled_to_double(result.real), scaled_to_double(result.imaginary)};
}

inline auto exact_scaled_complex_divide(const std::complex<double> numerator, const double argument_real,
                                        const double energy, const double argument_imaginary,
                                        const double multiplier = 1.0) {
  const auto real = NRG::Tools::detail::exact_rational(argument_real)
                    - NRG::Tools::detail::exact_rational(energy);
  const auto imaginary = NRG::Tools::detail::exact_rational(argument_imaginary);
  const auto norm = real * real + imaginary * imaginary;
  const auto exact_multiplier = NRG::Tools::detail::exact_rational(multiplier);
  const NRG::Tools::detail::ExactRational numerator_real =
    NRG::Tools::detail::exact_rational(numerator.real()) * exact_multiplier;
  const NRG::Tools::detail::ExactRational numerator_imaginary =
    NRG::Tools::detail::exact_rational(numerator.imag()) * exact_multiplier;
  const NRG::Tools::detail::ExactRational result_real =
    (numerator_real * real + numerator_imaginary * imaginary) / norm;
  const NRG::Tools::detail::ExactRational result_imaginary =
    (numerator_imaginary * real - numerator_real * imaginary) / norm;
  return NRG::Tools::detail::ExactComplex{result_real, result_imaginary};
}

inline auto logarithm_of_absolute_sum(const double first, const double second) {
  if (first == 0.0) return std::log(std::abs(second));
  if (second == 0.0) return std::log(std::abs(first));
  if (std::signbit(first) != std::signbit(second)) return std::log(std::abs(first + second));
  const auto scale = std::max(std::abs(first), std::abs(second));
  return std::log(scale) + std::log(std::abs(first / scale + second / scale));
}

// Result of Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}] (atg -> imQ).
inline auto imQ(const double x, const double y, const double B) { return std::atan((-B + x) / y) - std::atan((B + x) / y); }

// Result of Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}] (logs -> reQ).
inline auto reQ(const double x, const double y, const double B) {
  const long double scale = std::max({std::abs(static_cast<long double>(B)), std::abs(static_cast<long double>(x)),
                                      std::abs(static_cast<long double>(y))});
  const auto bandwidth = static_cast<long double>(B) / scale;
  const auto real = static_cast<long double>(x) / scale;
  const auto negative = bandwidth - real;
  const auto imaginary = static_cast<long double>(y) / scale;
  const auto denominator = negative * negative + imaginary * imaginary;
  const auto difference = 4.0L * bandwidth * real;
  const auto ratio = difference / denominator;
  // Rounding in a ratio near -1 is strongly amplified by log1p; use the log-space form below instead.
  if (std::isfinite(ratio) && ratio >= -0.5L) return static_cast<double>(0.5L * std::log1p(ratio));

  auto logarithmic_hypot = [y](const double first, const double second) {
    const auto real_logarithm = logarithm_of_absolute_sum(first, second);
    const auto imaginary_logarithm = std::log(std::abs(y));
    return 0.5L * NRG::Tools::detail::log_add_positive(2.0L * real_logarithm,
                                                       2.0L * imaginary_logarithm);
  };
  return static_cast<double>(logarithmic_hypot(B, x) - logarithmic_hypot(B, -x));
}

inline auto stable_interval_logarithm(const double left, const double right, const std::complex<double> argument) {
  if constexpr (NRG::Tools::detail::native_extended_precision) {
    const auto width = static_cast<long double>(right) - static_cast<long double>(left);
    const auto left_real_delta = static_cast<long double>(argument.real()) - static_cast<long double>(left);
    const auto right_real_delta = static_cast<long double>(argument.real()) - static_cast<long double>(right);
    const auto imaginary_delta = static_cast<long double>(argument.imag());
    const auto left_norm = left_real_delta * left_real_delta + imaginary_delta * imaginary_delta;
    const auto right_norm = right_real_delta * right_real_delta + imaginary_delta * imaginary_delta;
    const auto use_left = left_norm < right_norm;
    const auto real_delta = use_left ? left_real_delta : right_real_delta;
    const auto norm = real_delta * real_delta + imaginary_delta * imaginary_delta;
    const auto numerator = use_left ? -width : width;
    const auto quotient_real = numerator * real_delta / norm;
    const auto quotient_imaginary = -numerator * imaginary_delta / norm;
    const auto norm_difference = width * (2.0L * static_cast<long double>(argument.real())
                                          - (static_cast<long double>(left) + static_cast<long double>(right)));
    const auto smaller_norm = std::min(left_norm, right_norm);
    const auto larger_norm = std::max(left_norm, right_norm);
    const auto logarithm_real = smaller_norm >= larger_norm / 2.0L
                                  ? 0.5L * std::log1p(norm_difference / right_norm)
                                  : 0.5L * (std::log(left_norm) - std::log(right_norm));
    const auto logarithm_imaginary = std::atan2(quotient_imaginary, 1.0L + quotient_real);
    const auto direction = use_left ? -1.0L : 1.0L;
    return std::complex<long double>{logarithm_real, direction * logarithm_imaginary};
  } else {
    using Wide = NRG::Tools::PiecewiseWideFloat;
    const auto width = Wide{right} - Wide{left};
    const auto left_real_delta = Wide{argument.real()} - Wide{left};
    const auto right_real_delta = Wide{argument.real()} - Wide{right};
    const auto imaginary_delta = Wide{argument.imag()};
    const auto left_norm = left_real_delta * left_real_delta + imaginary_delta * imaginary_delta;
    const auto right_norm = right_real_delta * right_real_delta + imaginary_delta * imaginary_delta;
    const auto use_left = left_norm < right_norm;
    const auto real_delta = use_left ? left_real_delta : right_real_delta;
    const auto norm = real_delta * real_delta + imaginary_delta * imaginary_delta;
    const auto numerator = use_left ? -width : width;
    const auto quotient_real = numerator * real_delta / norm;
    const auto quotient_imaginary = -numerator * imaginary_delta / norm;
    const auto norm_difference = width * (Wide{2} * Wide{argument.real()} - (Wide{left} + Wide{right}));
    const auto smaller_norm = std::min(left_norm, right_norm);
    const auto larger_norm = std::max(left_norm, right_norm);
    const auto logarithm_real = smaller_norm >= larger_norm / Wide{2}
                                  ? Wide{0.5} * boost::multiprecision::log1p(norm_difference / right_norm)
                                  : Wide{0.5} * (boost::multiprecision::log(left_norm)
                                                 - boost::multiprecision::log(right_norm));
    const auto logarithm_imaginary = boost::multiprecision::atan2(quotient_imaginary, Wide{1} + quotient_real);
    const auto direction = use_left ? Wide{-1} : Wide{1};
    return NRG::Tools::PiecewiseWideComplex{logarithm_real, direction * logarithm_imaginary};
  }
}

template <bool NativeExtendedPrecision = NRG::Tools::detail::native_extended_precision>
inline auto scaled_interval_logarithm(const long double scale, const double left, const double right,
                                      const std::complex<double> argument) {
  const auto logarithm = stable_interval_logarithm(left, right, argument);
  if constexpr (NativeExtendedPrecision) {
    return scale * logarithm;
  } else {
    const NRG::Tools::PiecewiseWideFloat wide_scale{scale};
    const NRG::Tools::PiecewiseWideComplex wide_result{wide_scale * logarithm.real(),
                                                        wide_scale * logarithm.imag()};
    return std::complex<long double>{wide_result.real().template convert_to<long double>(),
                                     wide_result.imag().template convert_to<long double>()};
  }
}

struct zero_density {
  auto operator()([[maybe_unused]] const double energy) const { return 0.0; }
};

// Calculate the (half)bandwidth, i.e., the size B of the enclosing interval [-B:B].
inline auto bandwidth(const std::vector<double> &X) {
  assert(std::is_sorted(X.begin(), X.end()));
  const auto Xmin = X.front();
  const auto Xmax = X.back();
  return std::max(std::abs(Xmin), std::abs(Xmax));
}

inline void validate_hilbert_transform_argument(const double B, const std::complex<double> z, const int n) {
  if (n < 0) throw std::invalid_argument("Energy power must be nonnegative.");
  if (!std::isfinite(B) || B <= 0.0) throw std::invalid_argument("Half-bandwidth must be finite and positive.");
  if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) throw std::invalid_argument("Hilbert-transform argument must be finite.");
  if (std::abs(z.imag()) < minimum_safe_imaginary_part())
    throw std::invalid_argument("Absolute imaginary part is below the minimum safe value sqrt(numeric_limits<double>::min()).");
}

inline void validate_hilbert_transform_inputs(const double B, const std::complex<double> z, const double lim_direct,
                                               const int n, const double epsabs, const double epsrel) {
  validate_hilbert_transform_argument(B, z, n);
  validate_integration_settings(lim_direct, epsabs, epsrel);
}

template <typename Scalar>
auto hilbert_transform(const NRG::Tools::PiecewisePolynomial<Scalar> &density, const double B, const std::complex<double> z,
                        const double lim_direct = 1e-3, const int n = 0, const double epsabs = 1e-14,
                        const double epsrel = 1e-10) {
  validate_hilbert_transform_inputs(B, z, lim_direct, n, epsabs, epsrel);
  if (density.lower_bound() < -B || density.upper_bound() > B)
    throw std::invalid_argument("Piecewise-polynomial support must lie within [-B, B].");
  return NRG::Tools::cauchy_transform(density.multiply_by_monomial(static_cast<size_t>(n)), z);
}

/**
 * Calculate the Hilbert transform of a given spectral function at fixed complex value z. This is the low-level routine, called from
 * other interfaces. The general strategy is to perform a direct integration of the defining integral for cases where z has a sufficiently
 * large imaginary part, otherwise the singularity is subtracted out and the calculation of the transformed integrand is performed using
 * an integration-variable substitution to better handle small values.
 *
 * @param integration Reusable numerical integration workspace
 * @param rhor Real part of the spectral function
 * @param rhoi Imaginary part of the spectral function
 * @param B Half-bandwidth, i.e., the support of the spectral function is [-B:B]
 * @param z The complex value for which to evaluate the Hilbert transform
 * @param lim_direct value of y=Im(z) above which g(E)/(x+Iy-E) is directly integrated, and below which the singularity is removed
 * @param n Nonnegative integer power of the integration energy
 * @param epsabs Absolute integration tolerance
 * @param epsrel Relative integration tolerance
 */
template <typename FNCR, typename FNCI>
auto hilbert_transform_without_handler_guard(integrator &integration, FNCR &rhor, FNCI &rhoi, const double B,
                                              const std::complex<double> z, const double lim_direct = 1e-3,
                                              const int n = 0, const double epsabs = 1e-14,
                                              const double epsrel = 1e-10) {
  validate_hilbert_transform_inputs(B, z, lim_direct, n, epsabs, epsrel);
  auto integrate = [&integration](auto function, const double lower, const double upper,
                                  const double absolute_tolerance, const double relative_tolerance) {
    return integration(integrator::handler_already_disabled, std::move(function), lower, upper,
                       absolute_tolerance, relative_tolerance);
  };
  const auto x = real(z);
  const auto y = imag(z);
  auto g_r = [&rhor, n](const double omega) -> double {
    return n == 0 ? rhor(omega) : std::pow(omega, n) * rhor(omega);
  };
  auto g_i = [&rhoi, n](const double omega) -> double {
    return n == 0 ? rhoi(omega) : std::pow(omega, n) * rhoi(omega);
  };

  const bool inside = -B <= x && x <= B;
  const auto logarithmic_y = std::log(std::abs(y));
  const bool unusable_exterior_limits = !inside && [&]() {
    if (!std::isfinite(x + B) || !std::isfinite(x - B)) return true;
    const auto near = x > B ? logarithm_of_absolute_sum(x, -B) : logarithm_of_absolute_sum(x, B);
    const auto far = x > B ? logarithm_of_absolute_sum(x, B) : logarithm_of_absolute_sum(x, -B);
    return !(near < far);
  }();
  if (std::abs(y) >= lim_direct || unusable_exterior_limits) {
    auto direct_value = [x, y](const std::complex<double> numerator, const double omega) {
      if (x == 0.0 || omega == 0.0 || std::signbit(x) == std::signbit(omega))
        return robust_complex_divide_scaled(numerator, x - omega, y);

      // Opposite signs can make x-omega overflow, so form a scaled denominator first.
      const auto scale = std::max({std::abs(x), std::abs(omega), std::abs(y)});
      const auto quotient = robust_complex_divide_scaled(numerator, x / scale - omega / scale, y / scale);
      return scaled_complex_value{scale_value(quotient.real, 1.0, scale),
                                  scale_value(quotient.imaginary, 1.0, scale)};
    };
    auto integrate_density = [&integrate, B, epsabs, epsrel, &direct_value](auto &density,
                                                                            const bool imaginary) {
        auto integrate_component = [&](const bool output_imaginary) {
          auto integrand = [B, imaginary, output_imaginary, &density, &direct_value](const double unit) {
            const auto energy = B * unit;
            const auto density_value = density(energy);
            const auto numerator = imaginary ? std::complex<double>{0.0, density_value}
                                             : std::complex<double>{density_value, 0.0};
            const auto value = direct_value(numerator, energy);
            return scaled_to_double(scale_value(output_imaginary ? value.imaginary : value.real, B));
          };
          return integrate(integrand, -1.0, 1.0, epsabs, epsrel);
        };
        return std::complex<double>{integrate_component(false), integrate_component(true)};
    };
    const auto real_density_result = integrate_density(g_r, false);
    if constexpr (std::is_same_v<std::remove_cvref_t<FNCI>, zero_density>) return real_density_result;
    const auto imaginary_density_result = integrate_density(g_i, true);
    NRG::Tools::CompensatedComplexSum result;
    result.add(real_density_result);
    result.add(imaginary_density_result);
    auto combined = result.value();
    auto cancelled = [](const double first, const double second, const double value) {
        return first != 0.0 && second != 0.0 && std::signbit(first) != std::signbit(second)
               && std::abs(value) <= std::sqrt(std::numeric_limits<double>::epsilon())
                                        * (std::abs(first) + std::abs(second));
    };
    const auto recompute_real = cancelled(real_density_result.real(), imaginary_density_result.real(), combined.real());
    const auto recompute_imaginary = cancelled(real_density_result.imag(), imaginary_density_result.imag(),
                                                combined.imag());
    auto integrate_exact_component = [&](const bool imaginary) {
        auto integrand = [x, y, B, imaginary, &g_r, &g_i](const double unit) {
          const auto positive_energy = B * unit;
          const auto negative_energy = -positive_energy;
          const auto positive = exact_scaled_complex_divide(
            {g_r(positive_energy), g_i(positive_energy)}, x, positive_energy, y, B);
          const auto negative = exact_scaled_complex_divide(
            {g_r(negative_energy), g_i(negative_energy)}, x, negative_energy, y, B);
          const NRG::Tools::detail::ExactRational value = imaginary
                                                            ? positive.imaginary + negative.imaginary
                                                            : positive.real + negative.real;
          return value.convert_to<double>();
        };
        return integrate(integrand, 0.0, 1.0, epsabs, epsrel);
    };
    if (recompute_real) combined.real(integrate_exact_component(false));
    if (recompute_imaginary) combined.imag(integrate_exact_component(true));
    return combined;
  }

  const double g_r_x = inside ? g_r(x) : 0.0;
  const double g_i_x = inside ? g_i(x) : 0.0;

  const auto logarithmic_right_distance = logarithm_of_absolute_sum(B, -x) - logarithmic_y;
  const auto logarithmic_left_distance = logarithm_of_absolute_sum(B, x) - logarithmic_y;

  // Integrate the singularity-subtracted remainder after the logarithmic change of variables.
  auto calc_subtracted = [&integrate, x, B, inside, logarithmic_right_distance, logarithmic_left_distance,
                          logarithmic_y, epsabs, epsrel](auto positive, auto negative, auto analytic) -> double {
    double lim1down = 1.0, lim1up = -1.0, lim2down = 1.0, lim2up = -1.0;
    if (inside) {
      constexpr double log_zero_cutoff = -36.8; // approximately log(10^-16)
      if (x < B) {
        lim1down = log_zero_cutoff;
        lim1up = logarithmic_right_distance;
      }
      if (x > -B) {
        lim2down = log_zero_cutoff;
        lim2up = logarithmic_left_distance;
      }
    } else if (x > B) {
      lim2down = logarithm_of_absolute_sum(x, -B) - logarithmic_y;
      lim2up = logarithm_of_absolute_sum(x, B) - logarithmic_y;
    } else {
      lim1down = logarithm_of_absolute_sum(x, B) - logarithmic_y;
      lim1up = logarithm_of_absolute_sum(x, -B) - logarithmic_y;
    }
    const auto result1 = lim1down < lim1up ? integrate(positive, lim1down, lim1up, epsabs, epsrel) : 0.0;
    const auto result2 = lim2down < lim2up ? integrate(negative, lim2down, lim2up, epsabs, epsrel) : 0.0;
    const auto result3 = inside ? analytic : 0.0;
    return result1 + result2 + result3;
  };

  auto subtracted_value = [x, y, g_r_x, g_i_x, &g_r, &g_i](const double omega) {
    return robust_complex_divide_scaled({g_r(omega) - g_r_x, g_i(omega) - g_i_x}, x - omega, y);
  };
  auto transformed_value = [x, B, logarithmic_y, &subtracted_value](const double logarithmic_displacement,
                                                                    const bool positive,
                                                                    const bool imaginary) {
    const auto displacement = scaled_exponential(logarithmic_y + logarithmic_displacement);
    const auto narrow_displacement = scaled_to_double(displacement);
    double omega = 0.0;
    if (std::isfinite(narrow_displacement)) {
      omega = positive ? std::min(B, x + narrow_displacement) : std::max(-B, x - narrow_displacement);
    } else {
      const auto wide_displacement = boost::multiprecision::ldexp(NRG::Tools::detail::FarReal{displacement.fraction},
                                                                  displacement.exponent);
      const auto wide_omega = NRG::Tools::detail::FarReal{x}
                              + (positive ? wide_displacement : -wide_displacement);
      omega = std::clamp(wide_omega.convert_to<double>(), -B, B);
    }
    const auto value = subtracted_value(omega);
    return scaled_to_double(multiply_scaled_values(imaginary ? value.imaginary : value.real, displacement));
  };
  auto ref3p = [&transformed_value](const double r) { return transformed_value(r, true, false); };
  auto ref3m = [&transformed_value](const double r) { return transformed_value(r, false, false); };
  const auto red = g_r_x * reQ(x, y, B) - g_i_x * imQ(x, y, B);

  auto imf3p = [&transformed_value](const double r) { return transformed_value(r, true, true); };
  auto imf3m = [&transformed_value](const double r) { return transformed_value(r, false, true); };
  const auto imd = g_r_x * imQ(x, y, B) + g_i_x * reQ(x, y, B);

  const auto real_result = calc_subtracted(ref3p, ref3m, red);
  const auto imag_result = calc_subtracted(imf3p, imf3m, imd);
  return std::complex(real_result, imag_result);
}

template <typename FNCR, typename FNCI>
auto hilbert_transform(integrator &integration, FNCR &rhor, FNCI &rhoi, const double B,
                       const std::complex<double> z, const double lim_direct = 1e-3, const int n = 0,
                       const double epsabs = 1e-14, const double epsrel = 1e-10) {
  const NRG::Tools::GslErrorHandlerGuard error_handler;
  return hilbert_transform_without_handler_guard(integration, rhor, rhoi, B, z, lim_direct, n, epsabs, epsrel);
}

inline auto hilbert_transform(integrator &integration, const interpolator &rhor, const interpolator &rhoi, const double B,
                              const std::complex<double> z, const double lim_direct = 1e-3, const int n = 0,
                              const double epsabs = 1e-14, const double epsrel = 1e-10) {
  const auto &real_polynomial = rhor.piecewise_polynomial();
  const auto &imaginary_polynomial = rhoi.piecewise_polynomial();
  const auto exact_domain = rhor.out_of_bounds_value() == 0.0 && rhoi.out_of_bounds_value() == 0.0
                            && real_polynomial.knots() == imaginary_polynomial.knots()
                            && -B <= real_polynomial.lower_bound() && real_polynomial.upper_bound() <= B;
  if (exact_domain)
    return hilbert_transform(NRG::Tools::combine_piecewise_polynomials(real_polynomial, imaginary_polynomial), B, z,
                             lim_direct, n, epsabs, epsrel);
  auto real_callable = [&rhor](const double energy) { return rhor(energy); };
  auto imaginary_callable = [&rhoi](const double energy) { return rhoi(energy); };
  return hilbert_transform(integration, real_callable, imaginary_callable, B, z, lim_direct, n, epsabs, epsrel);
}

inline auto hilbert_transform(integrator &integration, interpolator &rhor, interpolator &rhoi, const double B,
                              const std::complex<double> z, const double lim_direct = 1e-3, const int n = 0,
                              const double epsabs = 1e-14, const double epsrel = 1e-10) {
  return hilbert_transform(integration, static_cast<const interpolator &>(rhor), static_cast<const interpolator &>(rhoi),
                           B, z, lim_direct, n, epsabs, epsrel);
}

inline auto hilbert_transform(integrator &integration, interpolator &rhor, const interpolator &rhoi, const double B,
                              const std::complex<double> z, const double lim_direct = 1e-3, const int n = 0,
                              const double epsabs = 1e-14, const double epsrel = 1e-10) {
  return hilbert_transform(integration, static_cast<const interpolator &>(rhor), rhoi, B, z, lim_direct, n, epsabs, epsrel);
}

inline auto hilbert_transform(integrator &integration, const interpolator &rhor, interpolator &rhoi, const double B,
                              const std::complex<double> z, const double lim_direct = 1e-3, const int n = 0,
                              const double epsabs = 1e-14, const double epsrel = 1e-10) {
  return hilbert_transform(integration, rhor, static_cast<const interpolator &>(rhoi), B, z, lim_direct, n, epsabs, epsrel);
}

template <typename FNCR, typename FNCI>
auto hilbert_transform(FNCR rhor, FNCI rhoi, const double B, const std::complex<double> z, const double lim_direct = 1e-3,
                        const int n = 0, const double epsabs = 1e-14, const double epsrel = 1e-10) {
  integrator integration;
  return hilbert_transform(integration, rhor, rhoi, B, z, lim_direct, n, epsabs, epsrel);
}

  /*
   * @param Xpts Frequency mesh
   * @param Rpts Real part of the spectral function
   * @param Ipts Imaginary part of the spectral function
   */ 
template <typename T>
auto hilbert_transform_with_interpolation(const T &Xpts, const T &Rpts, const T &Ipts, const std::complex<double> z,
                                          const NRG::Tools::InterpolationMethod interpolation_method,
                                          const double lim_direct = 1e-3, const int n = 0, const double epsabs = 1e-14,
                                          const double epsrel = 1e-10) {
  const auto rhor = NRG::Tools::make_gsl_piecewise_polynomial(Xpts, Rpts, interpolation_method);
  const auto rhoi = NRG::Tools::make_gsl_piecewise_polynomial(Xpts, Ipts, interpolation_method);
  const double B = bandwidth(Xpts);
  return hilbert_transform(NRG::Tools::combine_piecewise_polynomials(rhor, rhoi), B, z, lim_direct, n, epsabs, epsrel);
}

template <typename T>
auto hilbert_transform(const T &Xpts, const T &Rpts, const T &Ipts, const std::complex<double> z,
                       const double lim_direct = 1e-3, const int n = 0, const double epsabs = 1e-14,
                       const double epsrel = 1e-10) {
  return hilbert_transform_with_interpolation(Xpts, Rpts, Ipts, z, NRG::Tools::InterpolationMethod::cspline,
                                               lim_direct, n, epsabs, epsrel);
}

class Hilb {
  private:
  enum class BandwidthSource { defaults, scale, half_bandwidth };

  double scale = 1.0;         // scale factor
  double B     = 1.0 / scale; // half-bandwidth
  double Xmin = -B;
  double Xmax = +B;
  int verbosity    = 0;
  BandwidthSource bandwidth_source = BandwidthSource::defaults;
  bool G           = false; // G(z). Reports real and imaginary part.
  std::vector<double> Xpts, Ypts;
  std::optional<interpolator> tabulated_rho;
  std::optional<NRG::Tools::PiecewisePolynomial<double>> weighted_tabulated_rho;
  std::optional<double> dos_integral;
  bool tabulated = false; // Use tabulated DOS. If false, use rho_Bethe().
  double shiftx = 0.0;
  double shifty = 0.0;
  int n = 0; // power of the integration energy
  double clipping = minimum_safe_imaginary_part();
  double lim_direct = 1e-3;
  double epsabs = 1e-14;
  double epsrel = 1e-10;
  double frequency_tolerance = 1e-6;
  NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::cspline;
  Algorithm algorithm = Algorithm::qag;
  size_t workspace_limit = 1000;
  NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::warn;
  size_t jobs = 1;
  std::string jobs_source = "default";
  size_t processed_points = 0;
  size_t actual_workers = 0;
  gsl_failure_summary gsl_failures;
  std::chrono::steady_clock::time_point wall_start;
  std::clock_t cpu_start;

  struct transform_point {
    double x;
    double y;
  };

  struct transform_outcome {
    std::complex<double> value;
    gsl_failure_summary failures;
    std::exception_ptr failure;
  };

  void resolve_jobs(const std::optional<size_t> requested) {
    if (requested) {
      jobs = *requested;
      jobs_source = "--jobs";
    } else if (const auto *environment = std::getenv("OMP_NUM_THREADS")) {
      jobs = NRG::Tools::parse_worker_count(environment, "OMP_NUM_THREADS");
      jobs_source = "OMP_NUM_THREADS";
    } else {
      jobs = 1;
      jobs_source = "default";
    }
  }

  auto tabulated_qag(integrator &integration, const std::complex<double> z) const {
    if (!tabulated_rho) throw std::logic_error("Tabulated DOS interpolator is not initialized.");
    validate_hilbert_transform_inputs(B, z, lim_direct, n, epsabs, epsrel);
    const auto &polynomial = tabulated_rho->piecewise_polynomial();
    const auto &knots = polynomial.knots();
    const auto &all_coefficients = polynomial.coefficients();
    const auto support_width = static_cast<long double>(knots.back()) - static_cast<long double>(knots.front());
    NRG::Tools::CompensatedLongComplexSum result_sum;
    NRG::Tools::CompensatedLongComplexSum error_sum;
    auto failed_calls = size_t{0};
    int first_status = GSL_SUCCESS;
    auto first_lower = 0.0;
    auto first_upper = 0.0;

    for (size_t interval = 0; interval < all_coefficients.size(); ++interval) {
      const auto left = knots[interval];
      const auto right = knots[interval + 1];
      const auto extended_left = static_cast<long double>(left);
      const auto extended_width = static_cast<long double>(right) - extended_left;
      const auto width = static_cast<double>(extended_width);
      const auto &coefficients = all_coefficients[interval];
      const auto zero_interval = std::all_of(coefficients.begin(), coefficients.end(),
                                             [](const double coefficient) { return coefficient == 0.0; });
      if (zero_interval) continue;
      auto weighted_value = [&coefficients, extended_left, extended_width, power = n](const long double unit) {
        auto density = 0.0L;
        for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
          density = density * unit + static_cast<long double>(*coefficient);
        const auto energy = extended_left + extended_width * unit;
        return power == 0 ? density : density * std::pow(energy, power);
      };

      const auto regularized = std::abs(z.imag()) < lim_direct;
      const auto reference_unit = std::clamp((static_cast<long double>(z.real()) - extended_left) / extended_width,
                                             0.0L, 1.0L);
      const auto extended_reference = regularized ? weighted_value(reference_unit) : 0.0L;
      std::complex<long double> constant;
      if (regularized && extended_reference != 0.0L) {
        constant = scaled_interval_logarithm(extended_reference, left, right, z);
      }
      const auto constant_interval = n == 0
                                     && std::all_of(coefficients.begin() + 1, coefficients.end(),
                                                    [](const double coefficient) { return coefficient == 0.0; });
      if (regularized && constant_interval) {
        result_sum.add(constant);
        continue;
      }

      auto local_epsabs = epsabs;
      if (epsabs != 0.0) {
        const auto apportioned = static_cast<long double>(epsabs) * extended_width / support_width;
        local_epsabs = static_cast<double>(std::max(
          apportioned, static_cast<long double>(std::numeric_limits<double>::denorm_min())));
      }
      auto integrate_component = [&](const bool imaginary) {
        auto function = [&](const double unit) {
          const auto extended_unit = static_cast<long double>(unit);
          const auto energy = extended_left + extended_width * extended_unit;
          const auto numerator = static_cast<double>(weighted_value(extended_unit) - extended_reference);
          const auto narrow_energy = static_cast<double>(energy);
          scaled_complex_value value;
          if (z.real() == 0.0 || narrow_energy == 0.0
              || std::signbit(z.real()) == std::signbit(narrow_energy)) {
            const auto denominator_real = static_cast<double>(
              (static_cast<long double>(z.real()) - extended_left) - extended_width * extended_unit);
            value = robust_complex_divide_scaled({numerator, 0.0}, denominator_real, z.imag());
          } else {
            const auto scale = std::max({std::abs(z.real()), std::abs(narrow_energy), std::abs(z.imag())});
            const auto quotient = robust_complex_divide_scaled({numerator, 0.0},
                                                                 z.real() / scale - narrow_energy / scale,
                                                                 z.imag() / scale);
            value = scaled_complex_value{scale_value(quotient.real, 1.0, scale),
                                         scale_value(quotient.imaginary, 1.0, scale)};
          }
          return scaled_to_double(scale_value(imaginary ? value.imaginary : value.real, width));
        };
        const auto attempt = integration.attempt(integrator::handler_already_disabled, function, 0.0, 1.0,
                                                  local_epsabs, epsrel);
        if (NRG::Tools::gsl_integration_failed(attempt.status, attempt.result, attempt.estimated_error)) {
          if (failed_calls == 0) {
            first_status = attempt.status;
            first_lower = left;
            first_upper = right;
          }
          ++failed_calls;
        }
        return attempt;
      };

      const auto real = integrate_component(false);
      const auto imaginary = integrate_component(true);
      result_sum.add({static_cast<long double>(real.result) + constant.real(),
                      static_cast<long double>(imaginary.result) + constant.imag()});
      error_sum.add({std::abs(static_cast<long double>(real.estimated_error)),
                     std::abs(static_cast<long double>(imaginary.estimated_error))});
    }

    const auto extended_result = result_sum.value();
    const auto extended_error = error_sum.value();
    const std::complex<double> result{static_cast<double>(extended_result.real()),
                                      static_cast<double>(extended_result.imag())};
    const auto real_error = static_cast<double>(std::abs(extended_error.real()));
    const auto imaginary_error = static_cast<double>(std::abs(extended_error.imag()));
    const auto requested_real_error = std::max(epsabs, epsrel * std::abs(result.real()));
    const auto requested_imaginary_error = std::max(epsabs, epsrel * std::abs(result.imag()));
    const auto failed = !std::isfinite(result.real()) || !std::isfinite(result.imag())
                        || !std::isfinite(real_error) || !std::isfinite(imaginary_error)
                        || real_error > requested_real_error || imaginary_error > requested_imaginary_error;
    if (failed && integration.policy() != NRG::Tools::GslErrorPolicy::ignore) {
      std::ostringstream message;
      message << std::setprecision(OUTPUT_PRECISION)
              << "qag error: summed component error estimates (real=" << real_error
              << ", imaginary=" << imaginary_error << ") exceed requested global tolerances (real="
              << requested_real_error << ", imaginary=" << requested_imaginary_error << ')';
      if (failed_calls != 0)
        message << "; " << failed_calls << " local call(s) failed; first interval [" << first_lower << ", "
                << first_upper << "]: " << first_status << " -- " << gsl_strerror(first_status);
      if (integration.policy() == NRG::Tools::GslErrorPolicy::fail) throw std::runtime_error(message.str());
      integration.record_global_warning(message.str());
    }
    return result;
  }

  auto hilbert(const double x, const double y, integrator *integration) const {
    zero_density zero_fnc;
    const auto z = std::complex(x + shiftx, y + shifty); // shift here!
    if (algorithm == Algorithm::analytic) {
      if (!weighted_tabulated_rho) throw std::logic_error("Tabulated DOS polynomial is not initialized.");
      validate_hilbert_transform_argument(B, z, n);
      return NRG::Tools::cauchy_transform(*weighted_tabulated_rho, z);
    }
    if (!integration) throw std::logic_error("GSL integration workspace is not initialized.");
    if (tabulated) {
      return tabulated_qag(*integration, z);
    }
    auto density = [this](const auto energy) {
      return std::abs(energy * scale) < 1.0
               ? 2.0 / M_PI * scale * std::sqrt(1.0 - sqr(energy * scale))
               : 0.0;
    };
    return hilbert_transform_without_handler_guard(*integration, density, zero_fnc, B, z, lim_direct, n,
                                                   epsabs, epsrel);
  }

  auto calculate_points(const std::vector<transform_point> &points) {
    processed_points = points.size();
    actual_workers = std::min(jobs, points.size());
    std::vector<transform_outcome> outcomes(points.size());
    if (points.empty()) return std::vector<std::complex<double>>{};

    std::vector<integrator> integrations;
    if (algorithm == Algorithm::qag) {
      integrations.reserve(actual_workers);
      for (size_t worker = 0; worker < actual_workers; ++worker)
        integrations.emplace_back(integrator::configured, workspace_limit, quadrature_rule, gsl_error_policy);
    }

    std::atomic<size_t> next_index{0};
    std::atomic<bool> stopped{false};
    auto run_worker = [&](const size_t worker) {
      while (!stopped.load(std::memory_order_relaxed)) {
        const auto index = next_index.fetch_add(1, std::memory_order_relaxed);
        if (index >= points.size()) break;
        try {
          auto *integration = algorithm == Algorithm::qag ? &integrations[worker] : nullptr;
          if (integration) integration->set_failure_summary(&outcomes[index].failures);
          outcomes[index].value = hilbert(points[index].x, points[index].y, integration);
        } catch (...) {
          outcomes[index].failure = std::current_exception();
          stopped.store(true, std::memory_order_relaxed);
        }
      }
    };

    auto run_all_workers = [&] {
      if (actual_workers == 1) {
        run_worker(0);
        return;
      }
      std::vector<std::thread> workers;
      workers.reserve(actual_workers);
      try {
        for (size_t worker = 0; worker < actual_workers; ++worker) workers.emplace_back(run_worker, worker);
      } catch (...) {
        stopped.store(true, std::memory_order_relaxed);
        for (auto &worker : workers) worker.join();
        throw;
      }
      for (auto &worker : workers) worker.join();
    };

    if (algorithm == Algorithm::qag) {
      const NRG::Tools::GslErrorHandlerGuard error_handler;
      run_all_workers();
    } else {
      run_all_workers();
    }

    std::vector<std::complex<double>> results(points.size());
    for (size_t index = 0; index < outcomes.size(); ++index) {
      gsl_failures.merge(outcomes[index].failures);
      if (outcomes[index].failure) std::rethrow_exception(outcomes[index].failure);
      results[index] = outcomes[index].value;
    }
    return results;
  }

  void do_one(const double x, const double y, std::ostream &OUT) {
    const auto results = calculate_points({{x, y}});
    const auto res = results.front();
    if (!G)
      OUT << res.imag() << '\n';
    else
      OUT << res << '\n';
  }

  void do_stream(numeric_row_reader<3> &rows, std::ostream &OUT) {
    std::vector<std::array<double, 3>> input;
    while (const auto row = rows.next()) input.push_back(row->values);
    std::vector<transform_point> points;
    points.reserve(input.size());
    for (const auto &[label, x, y] : input) points.push_back({x, y});
    const auto results = calculate_points(points);
    for (size_t index = 0; index < input.size(); ++index) {
      const auto label = input[index][0];
      const auto res = results[index];
      if (!G)
        OUT << label << ' ' << res.imag() << '\n';
      else
        OUT << label << ' ' << res.real() << ' ' << res.imag() << '\n';
    }
  }

  void do_hilb(numeric_row_reader<2> &real_rows, numeric_row_reader<2> &imag_rows, std::ostream &Or, std::ostream &Oi) {
    size_t clipped = 0;
    double first_clipped_label = 0.0;
    double first_clipped_value = 0.0;
    auto report_clipping = [&]() {
      if (clipped == 0) return;
      const auto old_precision = std::cerr.precision();
      std::cerr << std::setprecision(OUTPUT_PRECISION)
                << "hilb: clipped " << clipped << " ImSigma value(s) to " << -clipping << "; first at omega=" << first_clipped_label
                << ": " << first_clipped_value << " -> " << -clipping << '\n';
      std::cerr.precision(old_precision);
    };
    std::vector<std::array<double, 2>> labels;
    std::vector<transform_point> points;
    try {
      while (true) {
        const auto real_row = real_rows.next();
        const auto imag_row = imag_rows.next();
        if (!real_row && !imag_row) break;
        if (!real_row)
          throw std::runtime_error(real_rows.source_name() + " ended after " + std::to_string(real_rows.records_read()) + " data rows; "
                                   + imag_rows.source_name() + ":" + std::to_string(imag_row->line_number) + " has an unpaired data row.");
        if (!imag_row)
          throw std::runtime_error(imag_rows.source_name() + " ended after " + std::to_string(imag_rows.records_read()) + " data rows; "
                                   + real_rows.source_name() + ":" + std::to_string(real_row->line_number) + " has an unpaired data row.");

        const auto [label1, real_sigma] = real_row->values;
        const auto [label2, imaginary_sigma] = imag_row->values;
        if (std::abs(label1 - label2) > frequency_tolerance)
          throw std::runtime_error("Frequency mismatch between " + real_rows.source_name() + ":" + std::to_string(real_row->line_number) + " ("
                                   + format_double(label1) + ") and " + imag_rows.source_name() + ":" + std::to_string(imag_row->line_number) + " ("
                                   + format_double(label2) + ").");

        auto y = imaginary_sigma;
        if (y > -clipping) {
          if (clipped == 0) {
            first_clipped_label = label2;
            first_clipped_value = y;
          }
          ++clipped;
          y = -clipping;
        }
        const auto x = label1 - real_sigma;
        y = -y;
        labels.push_back({label1, label2});
        points.push_back({x, y});
      }
      const auto results = calculate_points(points);
      for (size_t index = 0; index < results.size(); ++index) {
        const auto [label1, label2] = labels[index];
        const auto res = results[index];
        double resre = res.real();
        double resim = res.imag();
        if (!G) {
          resre /= -M_PI;
          resim /= -M_PI;
        }
        Or << label1 << ' ' << resre << '\n';
        Oi << label2 << ' ' << resim << '\n';
      }
    } catch (...) {
      report_clipping();
      throw;
    }
    report_clipping();
  }

  auto safe_open_rd(const std::string &filename) {
    std::ifstream F(filename);
    if (!F) throw std::runtime_error("Error opening file " + filename + " for reading.");
    return F;
  }

  static auto is_regular_output(const std::filesystem::path &filename) {
    if (filename == "/dev/stdout") return false;
    std::error_code error;
    const auto regular = std::filesystem::is_regular_file(filename, error);
    return !error && regular;
  }

  void load_dos(const std::string &filename) {
    if (verbosity >= 2) std::cerr << "Density of states filename: " << filename << '\n';
    auto F = safe_open_rd(filename);
    numeric_row_reader<2> rows(F, filename);
    std::vector<std::pair<double, double>> pts;
    while (const auto row = rows.next()) pts.emplace_back(row->values[0], row->values[1]);
    std::sort(pts.begin(), pts.end());
    Xpts.clear();
    Ypts.clear();
    for (const auto &[x, y] : pts) {
      Xpts.push_back(x);
      Ypts.push_back(y);
    }
  }

  void report_dos() {
    if (Xpts.empty()) throw std::runtime_error("DOS file contains no data points.");
    Xmin = Xpts.front();
    Xmax = Xpts.back();
    B = std::max(std::abs(Xmin), std::abs(Xmax));
    if (!std::isfinite(B) || B <= 0.0) throw std::runtime_error("DOS half-bandwidth must be finite and positive.");
    tabulated_rho.emplace(Xpts, Ypts, 0.0, interpolation_method);
    if (algorithm == Algorithm::analytic)
      weighted_tabulated_rho.emplace(
        tabulated_rho->piecewise_polynomial().multiply_by_monomial(static_cast<size_t>(n)));
    dos_integral = tabulated_rho->piecewise_polynomial().integral();
    if (!std::isfinite(*dos_integral)) throw std::runtime_error("Error: Integral is not a finite number.");
    if (verbosity >= 2) std::cerr << "Sum=" << *dos_integral << '\n';
  }

  void report_configuration(const int remaining, const std::vector<std::string> &args,
                            const std::optional<std::string> &output_filename,
                            const std::optional<std::string> &dos_filename,
                            const std::optional<std::pair<double, double>> &point) const {
    if (verbosity == 0) return;

    NRG::Tools::ConfigurationReport report("hilb");
    report.value("verbosity", verbosity);
    report.resolved("mode", remaining == 1 ? "argument-file" : remaining == 2 ? "single-point" : "dmft",
                    "positional argument count");
    report.value("output.mode", G ? "raw-complex" : remaining == 4 ? "minus-one-over-pi-scaled-components" : "imaginary");
    report.value("output.precision", OUTPUT_PRECISION);
    report.value("energy_power", n);
    report.value("argument.real_shift", shiftx);
    report.value("argument.imaginary_shift", shifty);
    report.value("algorithm", algorithm_name(algorithm));
    if (jobs_source == "--jobs")
      report.value("jobs", jobs);
    else
      report.resolved("jobs", jobs, jobs_source);
    report.value("jobs.source", jobs_source);

    if (tabulated) {
      report.value("dos.mode", "tabulated");
      report.value("dos.file", *dos_filename);
      report.value("dos.interpolation", NRG::Tools::interpolation_method_name(interpolation_method));
      report.value("dos.points", Xpts.size());
      report.resolved("dos.lower_bound", Xmin, "tabulated DOS");
      report.resolved("dos.upper_bound", Xmax, "tabulated DOS");
      report.resolved("dos.half_bandwidth", B, "tabulated DOS bounds");
      report.value("dos.integral", *dos_integral);
    } else {
      report.value("dos.mode", "semicircular");
      if (bandwidth_source == BandwidthSource::scale) {
        report.value("dos.scale", scale);
        report.resolved("dos.half_bandwidth", B, "reciprocal DOS scale");
      } else if (bandwidth_source == BandwidthSource::half_bandwidth) {
        report.value("dos.half_bandwidth", B);
        report.resolved("dos.scale", scale, "reciprocal half-bandwidth");
      } else {
        report.value("dos.scale", scale);
        report.value("dos.half_bandwidth", B);
      }
      report.resolved("dos.lower_bound", Xmin, "symmetric half-bandwidth");
      report.resolved("dos.upper_bound", Xmax, "symmetric half-bandwidth");
    }

    report.value("integration.qag.active", algorithm == Algorithm::qag);
    if (algorithm == Algorithm::qag) {
      report.value("integration.direct_threshold", lim_direct);
      report.value("integration.epsabs", epsabs);
      report.value("integration.epsrel", epsrel);
      report.value("integration.workspace_limit", workspace_limit);
      report.value("integration.quadrature_rule", static_cast<int>(quadrature_rule));
      report.value("integration.gsl_error_policy", NRG::Tools::gsl_error_policy_name(gsl_error_policy));
    }

    if (remaining == 1) {
      report.value("input.file", args[0]);
      if (output_filename)
        report.value("output.file", *output_filename);
      else
        report.resolved("output.file", "<stdout>", "-o not specified");
    } else if (remaining == 2) {
      report.value("argument.real", point->first);
      report.value("argument.imaginary", point->second);
      if (output_filename)
        report.value("output.file", *output_filename);
      else
        report.resolved("output.file", "<stdout>", "-o not specified");
    } else {
      report.value("input.real_self_energy", args[0]);
      report.value("input.imaginary_self_energy", args[1]);
      report.value("output.real", args[2]);
      report.value("output.imaginary", args[3]);
      report.value("dmft.clipping", clipping);
      report.value("dmft.frequency_tolerance", frequency_tolerance);
    }
    report.write(std::cerr);
  }

  void about() {
    std::cout << "# hilb -- Hilbert transformer for arbitrary density of states.\n";
  }

  void report_timing() const {
    const auto wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    const auto cpu_end = std::clock();
    std::optional<double> cpu_seconds;
    if (cpu_start != static_cast<std::clock_t>(-1) && cpu_end != static_cast<std::clock_t>(-1)
        && cpu_end >= cpu_start)
      cpu_seconds = static_cast<double>(cpu_end - cpu_start) / static_cast<double>(CLOCKS_PER_SEC);

    std::ostringstream report;
    report << std::setprecision(6) << "Time elapsed: " << wall_seconds << " s\n";
    if (verbosity != 0) {
      const auto throughput = wall_seconds > 0.0 ? static_cast<double>(processed_points) / wall_seconds : 0.0;
      report << "Performance: wall=" << wall_seconds << "s";
      if (cpu_seconds && wall_seconds > 0.0)
        report << " cpu=" << *cpu_seconds << "s effective_parallelism=" << *cpu_seconds / wall_seconds << "x";
      else
        report << " cpu=n/a effective_parallelism=n/a";
      report << " throughput=" << throughput << " points/s workers=" << actual_workers
             << " algorithm=" << algorithm_name(algorithm) << '\n';
    }
    std::cout << report.str();
    NRG::Tools::finish_output(std::cout, "<stdout>");
  }

  void usage() {
    std::cout << "Usage (1): hilb [options] <x> <y>\n";
    std::cout << "Usage (2): hilb [options] <inputfile>\n";
    std::cout << "Usage (3): hilb [options] <resigma.dat> <imsigma.dat> <reaw.dat> <imaw.dat>\n\n";
    std::cout << "Options:\n";
    std::cout << "-h, --help  show help\n";
    std::cout << "-v        Show resolved configuration on standard error\n";
    std::cout << "-vv       Also show detailed setup diagnostics\n";
    std::cout << "-V, --version\n";
    std::cout << "          Show project version\n";
    std::cout << "-d <dos>  Load the density of state data from file 'dos'\n";
    std::cout << "          If this option is not used, the Bethe lattice DOS is assumed.\n";
    std::cout << "--algorithm <algorithm>\n";
    std::cout << "          qag or analytic. Default is qag; analytic requires -d.\n";
    std::cout << "-i <method>, --interpolation <method>\n";
    std::cout << "          Interpolation method: linear, cspline, akima, or steffen. Default is cspline.\n";
    std::cout << "-j <n>, --jobs <n>\n";
    std::cout << "          Worker count. Default is the first OMP_NUM_THREADS value, or 1.\n";
    std::cout << "-s <s>    Rescale factor 'scale' for the DOS.\n";
    std::cout << "-B <B>    Half-bandwidth 'B' of the Bethe lattice DOS.\n";
    std::cout << "          Use either -s or -B. Default is scale=B=1.\n";
    std::cout << "-n <n>    Multiply the integrand by E^n; n must be a nonnegative integer.\n";
    std::cout << "          Default is n=0.\n";
    std::cout << "-G        Return both Re[H_n(z)] and Im[H_n(z)].\n";
    std::cout << "-o <file> Write output to file in single-point and input-file modes.\n";
    std::cout << "-x <dx>   Shift real part of the argument\n";
    std::cout << "-y <dy>   Shift imag part of the argument\n";
    std::cout << "-c <c>    Minimum magnitude used to clip ImSigma in DMFT mode.\n";
    std::cout << "          Default is sqrt(numeric_limits<double>::min()).\n";
    std::cout << "-t <t>    QAG direct-integration threshold. Default is 1e-3.\n";
    std::cout << "-a <a>, --epsabs <a>\n";
    std::cout << "          QAG absolute integration tolerance. Default is 1e-14.\n";
    std::cout << "-r <r>, --epsrel <r>\n";
    std::cout << "          QAG relative integration tolerance. Default is 1e-10.\n";
    std::cout << "--workspace-limit <n>\n";
    std::cout << "          QAG integration workspace size per worker. Default is 1000.\n";
    std::cout << "--quadrature-rule <rule>\n";
    std::cout << "          QAG Gauss-Kronrod rule: 15, 21, 31, 41, 51, or 61. Default is 15.\n";
    std::cout << "--gsl-error-policy <policy>\n";
    std::cout << "          QAG GSL error policy: ignore, warn, or fail. Default is warn.\n";
    std::cout << "          QAG numerical controls are ignored by --algorithm analytic.\n";
    std::cout << "-f <f>    Frequency-label tolerance in DMFT mode. Default is 1e-6.\n";
  }

  void parse_param_run(int argc, char *argv[]) {
    enum : int {
      workspace_limit_option = 256,
      quadrature_rule_option,
      gsl_error_policy_option,
      algorithm_option
    };
    static const option long_options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"algorithm", required_argument, nullptr, algorithm_option},
      {"interpolation", required_argument, nullptr, 'i'},
      {"jobs", required_argument, nullptr, 'j'},
      {"epsabs", required_argument, nullptr, 'a'},
      {"epsrel", required_argument, nullptr, 'r'},
      {"workspace-limit", required_argument, nullptr, workspace_limit_option},
      {"quadrature-rule", required_argument, nullptr, quadrature_rule_option},
      {"gsl-error-policy", required_argument, nullptr, gsl_error_policy_option},
      {nullptr, 0, nullptr, 0}
    };
    std::optional<std::string> output_filename;
    std::optional<std::string> dos_filename;
    std::optional<size_t> requested_jobs;
    int c;
    while ((c = getopt_long(argc, argv, "hGd:i:j:vs:B:o:x:y:n:c:t:a:r:f:", long_options, nullptr)) != -1) {
      switch (c) {
        case 'h': usage(); return;
        case 'G': G = true; break;
        case 'd':
          tabulated = true;
          dos_filename = optarg;
          break;
        case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
        case 'j': requested_jobs = NRG::Tools::parse_positive_size(optarg, "Worker count"); break;
        case 'v': ++verbosity; break;
        case 's':
          scale = parse_finite_double(optarg, "DOS scale");
          if (scale <= 0.0) throw std::runtime_error("DOS scale must be positive.");
          B = 1.0 / scale;
          if (!std::isfinite(B)) throw std::runtime_error("DOS scale reciprocal must be finite.");
          Xmin = -B;
          Xmax = +B;
          bandwidth_source = BandwidthSource::scale;
          break;
        case 'B':
          B = parse_finite_double(optarg, "Half-bandwidth");
          if (B <= 0.0) throw std::runtime_error("Half-bandwidth must be positive.");
          scale = 1.0 / B;
          if (!std::isfinite(scale)) throw std::runtime_error("Half-bandwidth reciprocal must be finite.");
          Xmin = -B;
          Xmax = +B;
          bandwidth_source = BandwidthSource::half_bandwidth;
          break;
        case 'o':
          output_filename = optarg;
          break;
        case 'x':
          shiftx = parse_finite_double(optarg, "Real shift");
          break;
        case 'y':
          shifty = parse_finite_double(optarg, "Imaginary shift");
          break;
        case 'n':
          n = parse_energy_power(optarg);
          break;
        case 'c':
          clipping = parse_finite_double(optarg, "Clipping magnitude");
          if (clipping < minimum_safe_imaginary_part())
            throw std::runtime_error("Clipping magnitude must be at least sqrt(numeric_limits<double>::min()).");
          break;
        case 't':
          lim_direct = parse_finite_double(optarg, "Direct-integration threshold");
          if (lim_direct < 0.0) throw std::runtime_error("Direct-integration threshold must be nonnegative.");
          break;
        case 'a':
          epsabs = parse_finite_double(optarg, "Absolute integration tolerance");
          if (epsabs < 0.0) throw std::runtime_error("Absolute integration tolerance must be nonnegative.");
          break;
        case 'r':
          epsrel = parse_finite_double(optarg, "Relative integration tolerance");
          if (epsrel < 0.0) throw std::runtime_error("Relative integration tolerance must be nonnegative.");
          break;
        case 'f':
          frequency_tolerance = parse_finite_double(optarg, "Frequency-label tolerance");
          if (frequency_tolerance < 0.0) throw std::runtime_error("Frequency-label tolerance must be nonnegative.");
          break;
        case workspace_limit_option:
          workspace_limit = NRG::Tools::parse_positive_size(optarg, "Workspace limit");
          break;
        case quadrature_rule_option:
          quadrature_rule = NRG::Tools::parse_qag_rule(optarg);
          break;
        case gsl_error_policy_option:
          gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg);
          break;
        case algorithm_option: algorithm = parse_algorithm(optarg); break;
        default: throw std::runtime_error("Invalid command-line option.");
      }
    }
    if (algorithm == Algorithm::qag) {
      try {
        validate_integration_settings(lim_direct, epsabs, epsrel);
        NRG::Tools::validate_qag_workspace_limit(workspace_limit);
      } catch (const std::invalid_argument &error) {
        throw std::runtime_error(error.what());
      }
    }
    resolve_jobs(requested_jobs);
    const auto remaining = argc - optind; // arguments left
    const std::vector<std::string> args(argv+optind, argv+argc); // NOLINT
    if (remaining != 1 && remaining != 2 && remaining != 4)
      throw std::runtime_error("Invalid number of positional arguments: got " + std::to_string(remaining) + "; expected 1, 2, or 4.");
    if (algorithm == Algorithm::analytic && !tabulated)
      throw std::runtime_error("The analytic algorithm requires a tabulated density selected with -d.");
    if (remaining == 4 && output_filename) throw std::runtime_error("The -o option is not valid in DMFT mode.");
    if (remaining == 1 && output_filename && files_refer_to_same_location(*output_filename, args[0]))
      throw std::runtime_error("Input and output files must be different.");
    if (remaining == 4
        && (files_refer_to_same_location(args[2], args[3]) || files_refer_to_same_location(args[2], args[0])
            || files_refer_to_same_location(args[2], args[1]) || files_refer_to_same_location(args[3], args[0])
            || files_refer_to_same_location(args[3], args[1])))
      throw std::runtime_error("DMFT output files must differ from each other and from both input files.");
    if (tabulated && output_filename && files_refer_to_same_location(*dos_filename, *output_filename))
      throw std::runtime_error("Density-of-states input and output files must be different.");
    if (tabulated && remaining == 4
        && (files_refer_to_same_location(*dos_filename, args[2]) || files_refer_to_same_location(*dos_filename, args[3])))
      throw std::runtime_error("DMFT output files must differ from the density-of-states input file.");
    std::optional<std::pair<double, double>> point;
    if (remaining == 2) {
      const auto x = parse_finite_double(args[0], "Real argument");
      const auto y = parse_finite_double(args[1], "Imaginary argument");
      point.emplace(x, y);
    }
    if (tabulated) {
      load_dos(*dos_filename);
      report_dos();
    }

    report_configuration(remaining, args, output_filename, dos_filename, point);

    if (remaining == 1) {
      about();
      auto F = safe_open_rd(args[0]);
      numeric_row_reader<3> rows(F, args[0]);
      if (output_filename) {
        if (verbosity >= 2) std::cerr << "Output file: " << *output_filename << '\n';
        std::ostringstream output;
        output << std::setprecision(OUTPUT_PRECISION);
        do_stream(rows, output);
        NRG::Tools::write_output_file(*output_filename, output.str());
        if (is_regular_output(*output_filename)) report_timing();
      } else {
        do_stream(rows, std::cout);
        finish_output(std::cout, "<stdout>");
      }
      return;
    }
    if (remaining == 2) {
      if (output_filename) {
        if (verbosity >= 2) std::cerr << "Output file: " << *output_filename << '\n';
        std::ostringstream output;
        output << std::setprecision(OUTPUT_PRECISION);
        do_one(point->first, point->second, output);
        NRG::Tools::write_output_file(*output_filename, output.str());
        if (is_regular_output(*output_filename)) report_timing();
      } else {
        do_one(point->first, point->second, std::cout);
        finish_output(std::cout, "<stdout>");
      }
      return;
    }

    about();
    auto Frs = safe_open_rd(args[0]);
    auto Fis = safe_open_rd(args[1]);
    numeric_row_reader<2> real_rows(Frs, args[0]);
    numeric_row_reader<2> imag_rows(Fis, args[1]);
    std::ostringstream real_output;
    std::ostringstream imaginary_output;
    real_output << std::setprecision(OUTPUT_PRECISION);
    imaginary_output << std::setprecision(OUTPUT_PRECISION);
    do_hilb(real_rows, imag_rows, real_output, imaginary_output);
    NRG::Tools::write_output_file(args[2], real_output.str());
    NRG::Tools::write_output_file(args[3], imaginary_output.str());
    if (is_regular_output(args[2]) && is_regular_output(args[3])) report_timing();
  }

  public:
  Hilb(int argc, char *argv[])
      : wall_start{std::chrono::steady_clock::now()}, cpu_start{std::clock()} {
    std::cout << std::setprecision(OUTPUT_PRECISION);
    try {
      parse_param_run(argc, argv);
      finish_output(std::cout, "<stdout>");
    } catch (...) {
      gsl_failures.report(std::cerr);
      throw;
    }
    gsl_failures.report(std::cerr);
  }
};

} // namespace

#endif
