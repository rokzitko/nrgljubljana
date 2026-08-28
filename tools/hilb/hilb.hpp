// Computes \int E^n rho(E)/(z-E) dE for given z and tabulated density of states rho(E).
//
// Legacy modes:
// Mode 1: <Re z> <Im z> as input. Returns Im part by default, or both Re and Im parts if the -G switch is used.
// Mode 2: read x and y from a file.
// Mode 3: convert imsigma/resigma.dat to imaw/reaw.dat files in the DMFT loop.
//
// Rok Zitko, rok.zitko@ijs.si, 2009-2026

#ifndef _hilb_hilb_hpp_
#define _hilb_hilb_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstddef>
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
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <charconv>
#include <cerrno>
#include <cfloat>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <getopt.h>
#include <unistd.h>

#include <gsl/gsl_errno.h> // GNU scientific library
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "../common/gsl_config.hpp"

namespace NRG::Hilb {

using std::size_t;

inline constexpr std::string_view HILB_VERSION = "2026.09";
inline constexpr int OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;

inline auto format_double(const double value) {
  std::ostringstream output;
  output << std::setprecision(OUTPUT_PRECISION) << value;
  return output.str();
}

struct gsl_accel_deleter {
  void operator()(gsl_interp_accel *acc) const {
    if (acc) gsl_interp_accel_free(acc);
  }
};

struct gsl_spline_deleter {
  void operator()(gsl_spline *spline) const {
    if (spline) gsl_spline_free(spline);
  }
};

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

// Unwrap a callable and evaluate it at x.
template <typename F>
inline double unwrap(const double x, void *p) {
  auto *f = static_cast<F *>(p);
  return (*f)(x);
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

class gsl_failure_summary {
  private:
  size_t total = 0;
  size_t nonfinite = 0;
  std::map<int, size_t> by_status;
  std::optional<gsl_failure> first;
  bool reported = false;

  public:
  void record(const gsl_failure &failure) {
    ++total;
    if (failure.status != GSL_SUCCESS) ++by_status[failure.status];
    if (failure.nonfinite) ++nonfinite;
    if (!first) first = failure;
  }

  auto count() const noexcept { return total; }

  void report(std::ostream &out) {
    if (reported || total == 0) return;
    reported = true;
    const auto old_precision = out.precision();
    out << std::setprecision(OUTPUT_PRECISION) << "hilb: warning: " << total << " GSL integration call(s) reported failure";
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
    if (!work) throw std::runtime_error("Cannot use a moved-from GSL integrator.");
    const NRG::Tools::GslErrorHandlerGuard error_handler;
    gsl_function function{&unwrap<F>, &f};
    double result = std::numeric_limits<double>::quiet_NaN();
    double error = std::numeric_limits<double>::quiet_NaN();
    const auto status = gsl_integration_qag(&function, a, b, epsabs, epsrel, limit, NRG::Tools::gsl_qag_rule(rule), work.get(), &result, &error);
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
};

// Wrap around GSL interpolation routines
class interpolator {
  private:
  size_t len = 0;                  // number of data points
  std::vector<double> X, Y;        // X and Y tables
  double Xmin = 0.0, Xmax = 0.0;   // boundary points
  NRG::Tools::InterpolationMethod method;
  double oob_value;                // out-of-boundary value
  std::unique_ptr<gsl_interp_accel, gsl_accel_deleter> acc; // workspace
  std::unique_ptr<gsl_spline, gsl_spline_deleter> spline;   // spline data

  void initialize() {
    if (X.size() != Y.size()) throw std::invalid_argument("Interpolation grids must have equal sizes.");
    const auto minimum_size = NRG::Tools::interpolation_minimum_size(method);
    if (X.size() < minimum_size) {
      if (method == NRG::Tools::InterpolationMethod::cspline)
        throw std::invalid_argument("Cubic interpolation requires at least " + std::to_string(minimum_size) + " points.");
      throw std::invalid_argument(std::string(NRG::Tools::interpolation_method_name(method)) + " interpolation requires at least "
                                  + std::to_string(minimum_size) + " points.");
    }
    for (size_t i = 0; i < X.size(); ++i) {
      if (!std::isfinite(X[i]) || !std::isfinite(Y[i])) throw std::invalid_argument("Interpolation data must be finite.");
      if (i != 0 && !(X[i - 1] < X[i])) throw std::invalid_argument("Interpolation energies must be strictly increasing.");
    }

    acc.reset(gsl_interp_accel_alloc());
    if (!acc) throw std::runtime_error("Failed to allocate GSL interpolation accelerator.");
    len = X.size();
    spline.reset(gsl_spline_alloc(NRG::Tools::gsl_interpolation_type(method), len));
    if (!spline) {
      if (method == NRG::Tools::InterpolationMethod::cspline) throw std::runtime_error("Failed to allocate GSL cubic spline.");
      throw std::runtime_error("Failed to allocate GSL interpolation spline.");
    }
    const auto status = gsl_spline_init(spline.get(), X.data(), Y.data(), len);
    if (status != GSL_SUCCESS) {
      if (method == NRG::Tools::InterpolationMethod::cspline)
        throw std::runtime_error("Failed to initialize GSL cubic spline: " + std::string(gsl_strerror(status)));
      throw std::runtime_error("Failed to initialize GSL interpolation spline: " + std::string(gsl_strerror(status)));
    }
    Xmin = X.front();
    Xmax = X.back();
  }

  public:
  interpolator(const std::vector<double> &_X, const std::vector<double> &_Y, const double _oob_value = 0.0)
      : X{_X}, Y{_Y}, method{NRG::Tools::InterpolationMethod::cspline}, oob_value{_oob_value} {
    initialize();
  }
  interpolator(const std::vector<double> &_X, const std::vector<double> &_Y, const double _oob_value,
               const NRG::Tools::InterpolationMethod method_)
      : X{_X}, Y{_Y}, method{method_}, oob_value{_oob_value} {
    initialize();
  }
  interpolator(const interpolator &I) : X{I.X}, Y{I.Y}, method{I.method}, oob_value{I.oob_value} {
    initialize();
  }
  interpolator(interpolator &&I) = default;
  interpolator &operator=(const interpolator &I) {
    if (this == &I) return *this;
    interpolator copy(I);
    *this = std::move(copy);
    return *this;
  }
  interpolator &operator=(interpolator &&I) = default;
  auto operator()(const double x) { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(spline.get(), x, acc.get()) : oob_value); }
};

// Square of x
inline auto sqr(const double x) { return x * x; }

// Result of Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}] (atg -> imQ).
inline auto imQ(const double x, const double y, const double B) { return std::atan((-B + x) / y) - std::atan((B + x) / y); }

// Result of Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}] (logs -> reQ).
inline auto reQ(const double x, const double y, const double B) { return (-std::log(sqr(B - x) + sqr(y)) + std::log(sqr(B + x) + sqr(y))) / 2.0; }

// Calculate the (half)bandwidth, i.e., the size B of the enclosing interval [-B:B].
inline auto bandwidth(const std::vector<double> &X) {
  assert(std::is_sorted(X.begin(), X.end()));
  const auto Xmin = X.front();
  const auto Xmax = X.back();
  return std::max(std::abs(Xmin), std::abs(Xmax));
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
auto hilbert_transform(integrator &integration, FNCR &rhor, FNCI &rhoi, const double B, const std::complex<double> z,
                        const double lim_direct = 1e-3, const int n = 0, const double epsabs = 1e-14, const double epsrel = 1e-10) {
  gsl_set_error_handler_off();
  if (n < 0) throw std::invalid_argument("Energy power must be nonnegative.");
  validate_integration_settings(lim_direct, epsabs, epsrel);
  if (!std::isfinite(B) || B <= 0.0) throw std::invalid_argument("Half-bandwidth must be finite and positive.");
  const auto x = real(z);
  const auto y = imag(z);
  if (!std::isfinite(x) || !std::isfinite(y)) throw std::invalid_argument("Hilbert-transform argument must be finite.");
  if (std::abs(y) < minimum_safe_imaginary_part())
    throw std::invalid_argument("Absolute imaginary part is below the minimum safe value sqrt(numeric_limits<double>::min()).");
  auto g_r = [&rhor, n](const double omega) -> double {
    return n == 0 ? rhor(omega) : std::pow(omega, n) * rhor(omega);
  };
  auto g_i = [&rhoi, n](const double omega) -> double {
    return n == 0 ? rhoi(omega) : std::pow(omega, n) * rhoi(omega);
  };

  if (std::abs(y) >= lim_direct) {
    auto ref0 = [x, y, &g_r, &g_i](double omega) -> double {
      return (g_r(omega) * (x - omega) + g_i(omega) * y) / (sqr(y) + sqr(x - omega));
    };
    auto imf0 = [x, y, &g_r, &g_i](double omega) -> double {
      return (g_r(omega) * (-y) + g_i(omega) * (x - omega)) / (sqr(y) + sqr(x - omega));
    };
    const auto real_result = integration(ref0, -B, B, epsabs, epsrel);
    const auto imag_result = integration(imf0, -B, B, epsabs, epsrel);
    return std::complex(real_result, imag_result);
  }

  const bool inside = -B <= x && x <= B;
  const double g_r_x = inside ? g_r(x) : 0.0;
  const double g_i_x = inside ? g_i(x) : 0.0;

  // Integrate the singularity-subtracted remainder after the logarithmic change of variables.
  auto calc_subtracted = [&integration, x, y, B, inside, epsabs, epsrel](auto f3p, auto f3m, auto analytic) -> double {
    const double W1 = (x - B) / std::abs(y); // Rescaled integration limits. Only the absolute value of y matters here.
    const double W2 = (B + x) / std::abs(y);
    assert(W2 >= W1);
    double lim1down = 1.0, lim1up = -1.0, lim2down = 1.0, lim2up = -1.0;
    if (inside) {
      constexpr double log_zero_cutoff = -36.8; // approximately log(10^-16)
      if (W1 < 0.0) {
        lim1down = log_zero_cutoff;
        lim1up = std::log(-W1);
      }
      if (W2 > 0.0) {
        lim2down = log_zero_cutoff;
        lim2up = std::log(W2);
      }
    } else if (W1 > 0 && W2 > 0) { // x above the band
      lim2down = std::log(W1);
      lim2up = std::log(W2);
    } else if (W1 < 0 && W2 < 0) { // x below the band
      lim1down = std::log(-W2);
      lim1up = std::log(-W1);
    }
    const auto result1 = lim1down < lim1up ? integration(f3p, lim1down, lim1up, epsabs, epsrel) : 0.0;
    const auto result2 = lim2down < lim2up ? integration(f3m, lim2down, lim2up, epsabs, epsrel) : 0.0;
    const auto result3 = inside ? analytic : 0.0;
    return result1 + result2 + result3;
  };

  // Re part of g(omega)/(z-omega) with the singularity subtracted out.
  auto ref1 = [x, y, g_r_x, g_i_x, &g_r, &g_i](double omega) -> double {
    return ((g_r(omega) - g_r_x) * (x - omega) + (g_i(omega) - g_i_x) * y) / (sqr(y) + sqr(x - omega));
  };
  auto ref2  = [x, y, ref1](double W) -> double { return std::abs(y) * ref1(std::abs(y) * W + x); };
  auto ref3p = [ref2](double r) -> double { return ref2(std::exp(r)) * std::exp(r); };
  auto ref3m = [ref2](double r) -> double { return ref2(-std::exp(r)) * std::exp(r); };
  const auto red = g_r_x * reQ(x, y, B) - g_i_x * imQ(x, y, B);

  // Im part of g(omega)/(z-omega) with the singularity subtracted out.
  auto imf1 = [x, y, g_r_x, g_i_x, &g_r, &g_i](double omega) -> double {
    return ((g_r(omega) - g_r_x) * (-y) + (g_i(omega) - g_i_x) * (x - omega)) / (sqr(y) + sqr(x - omega));
  };
  auto imf2  = [x, y, imf1](double W) -> double { return std::abs(y) * imf1(std::abs(y) * W + x); };
  auto imf3p = [imf2](double r) -> double { return imf2(std::exp(r)) * std::exp(r); };
  auto imf3m = [imf2](double r) -> double { return imf2(-std::exp(r)) * std::exp(r); };
  const auto imd = g_r_x * imQ(x, y, B) + g_i_x * reQ(x, y, B);

  const auto real_result = calc_subtracted(ref3p, ref3m, red);
  const auto imag_result = calc_subtracted(imf3p, imf3m, imd);
  return std::complex(real_result, imag_result);
}

template <typename FNCR, typename FNCI>
auto hilbert_transform(FNCR rhor, FNCI rhoi, const double B, const std::complex<double> z, const double lim_direct = 1e-3,
                        const int n = 0, const double epsabs = 1e-14, const double epsrel = 1e-10) {
  gsl_set_error_handler_off();
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
  gsl_set_error_handler_off();
  interpolator rhor(Xpts, Rpts, 0.0, interpolation_method);
  interpolator rhoi(Xpts, Ipts, 0.0, interpolation_method);
  const double B = bandwidth(Xpts);
  integrator integration;
  return hilbert_transform(integration, rhor, rhoi, B, z, lim_direct, n, epsabs, epsrel);
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
  double scale = 1.0;         // scale factor
  double B     = 1.0 / scale; // half-bandwidth
  double Xmin = -B;
  double Xmax = +B;
  bool verbose     = false;
  bool G           = false; // G(z). Reports real and imaginary part.
  std::vector<double> Xpts, Ypts;
  std::optional<interpolator> tabulated_rho;
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
  size_t workspace_limit = 1000;
  NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::warn;
  gsl_failure_summary gsl_failures;
  std::optional<integrator> integration;

  auto hilbert(const double x, const double y) {
    auto Bethe_fnc = [this](const auto w) { return std::abs(w*scale) < 1.0 ? 2.0 / M_PI * scale * std::sqrt(1 - sqr(w * scale)) : 0.0; };
    auto zero_fnc = []([[maybe_unused]] const auto w) { return 0.0; };
    const auto z = std::complex(x + shiftx, y + shifty); // shift here!
    if (tabulated) {
      if (!tabulated_rho) throw std::logic_error("Tabulated DOS spline is not initialized.");
      return hilbert_transform(*integration, *tabulated_rho, zero_fnc, B, z, lim_direct, n, epsabs, epsrel);
    }
    return hilbert_transform(*integration, Bethe_fnc, zero_fnc, B, z, lim_direct, n, epsabs, epsrel);
  }

  void do_one(const double x, const double y, std::ostream &OUT) {
    const auto res = hilbert(x, y);
    if (!G)
      OUT << res.imag() << '\n';
    else
      OUT << res << '\n';
  }

  void do_stream(numeric_row_reader<3> &rows, std::ostream &OUT) {
    while (const auto row = rows.next()) {
      const auto [label, x, y] = row->values;
      const auto res = hilbert(x, y);
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
        const auto res = hilbert(x, y);
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

  auto safe_open_wr(const std::string &filename) {
    std::ofstream F(filename);
    if (!F) throw std::runtime_error("Error opening file " + filename + " for writing.");
    F << std::setprecision(OUTPUT_PRECISION);
    return F;
  }

  void load_dos(const std::string &filename) {
    if (verbose) std::cout << "Density of states filename: " << filename << '\n';
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
    const auto sum = (*integration)(std::ref(*tabulated_rho), -B, B, epsabs, epsrel);
    if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
    if (verbose) std::cout << "Sum=" << sum << '\n';
  }

  void info() {
    if (tabulated)
      std::cout << "Xmin=" << Xmin << " Xmax=" << Xmax << '\n';
    else
      std::cout << "Semicircular DOS. scale=" << scale << '\n';
    std::cout << "B=" << B << '\n';
    if (shiftx != 0.0) std::cout << "shiftx=" << shiftx << '\n';
    if (shifty != 0.0) std::cout << "shifty=" << shifty << '\n';
    if (n != 0) std::cout << "n=" << n << '\n';
  }

  void about() {
    std::cout << "# hilb -- Hilbert transformer for arbitrary density of states.\n";
    std::cout << "# Rok Zitko, rok.zitko@ijs.si, 2009-2026\n";
  }

  void version() {
    std::cout << "hilb " << HILB_VERSION << '\n';
  }

  void usage() {
    std::cout << "Usage (1): hilb [options] <x> <y>\n";
    std::cout << "Usage (2): hilb [options] <inputfile>\n";
    std::cout << "Usage (3): hilb [options] <resigma.dat> <imsigma.dat> <reaw.dat> <imaw.dat>\n\n";
    std::cout << "Options:\n";
    std::cout << "-h, --help  show help\n";
    std::cout << "-V        show version\n";
    std::cout << "-d <dos>  Load the density of state data from file 'dos'\n";
    std::cout << "          If this option is not used, the Bethe lattice DOS is assumed.\n";
    std::cout << "-i <method>, --interpolation <method>\n";
    std::cout << "          Interpolation method: linear, cspline, or akima. Default is cspline.\n";
    std::cout << "-v        Increase verbosity\n";
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
    std::cout << "-t <t>    Direct-integration threshold. Default is 1e-3.\n";
    std::cout << "-a <a>, --epsabs <a>\n";
    std::cout << "          Absolute integration tolerance. Default is 1e-14.\n";
    std::cout << "-r <r>, --epsrel <r>\n";
    std::cout << "          Relative integration tolerance. Default is 1e-10.\n";
    std::cout << "--workspace-limit <n>\n";
    std::cout << "          GSL integration workspace size. Default is 1000.\n";
    std::cout << "--quadrature-rule <rule>\n";
    std::cout << "          Gauss-Kronrod rule: 15, 21, 31, 41, 51, or 61. Default is 15.\n";
    std::cout << "--gsl-error-policy <policy>\n";
    std::cout << "          GSL error policy: ignore, warn, or fail. Default is warn.\n";
    std::cout << "-f <f>    Frequency-label tolerance in DMFT mode. Default is 1e-6.\n";
  }

  void parse_param_run(int argc, char *argv[]) {
    enum : int {
      workspace_limit_option = 256,
      quadrature_rule_option,
      gsl_error_policy_option
    };
    static const option long_options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"interpolation", required_argument, nullptr, 'i'},
      {"epsabs", required_argument, nullptr, 'a'},
      {"epsrel", required_argument, nullptr, 'r'},
      {"workspace-limit", required_argument, nullptr, workspace_limit_option},
      {"quadrature-rule", required_argument, nullptr, quadrature_rule_option},
      {"gsl-error-policy", required_argument, nullptr, gsl_error_policy_option},
      {nullptr, 0, nullptr, 0}
    };
    std::optional<std::string> output_filename;
    std::optional<std::string> dos_filename;
    int c;
    while ((c = getopt_long(argc, argv, "hGd:i:vVs:B:o:x:y:n:c:t:a:r:f:", long_options, nullptr)) != -1) {
      switch (c) {
        case 'h': usage(); return;
        case 'V': version(); return;
        case 'G': G = true; break;
        case 'd':
          tabulated = true;
          dos_filename = optarg;
          break;
        case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
        case 'v': verbose = true; break;
        case 's':
          scale = parse_finite_double(optarg, "DOS scale");
          if (scale <= 0.0) throw std::runtime_error("DOS scale must be positive.");
          B = 1.0 / scale;
          if (!std::isfinite(B)) throw std::runtime_error("DOS scale reciprocal must be finite.");
          Xmin = -B;
          Xmax = +B;
          break;
        case 'B':
          B = parse_finite_double(optarg, "Half-bandwidth");
          if (B <= 0.0) throw std::runtime_error("Half-bandwidth must be positive.");
          scale = 1.0 / B;
          if (!std::isfinite(scale)) throw std::runtime_error("Half-bandwidth reciprocal must be finite.");
          Xmin = -B;
          Xmax = +B;
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
        case workspace_limit_option: workspace_limit = NRG::Tools::parse_positive_size(optarg, "Workspace limit"); break;
        case quadrature_rule_option: quadrature_rule = NRG::Tools::parse_qag_rule(optarg); break;
        case gsl_error_policy_option: gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg); break;
        default: throw std::runtime_error("Invalid command-line option.");
      }
    }
    try {
      validate_integration_settings(lim_direct, epsabs, epsrel);
    } catch (const std::invalid_argument &error) {
      throw std::runtime_error(error.what());
    }
    const auto remaining = argc - optind; // arguments left
    const std::vector<std::string> args(argv+optind, argv+argc); // NOLINT
    if (remaining != 1 && remaining != 2 && remaining != 4)
      throw std::runtime_error("Invalid number of positional arguments: got " + std::to_string(remaining) + "; expected 1, 2, or 4.");
    if (remaining == 4 && output_filename) throw std::runtime_error("The -o option is not valid in DMFT mode.");
    if (remaining == 1 && output_filename && *output_filename == args[0])
      throw std::runtime_error("Input and output files must be different.");
    if (remaining == 4
        && (args[2] == args[3] || args[2] == args[0] || args[2] == args[1] || args[3] == args[0] || args[3] == args[1]))
      throw std::runtime_error("DMFT output files must differ from each other and from both input files.");
    if (tabulated && output_filename && *dos_filename == *output_filename)
      throw std::runtime_error("Density-of-states input and output files must be different.");
    if (tabulated && remaining == 4 && (*dos_filename == args[2] || *dos_filename == args[3]))
      throw std::runtime_error("DMFT output files must differ from the density-of-states input file.");
    integration.emplace(integrator::configured, workspace_limit, quadrature_rule, gsl_error_policy, &gsl_failures);
    if (tabulated) {
      load_dos(*dos_filename);
      report_dos();
    }

    if (remaining == 1) {
      about();
      if (verbose) info();
      auto F = safe_open_rd(args[0]);
      numeric_row_reader<3> rows(F, args[0]);
      std::optional<std::ofstream> output;
      if (output_filename) {
        output = safe_open_wr(*output_filename);
        if (verbose) std::cout << "Output file: " << *output_filename << '\n';
      }
      do_stream(rows, output ? *output : std::cout);
      return;
    }
    if (remaining == 2) {
      const auto x = parse_finite_double(args[0], "Real argument");
      const auto y = parse_finite_double(args[1], "Imaginary argument");
      std::optional<std::ofstream> output;
      if (output_filename) {
        output = safe_open_wr(*output_filename);
        if (verbose) std::cout << "Output file: " << *output_filename << '\n';
      }
      do_one(x, y, output ? *output : std::cout);
      return;
    }

    about();
    if (verbose) info();
    auto Frs = safe_open_rd(args[0]);
    auto Fis = safe_open_rd(args[1]);
    numeric_row_reader<2> real_rows(Frs, args[0]);
    numeric_row_reader<2> imag_rows(Fis, args[1]);
    auto Fra = safe_open_wr(args[2]);
    auto Fia = safe_open_wr(args[3]);
    do_hilb(real_rows, imag_rows, Fra, Fia);
  }

  public:
  Hilb(int argc, char *argv[]) {
    std::cout << std::setprecision(OUTPUT_PRECISION);
    gsl_set_error_handler_off();
    try {
      parse_param_run(argc, argv);
    } catch (...) {
      gsl_failures.report(std::cerr);
      throw;
    }
    gsl_failures.report(std::cerr);
  }
};

} // namespace

#endif
