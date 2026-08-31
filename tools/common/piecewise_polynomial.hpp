#ifndef _tools_common_piecewise_polynomial_hpp_
#define _tools_common_piecewise_polynomial_hpp_

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace NRG::Tools {

using PiecewiseWideFloat = boost::multiprecision::cpp_bin_float_50;
using PiecewiseWideComplex = boost::multiprecision::cpp_complex_50;
using PiecewiseExactRational = boost::multiprecision::cpp_rational;

template <typename Scalar>
inline auto polynomial_scalar_is_finite(const Scalar value) {
  if constexpr (std::is_same_v<Scalar, double>)
    return std::isfinite(value);
  else
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

template <typename Scalar>
inline auto polynomial_scalar_as_complex(const Scalar value) {
  if constexpr (std::is_same_v<Scalar, double>)
    return std::complex<double>{value, 0.0};
  else
    return std::complex<double>{value};
}

template <typename Scalar>
inline auto polynomial_scalar_as_long_complex(const Scalar value) {
  if constexpr (std::is_same_v<Scalar, double>)
    return std::complex<long double>{static_cast<long double>(value), 0.0L};
  else
    return std::complex<long double>{static_cast<long double>(value.real()), static_cast<long double>(value.imag())};
}

namespace detail {

using ExactRational = PiecewiseExactRational;
inline constexpr bool native_extended_precision = std::numeric_limits<long double>::digits
                                                     > std::numeric_limits<double>::digits
                                                   && std::numeric_limits<long double>::max_exponent
                                                        >= 3 * std::numeric_limits<double>::max_exponent
                                                   && std::numeric_limits<long double>::min_exponent
                                                        <= 3 * std::numeric_limits<double>::min_exponent;
inline constexpr std::size_t double_exponent_span = std::numeric_limits<double>::max_exponent
                                                    - std::numeric_limits<double>::min_exponent
                                                    + std::numeric_limits<double>::digits + 8;

inline auto exact_binary_scale(unsigned exponent) {
  boost::multiprecision::cpp_int result{1};
  result <<= exponent;
  return result;
}

inline auto exact_rational(const double value) {
  if (value == 0.0) return ExactRational{};
  int exponent = 0;
  const auto fraction = std::frexp(value, &exponent);
  constexpr auto digits = std::numeric_limits<double>::digits;
  const auto mantissa = static_cast<std::int64_t>(std::ldexp(fraction, digits));
  ExactRational result{mantissa};
  if (exponent >= digits)
    result *= exact_binary_scale(static_cast<unsigned>(exponent - digits));
  else
    result /= exact_binary_scale(static_cast<unsigned>(digits - exponent));
  return result;
}

inline auto exact_power(ExactRational base, std::size_t exponent) {
  ExactRational result{1};
  while (exponent != 0) {
    if (exponent % 2 != 0) result *= base;
    exponent /= 2;
    if (exponent != 0) base *= base;
  }
  return result;
}

} // namespace detail

template <typename Real>
class CompensatedComplexSumType {
  Real real_sum = 0.0;
  Real real_correction = 0.0;
  Real imaginary_sum = 0.0;
  Real imaginary_correction = 0.0;
  Real real_absolute_sum = 0.0;
  Real imaginary_absolute_sum = 0.0;

  static void add_component(const Real value, Real &sum, Real &correction) {
    const auto updated = sum + value;
    if (std::abs(sum) >= std::abs(value))
      correction += (sum - updated) + value;
    else
      correction += (value - updated) + sum;
    sum = updated;
  }

  public:
  void add(const std::complex<Real> value) {
    add_component(value.real(), real_sum, real_correction);
    add_component(value.imag(), imaginary_sum, imaginary_correction);
    real_absolute_sum += std::abs(value.real());
    imaginary_absolute_sum += std::abs(value.imag());
  }

  auto value() const { return std::complex<Real>{real_sum + real_correction, imaginary_sum + imaginary_correction}; }
  auto real_absolute_bound() const { return real_absolute_sum; }
  auto imaginary_absolute_bound() const { return imaginary_absolute_sum; }
};

using CompensatedComplexSum = CompensatedComplexSumType<double>;
using CompensatedLongComplexSum = CompensatedComplexSumType<long double>;

namespace detail {

inline auto normalized_interval_coordinate(const double value, const double left, const double right) {
  if (value == left) return 0.0L;
  if (value == right) return 1.0L;
  const auto extended_value = static_cast<long double>(value);
  const auto extended_left = static_cast<long double>(left);
  const auto extended_right = static_cast<long double>(right);
  const auto width = extended_right - extended_left;
  if (extended_value - extended_left <= extended_right - extended_value)
    return (extended_value - extended_left) / width;
  return 1.0L + (extended_value - extended_right) / width;
}

inline auto wide_normalized_interval_coordinate(const double value, const double left, const double right) {
  if (value == left) return PiecewiseWideFloat{0};
  if (value == right) return PiecewiseWideFloat{1};
  const PiecewiseWideFloat wide_value{value};
  const PiecewiseWideFloat wide_left{left};
  const PiecewiseWideFloat wide_right{right};
  const auto width = wide_right - wide_left;
  if (wide_value - wide_left <= wide_right - wide_value) return (wide_value - wide_left) / width;
  return PiecewiseWideFloat{1} + (wide_value - wide_right) / width;
}

inline auto normalized_primitive_average(const std::size_t power, const long double lower, const long double upper) {
  const auto exponent = power + 1;
  auto lower_power = 1.0L;
  auto divided_difference = 1.0L;
  for (std::size_t term = 1; term < exponent; ++term) {
    lower_power *= lower;
    divided_difference = upper * divided_difference + lower_power;
  }
  return divided_difference / static_cast<double>(exponent);
}

template <typename Scalar>
inline auto normalized_polynomial_integral(const std::vector<Scalar> &coefficients, const long double lower,
                                            const long double upper, const long double physical_span) {
  CompensatedLongComplexSum sum;
  for (std::size_t power = 0; power < coefficients.size(); ++power)
    sum.add(polynomial_scalar_as_long_complex(coefficients[power])
            * (physical_span * normalized_primitive_average(power, lower, upper)));
  return sum.value();
}

inline auto wide_normalized_primitive_average(const std::size_t power, const PiecewiseWideFloat &lower,
                                               const PiecewiseWideFloat &upper) {
  const auto exponent = power + 1;
  PiecewiseWideFloat lower_power{1};
  PiecewiseWideFloat divided_difference{1};
  for (std::size_t term = 1; term < exponent; ++term) {
    lower_power *= lower;
    divided_difference = upper * divided_difference + lower_power;
  }
  return divided_difference / PiecewiseWideFloat{exponent};
}

template <typename Scalar>
inline auto wide_normalized_polynomial_integral(const std::vector<Scalar> &coefficients,
                                                 const PiecewiseWideFloat &lower,
                                                 const PiecewiseWideFloat &upper,
                                                 const PiecewiseWideFloat &physical_span) {
  PiecewiseWideComplex sum;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
    sum += PiecewiseWideComplex{coefficient.real(), coefficient.imag()}
           * (physical_span * wide_normalized_primitive_average(power, lower, upper));
  }
  return sum;
}

} // namespace detail

struct PiecewiseNativeMoment {
  std::complex<long double> value;
  long double real_absolute_bound = 0.0L;
  long double imaginary_absolute_bound = 0.0L;
};

struct PiecewiseFarMomentCache {
  std::mutex mutex;
  std::vector<PiecewiseNativeMoment> native_moments;
  std::vector<std::vector<long double>> native_powers;
  std::vector<std::pair<PiecewiseWideFloat, PiecewiseWideFloat>> moments;
  std::vector<std::vector<PiecewiseWideFloat>> powers;
  std::vector<std::pair<PiecewiseExactRational, PiecewiseExactRational>> exact_moments;
  std::vector<std::vector<PiecewiseExactRational>> exact_powers;
};

template <typename Scalar>
class PiecewisePolynomial {
  static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, std::complex<double>>);

  std::vector<double> knots_;
  std::vector<std::vector<Scalar>> coefficients_;
  mutable std::shared_ptr<PiecewiseFarMomentCache> far_moment_cache_ = std::make_shared<PiecewiseFarMomentCache>();

  auto interval_index(const double x) const {
    if (!std::isfinite(x) || x < knots_.front() || x > knots_.back())
      throw std::domain_error("Polynomial argument is outside the interpolation domain.");
    if (x == knots_.back()) return coefficients_.size() - 1;
    const auto upper = std::upper_bound(knots_.begin(), knots_.end(), x);
    return static_cast<std::size_t>(std::distance(knots_.begin(), upper) - 1);
  }

  public:
  using scalar_type = Scalar;

  PiecewisePolynomial(std::vector<double> knots, std::vector<std::vector<Scalar>> coefficients)
      : knots_{std::move(knots)}, coefficients_{std::move(coefficients)} {
    if (knots_.size() < 2) throw std::invalid_argument("A piecewise polynomial requires at least two knots.");
    if (coefficients_.size() + 1 != knots_.size())
      throw std::invalid_argument("Piecewise-polynomial intervals must match the knot grid.");
    for (std::size_t index = 0; index < knots_.size(); ++index) {
      if (!std::isfinite(knots_[index])) throw std::invalid_argument("Piecewise-polynomial knots must be finite.");
      if (index != 0 && !(knots_[index - 1] < knots_[index]))
        throw std::invalid_argument("Piecewise-polynomial knots must be strictly increasing.");
    }
    if (!std::isfinite(knots_.back() - knots_.front()))
      throw std::invalid_argument("Piecewise-polynomial domain width must be finite.");
    for (const auto &interval : coefficients_) {
      if (interval.empty()) throw std::invalid_argument("Each polynomial interval requires at least one coefficient.");
      for (const auto coefficient : interval)
        if (!polynomial_scalar_is_finite(coefficient)) throw std::invalid_argument("Polynomial coefficients must be finite.");
    }
  }

  const auto &knots() const noexcept { return knots_; }
  const auto &coefficients() const noexcept { return coefficients_; }
  auto interval_count() const noexcept { return coefficients_.size(); }
  auto lower_bound() const noexcept { return knots_.front(); }
  auto upper_bound() const noexcept { return knots_.back(); }

  auto far_moment_cache() const { return far_moment_cache_; }

  auto evaluate(const double x) const {
    const auto index = interval_index(x);
    const auto width = knots_[index + 1] - knots_[index];
    const auto u = (x - knots_[index]) / width;
    Scalar result{};
    for (auto coefficient = coefficients_[index].rbegin(); coefficient != coefficients_[index].rend(); ++coefficient)
      result = result * u + *coefficient;
    return result;
  }

  auto integral(const double lower, const double upper) const {
    if (!std::isfinite(lower) || !std::isfinite(upper))
      throw std::invalid_argument("Polynomial integration bounds must be finite.");
    if (lower > upper) throw std::invalid_argument("Polynomial integration bounds must be ordered.");
    if (lower < knots_.front() || upper > knots_.back())
      throw std::domain_error("Polynomial integration bounds are outside the interpolation domain.");
    if (lower == upper) return Scalar{};

    if constexpr (detail::native_extended_precision) {
      CompensatedLongComplexSum sum;
      for (std::size_t interval = 0; interval < coefficients_.size(); ++interval) {
        const auto left = std::max(lower, knots_[interval]);
        const auto right = std::min(upper, knots_[interval + 1]);
        if (!(left < right)) continue;
        const auto interval_left = knots_[interval];
        const auto interval_right = knots_[interval + 1];
        const auto normalized_left = detail::normalized_interval_coordinate(left, interval_left, interval_right);
        const auto normalized_right = detail::normalized_interval_coordinate(right, interval_left, interval_right);
        const auto physical_span = static_cast<long double>(right) - static_cast<long double>(left);
        sum.add(detail::normalized_polynomial_integral(coefficients_[interval], normalized_left, normalized_right,
                                                       physical_span));
      }
      if constexpr (std::is_same_v<Scalar, double>)
        return static_cast<double>(sum.value().real());
      else
        return std::complex<double>{static_cast<double>(sum.value().real()), static_cast<double>(sum.value().imag())};
    } else {
      PiecewiseWideComplex sum;
      for (std::size_t interval = 0; interval < coefficients_.size(); ++interval) {
        const auto left = std::max(lower, knots_[interval]);
        const auto right = std::min(upper, knots_[interval + 1]);
        if (!(left < right)) continue;
        const auto interval_left = knots_[interval];
        const auto interval_right = knots_[interval + 1];
        const auto normalized_left = detail::wide_normalized_interval_coordinate(left, interval_left, interval_right);
        const auto normalized_right = detail::wide_normalized_interval_coordinate(right, interval_left, interval_right);
        const auto physical_span = PiecewiseWideFloat{right} - PiecewiseWideFloat{left};
        sum += detail::wide_normalized_polynomial_integral(coefficients_[interval], normalized_left, normalized_right,
                                                           physical_span);
      }
      if constexpr (std::is_same_v<Scalar, double>)
        return static_cast<double>(sum.real());
      else
        return std::complex<double>{static_cast<double>(sum.real()), static_cast<double>(sum.imag())};
    }
  }

  auto integral() const { return integral(lower_bound(), upper_bound()); }

  auto multiply_by_monomial(const std::size_t exponent) const {
    if (exponent == 0) return *this;
    if (exponent == std::numeric_limits<std::size_t>::max())
      throw std::length_error("Energy-weighted polynomial degree is too large.");
    std::vector<std::vector<Scalar>> weighted;
    weighted.reserve(coefficients_.size());
    for (std::size_t interval = 0; interval < coefficients_.size(); ++interval) {
      if (coefficients_[interval].size() > std::numeric_limits<std::size_t>::max() - exponent)
        throw std::length_error("Energy-weighted polynomial degree is too large.");
      if (std::all_of(coefficients_[interval].begin(), coefficients_[interval].end(),
                      [](const auto coefficient) { return coefficient == Scalar{}; })) {
        weighted.push_back({Scalar{}});
        continue;
      }

      const auto left = knots_[interval];
      const auto width = knots_[interval + 1] - left;
      const auto extended_left = static_cast<long double>(left);
      const auto extended_width = static_cast<long double>(width);
      std::vector<long double> monomial(exponent + 1, 0.0L);
      if (left == 0.0) {
        monomial[exponent] = std::pow(extended_width, static_cast<long double>(exponent));
      } else if (std::abs(left) >= width) {
        monomial[0] = std::pow(extended_left, static_cast<long double>(exponent));
        for (std::size_t power = 0; power < exponent; ++power)
          monomial[power + 1] = monomial[power] * static_cast<long double>(exponent - power)
                                / static_cast<long double>(power + 1) * extended_width / extended_left;
      } else {
        monomial[exponent] = std::pow(extended_width, static_cast<long double>(exponent));
        for (std::size_t power = exponent; power > 0; --power)
          monomial[power - 1] = monomial[power] * static_cast<long double>(power)
                                / static_cast<long double>(exponent - power + 1) * extended_left / extended_width;
      }

      auto requires_exact_arithmetic = false;
      for (std::size_t power = 0; power < monomial.size(); ++power) {
        const auto expected_nonzero = left != 0.0 || power == exponent;
        if (!std::isfinite(monomial[power]) || (expected_nonzero && monomial[power] == 0.0L)) {
          requires_exact_arithmetic = true;
          break;
        }
      }
      std::vector<std::complex<long double>> extended_product(coefficients_[interval].size() + exponent);
      if (!requires_exact_arithmetic) {
        std::vector<CompensatedLongComplexSum> product_sums(extended_product.size());
        std::vector<long double> real_logarithmic_bounds(extended_product.size(),
                                                         -std::numeric_limits<long double>::infinity());
        std::vector<long double> imaginary_logarithmic_bounds(extended_product.size(),
                                                              -std::numeric_limits<long double>::infinity());
        auto add_logarithm = [](long double &sum, const long double value) {
          if (value == 0.0L) return;
          const auto logarithm = std::log(std::abs(value));
          if (sum == -std::numeric_limits<long double>::infinity()) {
            sum = logarithm;
            return;
          }
          const auto larger = std::max(sum, logarithm);
          sum = larger + std::log1p(std::exp(std::min(sum, logarithm) - larger));
        };
        for (std::size_t first = 0; first < coefficients_[interval].size(); ++first) {
          const auto coefficient = polynomial_scalar_as_long_complex(coefficients_[interval][first]);
          for (std::size_t second = 0; second < monomial.size(); ++second) {
            const auto contribution = coefficient * monomial[second];
            if (monomial[second] != 0.0L
                && ((coefficient.real() != 0.0L && contribution.real() == 0.0L)
                    || (coefficient.imag() != 0.0L && contribution.imag() == 0.0L)))
              requires_exact_arithmetic = true;
            product_sums[first + second].add(contribution);
            add_logarithm(real_logarithmic_bounds[first + second], contribution.real());
            add_logarithm(imaginary_logarithmic_bounds[first + second], contribution.imag());
          }
        }
        const auto logarithmic_roundoff = std::log(8.0L * std::numeric_limits<long double>::epsilon()
                                                    / std::numeric_limits<double>::epsilon());
        for (std::size_t power = 0; power < extended_product.size(); ++power) {
          extended_product[power] = product_sums[power].value();
          const auto uncertain = [logarithmic_roundoff](const long double value, const long double bound) {
            if (!std::isfinite(value)) return true;
            if (bound == -std::numeric_limits<long double>::infinity()) return false;
            return value == 0.0L || std::log(std::abs(value)) <= logarithmic_roundoff + bound;
          };
          if (uncertain(extended_product[power].real(), real_logarithmic_bounds[power])
              || uncertain(extended_product[power].imag(), imaginary_logarithmic_bounds[power])) {
            requires_exact_arithmetic = true;
            break;
          }
        }
      }

      std::vector<Scalar> product;
      product.reserve(extended_product.size());
      if (!requires_exact_arithmetic) {
        for (const auto coefficient : extended_product) {
          const auto real_part = static_cast<double>(coefficient.real());
          const auto imaginary_part = static_cast<double>(coefficient.imag());
          if (!std::isfinite(real_part) || !std::isfinite(imaginary_part)) {
            requires_exact_arithmetic = true;
            break;
          }
          if constexpr (std::is_same_v<Scalar, double>)
            product.push_back(real_part);
          else
            product.emplace_back(real_part, imaginary_part);
        }
      }

      if (requires_exact_arithmetic) {
        std::vector<detail::ExactRational> exact_monomial(exponent + 1);
        const auto exact_left = detail::exact_rational(left);
        const auto exact_width = detail::exact_rational(width);
        if (left == 0.0) {
          exact_monomial[exponent] = detail::exact_power(exact_width, exponent);
        } else if (std::abs(left) >= width) {
          exact_monomial[0] = detail::exact_power(exact_left, exponent);
          for (std::size_t power = 0; power < exponent; ++power)
            exact_monomial[power + 1] = exact_monomial[power] * (exponent - power) / (power + 1)
                                          * exact_width / exact_left;
        } else {
          exact_monomial[exponent] = detail::exact_power(exact_width, exponent);
          for (std::size_t power = exponent; power > 0; --power)
            exact_monomial[power - 1] = exact_monomial[power] * power / (exponent - power + 1)
                                        * exact_left / exact_width;
        }

        std::vector<detail::ExactRational> exact_real(extended_product.size());
        std::vector<detail::ExactRational> exact_imaginary(extended_product.size());
        for (std::size_t first = 0; first < coefficients_[interval].size(); ++first) {
          const auto coefficient = polynomial_scalar_as_complex(coefficients_[interval][first]);
          for (std::size_t second = 0; second < exact_monomial.size(); ++second) {
            exact_real[first + second] += detail::exact_rational(coefficient.real()) * exact_monomial[second];
            exact_imaginary[first + second] += detail::exact_rational(coefficient.imag()) * exact_monomial[second];
          }
        }

        product.clear();
        for (std::size_t power = 0; power < exact_real.size(); ++power) {
          const auto real_part = exact_real[power].convert_to<double>();
          const auto imaginary_part = exact_imaginary[power].convert_to<double>();
          if (!std::isfinite(real_part) || !std::isfinite(imaginary_part))
            throw std::overflow_error("Energy weighting produced a nonfinite polynomial coefficient.");
          if constexpr (std::is_same_v<Scalar, double>)
            product.push_back(real_part);
          else
            product.emplace_back(real_part, imaginary_part);
        }
      }
      weighted.push_back(std::move(product));
    }
    return PiecewisePolynomial<Scalar>{knots_, std::move(weighted)};
  }
};

namespace detail {

inline auto evaluate_real_polynomial(const std::vector<long double> &coefficients, const long double argument) {
  auto result = 0.0L;
  for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
    result = std::fma(result, argument, *coefficient);
  return result;
}

inline auto isolate_real_polynomial_roots(std::vector<long double> coefficients, const long double lower,
                                          const long double upper) {
  while (!coefficients.empty() && coefficients.back() == 0.0L) coefficients.pop_back();
  std::vector<long double> roots;
  if (coefficients.size() <= 1 || lower > upper) return roots;

  std::vector<long double> derivative(coefficients.size() - 1);
  for (std::size_t power = 1; power < coefficients.size(); ++power)
    derivative[power - 1] = static_cast<long double>(power) * coefficients[power];
  auto critical_points = isolate_real_polynomial_roots(std::move(derivative), lower, upper);
  std::sort(critical_points.begin(), critical_points.end());
  critical_points.erase(std::unique(critical_points.begin(), critical_points.end()), critical_points.end());

  std::vector<long double> partition{lower};
  for (const auto point : critical_points)
    if (lower < point && point < upper) partition.push_back(point);
  if (lower < upper) partition.push_back(upper);

  std::vector<long double> values;
  values.reserve(partition.size());
  for (const auto point : partition) values.push_back(evaluate_real_polynomial(coefficients, point));
  for (std::size_t point = 0; point < partition.size(); ++point) {
    if (values[point] == 0.0L) roots.push_back(partition[point]);
    if (point + 1 == partition.size() || values[point] == 0.0L || values[point + 1] == 0.0L
        || std::signbit(values[point]) == std::signbit(values[point + 1]))
      continue;

    auto left = partition[point];
    auto right = partition[point + 1];
    auto left_value = values[point];
    auto right_value = values[point + 1];
    for (int iteration = 0; iteration < std::numeric_limits<long double>::digits + 4; ++iteration) {
      const auto midpoint = left + (right - left) / 2.0L;
      if (midpoint == left || midpoint == right) break;
      const auto midpoint_value = evaluate_real_polynomial(coefficients, midpoint);
      if (midpoint_value == 0.0L) {
        left = midpoint;
        right = midpoint;
        left_value = 0.0L;
        right_value = 0.0L;
        break;
      }
      if (std::signbit(left_value) == std::signbit(midpoint_value)) {
        left = midpoint;
        left_value = midpoint_value;
      } else {
        right = midpoint;
        right_value = midpoint_value;
      }
    }
    roots.push_back(std::abs(left_value) <= std::abs(right_value) ? left : right);
  }
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
  return roots;
}

inline auto real_polynomial_roots(const std::vector<double> &coefficients, const long double lower,
                                  const long double upper) {
  auto scale = 0.0L;
  for (const auto coefficient : coefficients)
    scale = std::max(scale, std::abs(static_cast<long double>(coefficient)));
  if (scale == 0.0L) return std::vector<long double>{};

  std::vector<long double> scaled_coefficients;
  scaled_coefficients.reserve(coefficients.size());
  for (const auto coefficient : coefficients)
    scaled_coefficients.push_back(static_cast<long double>(coefficient) / scale);
  const auto extended_roots = isolate_real_polynomial_roots(std::move(scaled_coefficients), lower, upper);
  std::vector<long double> roots;
  roots.reserve(extended_roots.size());
  for (const auto root : extended_roots) roots.push_back(std::clamp(root, lower, upper));
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
  return roots;
}

inline auto evaluate_wide_real_polynomial(const std::vector<PiecewiseWideFloat> &coefficients,
                                           const PiecewiseWideFloat &argument) {
  PiecewiseWideFloat result{};
  for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
    result = result * argument + *coefficient;
  return result;
}

inline auto isolate_wide_real_polynomial_roots(std::vector<PiecewiseWideFloat> coefficients,
                                                const PiecewiseWideFloat &lower,
                                                const PiecewiseWideFloat &upper) {
  while (!coefficients.empty() && coefficients.back() == 0) coefficients.pop_back();
  std::vector<PiecewiseWideFloat> roots;
  if (coefficients.size() <= 1 || lower > upper) return roots;

  std::vector<PiecewiseWideFloat> derivative(coefficients.size() - 1);
  for (std::size_t power = 1; power < coefficients.size(); ++power)
    derivative[power - 1] = PiecewiseWideFloat{power} * coefficients[power];
  auto critical_points = isolate_wide_real_polynomial_roots(std::move(derivative), lower, upper);
  std::sort(critical_points.begin(), critical_points.end());
  critical_points.erase(std::unique(critical_points.begin(), critical_points.end()), critical_points.end());

  std::vector<PiecewiseWideFloat> partition{lower};
  for (const auto &point : critical_points)
    if (lower < point && point < upper) partition.push_back(point);
  if (lower < upper) partition.push_back(upper);

  std::vector<PiecewiseWideFloat> values;
  values.reserve(partition.size());
  for (const auto &point : partition) values.push_back(evaluate_wide_real_polynomial(coefficients, point));
  for (std::size_t point = 0; point < partition.size(); ++point) {
    if (values[point] == 0) roots.push_back(partition[point]);
    if (point + 1 == partition.size() || values[point] == 0 || values[point + 1] == 0
        || (values[point] < 0) == (values[point + 1] < 0))
      continue;

    auto left = partition[point];
    auto right = partition[point + 1];
    auto left_value = values[point];
    auto right_value = values[point + 1];
    for (int iteration = 0; iteration < std::numeric_limits<PiecewiseWideFloat>::digits + 4; ++iteration) {
      const auto midpoint = left + (right - left) / 2;
      if (midpoint == left || midpoint == right) break;
      const auto midpoint_value = evaluate_wide_real_polynomial(coefficients, midpoint);
      if (midpoint_value == 0) {
        left = midpoint;
        right = midpoint;
        left_value = 0;
        right_value = 0;
        break;
      }
      if ((left_value < 0) == (midpoint_value < 0)) {
        left = midpoint;
        left_value = midpoint_value;
      } else {
        right = midpoint;
        right_value = midpoint_value;
      }
    }
    roots.push_back(boost::multiprecision::abs(left_value) <= boost::multiprecision::abs(right_value) ? left : right);
  }
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
  return roots;
}

inline auto wide_real_polynomial_roots(const std::vector<double> &coefficients, const PiecewiseWideFloat &lower,
                                       const PiecewiseWideFloat &upper) {
  PiecewiseWideFloat scale{};
  for (const auto coefficient : coefficients)
    scale = std::max(scale, boost::multiprecision::abs(PiecewiseWideFloat{coefficient}));
  if (scale == 0) return std::vector<PiecewiseWideFloat>{};

  std::vector<PiecewiseWideFloat> scaled_coefficients;
  scaled_coefficients.reserve(coefficients.size());
  for (const auto coefficient : coefficients) scaled_coefficients.push_back(PiecewiseWideFloat{coefficient} / scale);
  return isolate_wide_real_polynomial_roots(std::move(scaled_coefficients), lower, upper);
}

} // namespace detail

inline auto absolute_integral(const PiecewisePolynomial<double> &polynomial, const double lower, const double upper) {
  if (!std::isfinite(lower) || !std::isfinite(upper))
    throw std::invalid_argument("Absolute-integration bounds must be finite.");
  if (lower > upper) throw std::invalid_argument("Absolute-integration bounds must be ordered.");
  if (lower < polynomial.lower_bound() || upper > polynomial.upper_bound())
    throw std::domain_error("Absolute-integration bounds are outside the interpolation domain.");
  for (const auto &coefficients : polynomial.coefficients())
    if (coefficients.size() > 4)
      throw std::invalid_argument("Absolute integration supports polynomial intervals of degree at most three.");
  if (lower == upper) return 0.0;

  if constexpr (detail::native_extended_precision) {
    CompensatedLongComplexSum sum;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left = std::max(lower, polynomial.knots()[interval]);
      const auto right = std::min(upper, polynomial.knots()[interval + 1]);
      if (!(left < right)) continue;
      const auto interval_left = polynomial.knots()[interval];
      const auto interval_right = polynomial.knots()[interval + 1];
      const auto width = static_cast<long double>(interval_right) - static_cast<long double>(interval_left);
      const auto normalized_left = detail::normalized_interval_coordinate(left, interval_left, interval_right);
      const auto normalized_right = detail::normalized_interval_coordinate(right, interval_left, interval_right);
      const auto roots = detail::real_polynomial_roots(polynomial.coefficients()[interval], normalized_left,
                                                       normalized_right);
      auto segment_left = normalized_left;
      for (const auto root : roots) {
        if (!(segment_left < root && root < normalized_right)) continue;
        const auto value = detail::normalized_polynomial_integral(polynomial.coefficients()[interval], segment_left,
                                                                  root, width * (root - segment_left)).real();
        sum.add({std::abs(value), 0.0L});
        segment_left = root;
      }
      const auto physical_span = segment_left == normalized_left
                                   ? static_cast<long double>(right) - static_cast<long double>(left)
                                   : width * (normalized_right - segment_left);
      const auto value = detail::normalized_polynomial_integral(polynomial.coefficients()[interval], segment_left,
                                                                normalized_right, physical_span).real();
      sum.add({std::abs(value), 0.0L});
    }
    return static_cast<double>(sum.value().real());
  } else {
    PiecewiseWideFloat sum{};
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left = std::max(lower, polynomial.knots()[interval]);
      const auto right = std::min(upper, polynomial.knots()[interval + 1]);
      if (!(left < right)) continue;
      const auto interval_left = polynomial.knots()[interval];
      const auto interval_right = polynomial.knots()[interval + 1];
      const auto width = PiecewiseWideFloat{interval_right} - PiecewiseWideFloat{interval_left};
      const auto normalized_left = detail::wide_normalized_interval_coordinate(left, interval_left, interval_right);
      const auto normalized_right = detail::wide_normalized_interval_coordinate(right, interval_left, interval_right);
      const auto roots = detail::wide_real_polynomial_roots(polynomial.coefficients()[interval], normalized_left,
                                                            normalized_right);
      auto segment_left = normalized_left;
      for (const auto &root : roots) {
        if (!(segment_left < root && root < normalized_right)) continue;
        const auto value = detail::wide_normalized_polynomial_integral(
          polynomial.coefficients()[interval], segment_left, root, width * (root - segment_left)).real();
        sum += boost::multiprecision::abs(value);
        segment_left = root;
      }
      const auto physical_span = segment_left == normalized_left
                                   ? PiecewiseWideFloat{right} - PiecewiseWideFloat{left}
                                   : width * (normalized_right - segment_left);
      const auto value = detail::wide_normalized_polynomial_integral(
        polynomial.coefficients()[interval], segment_left, normalized_right, physical_span).real();
      sum += boost::multiprecision::abs(value);
    }
    return static_cast<double>(sum);
  }
}

inline auto absolute_integral(const PiecewisePolynomial<double> &polynomial) {
  return absolute_integral(polynomial, polynomial.lower_bound(), polynomial.upper_bound());
}

inline auto combine_piecewise_polynomials(const PiecewisePolynomial<double> &real_part,
                                          const PiecewisePolynomial<double> &imaginary_part) {
  if (real_part.knots() != imaginary_part.knots())
    throw std::invalid_argument("Real and imaginary piecewise polynomials must use the same knots.");
  std::vector<std::vector<std::complex<double>>> coefficients;
  coefficients.reserve(real_part.interval_count());
  for (std::size_t interval = 0; interval < real_part.interval_count(); ++interval) {
    const auto &real_coefficients = real_part.coefficients()[interval];
    const auto &imaginary_coefficients = imaginary_part.coefficients()[interval];
    const auto size = std::max(real_coefficients.size(), imaginary_coefficients.size());
    std::vector<std::complex<double>> combined(size);
    for (std::size_t power = 0; power < size; ++power) {
      const auto real_value = power < real_coefficients.size() ? real_coefficients[power] : 0.0;
      const auto imaginary_value = power < imaginary_coefficients.size() ? imaginary_coefficients[power] : 0.0;
      combined[power] = {real_value, imaginary_value};
    }
    coefficients.push_back(std::move(combined));
  }
  return PiecewisePolynomial<std::complex<double>>{real_part.knots(), std::move(coefficients)};
}

enum class CauchyEndpointPolicy { reject, subtracted };

namespace detail {

using FarReal = PiecewiseWideFloat;

struct FarComplex {
  FarReal real;
  FarReal imaginary;
};

struct ExactComplex {
  ExactRational real;
  ExactRational imaginary;
};

struct ResolvedMoment {
  long double real = 0.0L;
  long double imaginary = 0.0L;
  bool wide_real = false;
  bool wide_imaginary = false;
  FarReal real_wide;
  FarReal imaginary_wide;
};

inline void add(FarComplex &sum, const FarComplex &value) {
  sum.real += value.real;
  sum.imaginary += value.imaginary;
}

inline auto multiply(const FarComplex &first, const FarComplex &second) {
  return FarComplex{first.real * second.real - first.imaginary * second.imaginary,
                    first.real * second.imaginary + first.imaginary * second.real};
}

inline auto multiply(const FarComplex &value, const FarReal &scale) {
  return FarComplex{value.real * scale, value.imaginary * scale};
}

inline auto far_inverse_root_from_argument(const double right, const double left, const double argument_real,
                                           const double argument_imaginary) {
  const FarReal real = FarReal{argument_real} - FarReal{left};
  const FarReal imaginary{argument_imaginary};
  if (real == 0 && imaginary == 0) return FarComplex{};
  const auto norm = real * real + imaginary * imaginary;
  const auto numerator = FarReal{right} - FarReal{left};
  return FarComplex{numerator * real / norm, -numerator * imaginary / norm};
}

inline auto exact_inverse_root_from_argument(const double right, const double left, const double argument_real,
                                             const double argument_imaginary) {
  const auto real = exact_rational(argument_real) - exact_rational(left);
  const auto imaginary = exact_rational(argument_imaginary);
  const auto norm = real * real + imaginary * imaginary;
  const auto numerator = exact_rational(right) - exact_rational(left);
  return ExactComplex{numerator * real / norm, -numerator * imaginary / norm};
}

inline auto multiply(const ExactComplex &first, const ExactComplex &second) {
  return ExactComplex{first.real * second.real - first.imaginary * second.imaginary,
                      first.real * second.imaginary + first.imaginary * second.real};
}

inline auto far_power(FarComplex base, std::size_t exponent) {
  FarComplex result{FarReal{1}, FarReal{0}};
  while (exponent != 0) {
    if (exponent % 2 != 0) result = multiply(result, base);
    exponent /= 2;
    if (exponent != 0) base = multiply(base, base);
  }
  return result;
}

inline auto log_add_positive(const long double first, const long double second) {
  if (first == -std::numeric_limits<long double>::infinity()) return second;
  if (second == -std::numeric_limits<long double>::infinity()) return first;
  const auto larger = std::max(first, second);
  return larger + std::log1p(std::exp(std::min(first, second) - larger));
}

inline auto logarithmic_absolute(const FarReal &value) {
  if (value == 0) return -std::numeric_limits<long double>::infinity();
  int exponent = 0;
  const auto fraction = boost::multiprecision::frexp(boost::multiprecision::abs(value), &exponent);
  return static_cast<long double>(exponent) * std::log(2.0L)
         + std::log(fraction.convert_to<long double>());
}

inline auto logarithmic_magnitude(const FarComplex &value) {
  const auto real_logarithm = logarithmic_absolute(value.real);
  const auto imaginary_logarithm = logarithmic_absolute(value.imaginary);
  return 0.5L * log_add_positive(2.0L * real_logarithm, 2.0L * imaginary_logarithm);
}

inline auto logarithmic_magnitude(const std::complex<long double> value) {
  const auto scale = std::max(std::abs(value.real()), std::abs(value.imag()));
  if (scale == 0.0L) return -std::numeric_limits<long double>::infinity();
  const auto real = value.real() / scale;
  const auto imaginary = value.imag() / scale;
  return std::log(scale) + 0.5L * std::log(real * real + imaginary * imaginary);
}

inline auto log_absolute(const long double value) {
  return value == 0.0L ? -std::numeric_limits<long double>::infinity() : std::log(std::abs(value));
}

inline auto log_absolute_difference(const double first, const double second) {
  if (first == second) return -std::numeric_limits<long double>::infinity();
  if (std::signbit(first) == std::signbit(second))
    return log_absolute(static_cast<long double>(first) - static_cast<long double>(second));
  return log_add_positive(log_absolute(static_cast<long double>(first)), log_absolute(static_cast<long double>(second)));
}

inline auto moment_component_is_uncertain(const long double value, const long double logarithmic_term_bound,
                                          const bool require_double_accuracy) {
  if (!std::isfinite(value)) return true;
  if (logarithmic_term_bound == -std::numeric_limits<long double>::infinity()) return false;
  if (value == 0.0L) return true;
  const auto target = require_double_accuracy ? std::numeric_limits<double>::epsilon() : 1.0L;
  return std::log(std::abs(value))
         <= std::log(8.0L * std::numeric_limits<long double>::epsilon() / target) + logarithmic_term_bound;
}

inline auto moment_component_is_uncertain(const FarReal &value, const long double logarithmic_term_bound) {
  if (logarithmic_term_bound == -std::numeric_limits<long double>::infinity()) return false;
  if (value == 0) return true;
  constexpr auto precision = std::numeric_limits<FarReal>::digits;
  const auto logarithmic_roundoff = std::log(8.0L / std::numeric_limits<double>::epsilon())
                                     + static_cast<long double>(1 - precision) * std::log(2.0L);
  return logarithmic_absolute(value) <= logarithmic_roundoff + logarithmic_term_bound;
}

template <typename HighPrecisionMoment, typename ExactMoment>
inline auto resolve_uncertain_moment(const std::complex<long double> value, const long double real_logarithmic_bound,
                                     const long double imaginary_logarithmic_bound, const bool require_double_accuracy,
                                     HighPrecisionMoment high_precision_moment, ExactMoment exact_moment) {
  const auto uncertain_real = moment_component_is_uncertain(value.real(), real_logarithmic_bound,
                                                             require_double_accuracy);
  const auto uncertain_imaginary = moment_component_is_uncertain(value.imag(), imaginary_logarithmic_bound,
                                                                  require_double_accuracy);
  ResolvedMoment result;
  result.real = value.real();
  result.imaginary = value.imag();
  if (!uncertain_real && !uncertain_imaginary) return result;
  const auto high_precision = high_precision_moment();
  auto selected_real = uncertain_real ? high_precision.first : FarReal{value.real()};
  auto selected_imaginary = uncertain_imaginary ? high_precision.second : FarReal{value.imag()};
  const auto exact_real = uncertain_real && moment_component_is_uncertain(selected_real, real_logarithmic_bound);
  const auto exact_imaginary = uncertain_imaginary
                                && moment_component_is_uncertain(selected_imaginary, imaginary_logarithmic_bound);
  if (exact_real || exact_imaginary) {
    const auto exact = exact_moment();
    if (exact_real) selected_real = exact.first.template convert_to<FarReal>();
    if (exact_imaginary) selected_imaginary = exact.second.template convert_to<FarReal>();
  }
  auto resolve_component = [](const FarReal &high_precision_value, long double &narrow, bool &wide,
                              FarReal &wide_value) {
    const auto converted = high_precision_value.convert_to<long double>();
    if (std::isfinite(converted) && !(converted == 0.0L && high_precision_value != 0)) {
      narrow = converted;
      return;
    }
    wide = true;
    wide_value = high_precision_value;
  };
  if (uncertain_real) resolve_component(selected_real, result.real, result.wide_real, result.real_wide);
  if (uncertain_imaginary)
    resolve_component(selected_imaginary, result.imaginary, result.wide_imaginary, result.imaginary_wide);
  return result;
}

template <typename HighPrecisionMoment, typename ExactMoment>
inline auto resolve_high_precision_moment(const long double real_logarithmic_bound,
                                          const long double imaginary_logarithmic_bound,
                                          HighPrecisionMoment high_precision_moment, ExactMoment exact_moment) {
  auto selected = high_precision_moment();
  const auto exact_real = moment_component_is_uncertain(selected.first, real_logarithmic_bound);
  const auto exact_imaginary = moment_component_is_uncertain(selected.second, imaginary_logarithmic_bound);
  if (exact_real || exact_imaginary) {
    const auto exact = exact_moment();
    if (exact_real) selected.first = exact.first.template convert_to<FarReal>();
    if (exact_imaginary) selected.second = exact.second.template convert_to<FarReal>();
  }
  ResolvedMoment result;
  result.wide_real = true;
  result.wide_imaginary = true;
  result.real_wide = std::move(selected.first);
  result.imaginary_wide = std::move(selected.second);
  return result;
}

template <typename Sum>
inline auto far_series_converged(const long double logarithmic_bound, const std::size_t expansion,
                                 const long double inverse_logarithm, const long double inverse_magnitude,
                                 const Sum &sum) {
  if (logarithmic_bound == -std::numeric_limits<long double>::infinity()) return true;
  const auto sum_logarithm = logarithmic_magnitude(sum);
  if (sum_logarithm == -std::numeric_limits<long double>::infinity()) return false;
  const auto logarithmic_remaining_bound = logarithmic_bound
                                           + static_cast<long double>(expansion + 2) * inverse_logarithm
                                           - std::log1p(-inverse_magnitude);
  return logarithmic_remaining_bound
         <= std::log(static_cast<long double>(std::numeric_limits<double>::epsilon())) + sum_logarithm;
}

class FarSeriesAccumulator {
  CompensatedLongComplexSum real_narrow;
  CompensatedLongComplexSum imaginary_narrow;
  FarComplex real_wide;
  FarComplex imaginary_wide;
  FarComplex inverse_root_wide;
  FarComplex inverse_power_wide;
  FarReal wide_output_real_absolute_bound;
  FarReal wide_output_imaginary_absolute_bound;
  bool wide = false;

  void activate_wide(const FarComplex &inverse_root, const std::size_t expansion) {
    wide = true;
    const auto real_value = real_narrow.value();
    const auto imaginary_value = imaginary_narrow.value();
    real_wide = {FarReal{real_value.real()}, FarReal{real_value.imag()}};
    imaginary_wide = {FarReal{imaginary_value.real()}, FarReal{imaginary_value.imag()}};
    wide_output_real_absolute_bound = FarReal{real_narrow.real_absolute_bound()}
                                      + FarReal{imaginary_narrow.imaginary_absolute_bound()};
    wide_output_imaginary_absolute_bound = FarReal{real_narrow.imaginary_absolute_bound()}
                                           + FarReal{imaginary_narrow.real_absolute_bound()};
    inverse_root_wide = inverse_root;
    inverse_power_wide = far_power(inverse_root_wide, expansion + 1);
  }

  public:
  void force_wide(const FarComplex &inverse_root) {
    wide = true;
    inverse_root_wide = inverse_root;
    inverse_power_wide = inverse_root;
  }

  void add_moment(const ResolvedMoment &moment, const std::complex<long double> inverse_root,
                  const std::optional<FarComplex> &inverse_root_extended,
                  const std::complex<long double> inverse_power,
                  const std::size_t expansion) {
    if (!wide && !moment.wide_real && !moment.wide_imaginary) {
      const auto real_term = inverse_power * moment.real;
      const auto imaginary_term = inverse_power * moment.imaginary;
      if (std::isfinite(real_term.real()) && std::isfinite(real_term.imag())
          && std::isfinite(imaginary_term.real()) && std::isfinite(imaginary_term.imag())) {
        real_narrow.add(real_term);
        imaginary_narrow.add(imaginary_term);
        return;
      }
    }
    if (!wide) {
      const auto root = inverse_root_extended
                          ? *inverse_root_extended
                          : FarComplex{FarReal{inverse_root.real()}, FarReal{inverse_root.imag()}};
      activate_wide(root, expansion);
    }
    const auto real_term = multiply(inverse_power_wide, moment.wide_real ? moment.real_wide : FarReal{moment.real});
    const auto imaginary_term = multiply(inverse_power_wide,
                                         moment.wide_imaginary ? moment.imaginary_wide : FarReal{moment.imaginary});
    add(real_wide, real_term);
    add(imaginary_wide, imaginary_term);
    wide_output_real_absolute_bound += boost::multiprecision::abs(real_term.real)
                                       + boost::multiprecision::abs(imaginary_term.imaginary);
    wide_output_imaginary_absolute_bound += boost::multiprecision::abs(real_term.imaginary)
                                            + boost::multiprecision::abs(imaginary_term.real);
  }

  void advance() {
    if (wide) inverse_power_wide = multiply(inverse_power_wide, inverse_root_wide);
  }

  void advance_or_promote(const std::complex<long double> inverse_root,
                          const std::complex<long double> next_inverse_power,
                          const std::size_t next_expansion) {
    if (wide) {
      advance();
      return;
    }
    if (next_inverse_power == std::complex<long double>{} && inverse_root != std::complex<long double>{}) {
      const FarComplex root{FarReal{inverse_root.real()}, FarReal{inverse_root.imag()}};
      activate_wide(root, next_expansion);
    }
  }

  auto output_is_ill_conditioned() const {
    if (wide) {
      const auto real = real_wide.real - imaginary_wide.imaginary;
      const auto imaginary = real_wide.imaginary + imaginary_wide.real;
      const auto threshold = FarReal{8.0 / std::numeric_limits<double>::epsilon()}
                             * std::numeric_limits<FarReal>::epsilon();
      return (wide_output_real_absolute_bound != 0
              && boost::multiprecision::abs(real) <= threshold * wide_output_real_absolute_bound)
             || (wide_output_imaginary_absolute_bound != 0
                 && boost::multiprecision::abs(imaginary)
                      <= threshold * wide_output_imaginary_absolute_bound);
    }
    const auto real_value = real_narrow.value();
    const auto imaginary_value = imaginary_narrow.value();
    constexpr auto threshold = 8.0L * std::numeric_limits<long double>::epsilon()
                               / std::numeric_limits<double>::epsilon();
    auto component_is_ill_conditioned = [](const long double value, const long double absolute_bound) {
      return absolute_bound > 0.0L && std::abs(value) <= threshold * absolute_bound;
    };
    return component_is_ill_conditioned(
             real_value.real() - imaginary_value.imag(),
             real_narrow.real_absolute_bound() + imaginary_narrow.imaginary_absolute_bound())
           || component_is_ill_conditioned(
             real_value.imag() + imaginary_value.real(),
             real_narrow.imaginary_absolute_bound() + imaginary_narrow.real_absolute_bound());
  }

  auto converged(const long double real_logarithmic_bound, const long double imaginary_logarithmic_bound,
                 const long double output_real_logarithmic_bound,
                 const long double output_imaginary_logarithmic_bound,
                 const std::size_t expansion, const long double inverse_logarithm,
                 const long double inverse_magnitude) const {
    const auto independent = wide
                               ? far_series_converged(real_logarithmic_bound, expansion, inverse_logarithm,
                                                      inverse_magnitude, real_wide)
                                   && far_series_converged(imaginary_logarithmic_bound, expansion, inverse_logarithm,
                                                           inverse_magnitude, imaginary_wide)
                               : far_series_converged(real_logarithmic_bound, expansion, inverse_logarithm,
                                                      inverse_magnitude, real_narrow.value())
                                   && far_series_converged(imaginary_logarithmic_bound, expansion, inverse_logarithm,
                                                           inverse_magnitude, imaginary_narrow.value());
    if (!independent) return false;

    long double real_logarithm = 0.0L;
    long double imaginary_logarithm = 0.0L;
    if (wide) {
      real_logarithm = logarithmic_absolute(real_wide.real - imaginary_wide.imaginary);
      imaginary_logarithm = logarithmic_absolute(real_wide.imaginary + imaginary_wide.real);
    } else {
      const auto real_value = real_narrow.value();
      const auto imaginary_value = imaginary_narrow.value();
      real_logarithm = log_absolute(real_value.real() - imaginary_value.imag());
      imaginary_logarithm = log_absolute(real_value.imag() + imaginary_value.real());
    }
    auto component_converged = [expansion, inverse_logarithm, inverse_magnitude](
                                 const long double value_logarithm, const long double logarithmic_bound) {
      if (logarithmic_bound == -std::numeric_limits<long double>::infinity()) return true;
      const auto logarithmic_remaining_bound = logarithmic_bound
                                               + static_cast<long double>(expansion + 2) * inverse_logarithm
                                               - std::log1p(-inverse_magnitude);
      if (value_logarithm == -std::numeric_limits<long double>::infinity())
        return logarithmic_remaining_bound
               <= std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min())) - std::log(2.0L);
      return logarithmic_remaining_bound
             <= std::log(static_cast<long double>(std::numeric_limits<double>::epsilon())) + value_logarithm;
    };
    return component_converged(real_logarithm, output_real_logarithmic_bound)
           && component_converged(imaginary_logarithm, output_imaginary_logarithmic_bound);
  }

  auto value() const {
    if (wide) {
      const FarReal real = real_wide.real - imaginary_wide.imaginary;
      const FarReal imaginary = real_wide.imaginary + imaginary_wide.real;
      return std::complex<double>{real.convert_to<double>(), imaginary.convert_to<double>()};
    }
    const auto real_value = real_narrow.value();
    const auto imaginary_value = imaginary_narrow.value();
    return std::complex<double>{static_cast<double>(real_value.real() - imaginary_value.imag()),
                                static_cast<double>(real_value.imag() + imaginary_value.real())};
  }
};

template <typename Scalar>
inline auto high_precision_local_moment(const std::vector<Scalar> &coefficients,
                                        const std::complex<double> subtraction, const std::size_t expansion) {
  FarReal real;
  FarReal imaginary;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
    const FarReal denominator{power + expansion + 1};
    real += (FarReal{coefficient.real()} - (power == 0 ? FarReal{subtraction.real()} : FarReal{})) / denominator;
    imaginary += (FarReal{coefficient.imag()} - (power == 0 ? FarReal{subtraction.imag()} : FarReal{})) / denominator;
  }
  return std::pair<FarReal, FarReal>{std::move(real), std::move(imaginary)};
}

template <typename Scalar>
inline auto exact_local_moment(const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
                               const std::size_t expansion) {
  ExactRational real;
  ExactRational imaginary;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
    const auto denominator = power + expansion + 1;
    real += (exact_rational(coefficient.real())
             - (power == 0 ? exact_rational(subtraction.real()) : ExactRational{})) / denominator;
    imaginary += (exact_rational(coefficient.imag())
                  - (power == 0 ? exact_rational(subtraction.imag()) : ExactRational{})) / denominator;
  }
  return std::pair<ExactRational, ExactRational>{std::move(real), std::move(imaginary)};
}

template <typename Scalar>
inline auto native_global_moment(const PiecewisePolynomial<Scalar> &polynomial, const std::size_t expansion) {
  const auto cache = polynomial.far_moment_cache();
  const std::lock_guard lock(cache->mutex);
  if (cache->native_powers.empty())
    cache->native_powers.assign(polynomial.interval_count(), std::vector<long double>{1.0L});
  const auto global_left = static_cast<long double>(polynomial.lower_bound());
  const auto global_width = static_cast<long double>(polynomial.upper_bound()) - global_left;
  while (cache->native_moments.size() <= expansion) {
    CompensatedLongComplexSum sum;
    auto real_absolute_bound = 0.0L;
    auto imaginary_absolute_bound = 0.0L;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto width_fraction = (static_cast<long double>(polynomial.knots()[interval + 1])
                                   - static_cast<long double>(polynomial.knots()[interval]))
                                  / global_width;
      const auto &power_coefficients = cache->native_powers[interval];
      for (std::size_t global_power = 0; global_power < power_coefficients.size(); ++global_power) {
        for (std::size_t local_power = 0; local_power < polynomial.coefficients()[interval].size(); ++local_power) {
          const auto coefficient = polynomial_scalar_as_long_complex(polynomial.coefficients()[interval][local_power]);
          const auto factor = width_fraction * power_coefficients[global_power]
                              / static_cast<long double>(local_power + global_power + 1);
          const auto contribution = coefficient * factor;
          sum.add(contribution);
          real_absolute_bound += std::abs(contribution.real());
          imaginary_absolute_bound += std::abs(contribution.imag());
        }
      }
    }
    cache->native_moments.push_back({sum.value(), real_absolute_bound, imaginary_absolute_bound});

    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left_fraction = (static_cast<long double>(polynomial.knots()[interval]) - global_left) / global_width;
      const auto width_fraction = (static_cast<long double>(polynomial.knots()[interval + 1])
                                   - static_cast<long double>(polynomial.knots()[interval]))
                                  / global_width;
      const auto &power_coefficients = cache->native_powers[interval];
      std::vector<long double> next_power(power_coefficients.size() + 1);
      for (std::size_t coefficient = 0; coefficient < power_coefficients.size(); ++coefficient) {
        next_power[coefficient] += left_fraction * power_coefficients[coefficient];
        next_power[coefficient + 1] += width_fraction * power_coefficients[coefficient];
      }
      cache->native_powers[interval] = std::move(next_power);
    }
  }
  return cache->native_moments[expansion];
}

template <typename Scalar>
inline auto high_precision_global_moment(const PiecewisePolynomial<Scalar> &polynomial, const std::size_t expansion) {
  const auto cache = polynomial.far_moment_cache();
  const std::lock_guard lock(cache->mutex);
  if (cache->powers.empty())
    cache->powers.assign(polynomial.interval_count(), std::vector<FarReal>{FarReal{1}});
  const auto global_left = polynomial.lower_bound();
  const auto global_width = FarReal{polynomial.upper_bound()} - FarReal{global_left};
  while (cache->moments.size() <= expansion) {
    FarReal real;
    FarReal imaginary;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto width_fraction = (FarReal{polynomial.knots()[interval + 1]} - FarReal{polynomial.knots()[interval]})
                                  / global_width;
      const auto &power_coefficients = cache->powers[interval];
      for (std::size_t global_power = 0; global_power < power_coefficients.size(); ++global_power) {
        for (std::size_t local_power = 0; local_power < polynomial.coefficients()[interval].size(); ++local_power) {
          const auto coefficient = polynomial_scalar_as_complex(polynomial.coefficients()[interval][local_power]);
          const auto factor = width_fraction * power_coefficients[global_power]
                              / FarReal{local_power + global_power + 1};
          real += FarReal{coefficient.real()} * factor;
          imaginary += FarReal{coefficient.imag()} * factor;
        }
      }
    }
    cache->moments.emplace_back(std::move(real), std::move(imaginary));

    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left_fraction = (FarReal{polynomial.knots()[interval]} - FarReal{global_left}) / global_width;
      const auto width_fraction = (FarReal{polynomial.knots()[interval + 1]} - FarReal{polynomial.knots()[interval]})
                                  / global_width;
      const auto &power_coefficients = cache->powers[interval];
      std::vector<FarReal> next_power(power_coefficients.size() + 1);
      for (std::size_t coefficient = 0; coefficient < power_coefficients.size(); ++coefficient) {
        next_power[coefficient] += left_fraction * power_coefficients[coefficient];
        next_power[coefficient + 1] += width_fraction * power_coefficients[coefficient];
      }
      cache->powers[interval] = std::move(next_power);
    }
  }
  return cache->moments[expansion];
}

template <typename Scalar>
inline auto exact_global_moment(const PiecewisePolynomial<Scalar> &polynomial, const std::size_t expansion) {
  const auto cache = polynomial.far_moment_cache();
  const std::lock_guard lock(cache->mutex);
  if (cache->exact_powers.empty())
    cache->exact_powers.assign(polynomial.interval_count(), std::vector<ExactRational>{ExactRational{1}});
  const auto global_left = polynomial.lower_bound();
  const auto global_width = exact_rational(polynomial.upper_bound()) - exact_rational(global_left);
  while (cache->exact_moments.size() <= expansion) {
    ExactRational real;
    ExactRational imaginary;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto width_fraction = (exact_rational(polynomial.knots()[interval + 1])
                                   - exact_rational(polynomial.knots()[interval]))
                                  / global_width;
      const auto &power_coefficients = cache->exact_powers[interval];
      for (std::size_t global_power = 0; global_power < power_coefficients.size(); ++global_power) {
        for (std::size_t local_power = 0; local_power < polynomial.coefficients()[interval].size(); ++local_power) {
          const auto coefficient = polynomial_scalar_as_complex(polynomial.coefficients()[interval][local_power]);
          const auto factor = width_fraction * power_coefficients[global_power]
                              / (local_power + global_power + 1);
          real += exact_rational(coefficient.real()) * factor;
          imaginary += exact_rational(coefficient.imag()) * factor;
        }
      }
    }
    cache->exact_moments.emplace_back(std::move(real), std::move(imaginary));

    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left_fraction = (exact_rational(polynomial.knots()[interval]) - exact_rational(global_left))
                                 / global_width;
      const auto width_fraction = (exact_rational(polynomial.knots()[interval + 1])
                                   - exact_rational(polynomial.knots()[interval]))
                                  / global_width;
      const auto &power_coefficients = cache->exact_powers[interval];
      std::vector<ExactRational> next_power(power_coefficients.size() + 1);
      for (std::size_t coefficient = 0; coefficient < power_coefficients.size(); ++coefficient) {
        next_power[coefficient] += left_fraction * power_coefficients[coefficient];
        next_power[coefficient + 1] += width_fraction * power_coefficients[coefficient];
      }
      cache->exact_powers[interval] = std::move(next_power);
    }
  }
  return cache->exact_moments[expansion];
}

template <typename ExactMoment>
inline auto exact_far_series(const ExactComplex inverse_root, const long double inverse_logarithm,
                             const long double inverse_magnitude, const long double output_real_logarithmic_bound,
                             const long double output_imaginary_logarithmic_bound,
                             const std::size_t maximum_expansions, ExactMoment exact_moment) {
  ExactComplex sum;
  auto inverse_power = inverse_root;
  auto component_converged = [inverse_logarithm, inverse_magnitude](const ExactRational &value,
                                                                    const long double logarithmic_bound,
                                                                    const std::size_t expansion) {
    if (logarithmic_bound == -std::numeric_limits<long double>::infinity()) return std::pair{true, true};
    const auto logarithmic_remaining_bound = logarithmic_bound
                                             + static_cast<long double>(expansion + 2) * inverse_logarithm
                                             - std::log1p(-inverse_magnitude);
    const auto logarithmic_half_denormal =
      std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min())) - std::log(2.0L);
    const auto value_logarithm = value == 0
                                  ? -std::numeric_limits<long double>::infinity()
                                  : logarithmic_absolute(
                                      boost::multiprecision::abs(value).template convert_to<FarReal>());
    const auto rounds_to_zero = log_add_positive(value_logarithm, logarithmic_remaining_bound)
                                <= logarithmic_half_denormal;
    const auto relatively_converged = value != 0
                                      && logarithmic_remaining_bound
                                           <= std::log(static_cast<long double>(std::numeric_limits<double>::epsilon()))
                                                + value_logarithm;
    return std::pair{rounds_to_zero || relatively_converged, rounds_to_zero};
  };

  for (std::size_t expansion = 0; expansion < maximum_expansions; ++expansion) {
    const auto moment = exact_moment(expansion);
    const auto term = multiply(ExactComplex{moment.first, moment.second}, inverse_power);
    sum.real += term.real;
    sum.imaginary += term.imaginary;
    const auto real_status = component_converged(sum.real, output_real_logarithmic_bound, expansion);
    const auto imaginary_status = component_converged(sum.imaginary, output_imaginary_logarithmic_bound, expansion);
    if (real_status.first && imaginary_status.first)
      return std::complex<double>{real_status.second ? 0.0 : sum.real.convert_to<double>(),
                                  imaginary_status.second ? 0.0 : sum.imaginary.convert_to<double>()};
    inverse_power = multiply(inverse_power, inverse_root);
  }
  throw std::overflow_error("Exact far-field Cauchy-transform series did not converge.");
}

template <typename Argument>
inline auto normalized_root(const Argument argument, const double left, const double right) {
  using Real = decltype(std::abs(argument));
  const auto extended_left = static_cast<Real>(left);
  const auto extended_right = static_cast<Real>(right);
  const auto width = extended_right - extended_left;
  if (std::abs(argument - extended_left) <= std::abs(argument - extended_right))
    return (argument - extended_left) / width;
  return Argument{1.0} + (argument - extended_right) / width;
}

template <typename Scalar, typename Argument>
inline auto evaluate_polynomial(const std::vector<Scalar> &coefficients, const Argument argument) {
  using Real = decltype(std::abs(argument));
  using Result = std::complex<Real>;
  Result result{};
  for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient) {
    const auto value = polynomial_scalar_as_complex(*coefficient);
    result = result * argument + Result{static_cast<Real>(value.real()), static_cast<Real>(value.imag())};
  }
  return result;
}

template <typename Scalar, typename Argument>
inline auto quotient_integral(const std::vector<Scalar> &coefficients, const Argument root) {
  using Real = decltype(std::abs(root));
  using Result = std::complex<Real>;
  if (coefficients.size() == 1) return Result{};
  const auto last = polynomial_scalar_as_complex(coefficients.back());
  auto quotient = Result{static_cast<Real>(last.real()), static_cast<Real>(last.imag())};
  auto result = quotient / static_cast<Real>(coefficients.size() - 1);
  for (std::size_t power = coefficients.size() - 1; power > 1; --power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power - 1]);
    quotient = Result{static_cast<Real>(coefficient.real()), static_cast<Real>(coefficient.imag())} + root * quotient;
    result += quotient / static_cast<Real>(power - 1);
  }
  return result;
}

template <typename Scalar>
inline auto wide_coefficient(const Scalar coefficient) {
  const auto value = polynomial_scalar_as_complex(coefficient);
  return PiecewiseWideComplex{value.real(), value.imag()};
}

template <typename Scalar, typename Root>
inline auto wide_polynomial_value(const std::vector<Scalar> &coefficients, const Root &root) {
  PiecewiseWideComplex result;
  for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
    result = result * root + wide_coefficient(*coefficient);
  return result;
}

template <typename Scalar, typename Root>
inline auto wide_quotient_integral(const std::vector<Scalar> &coefficients, const Root &root) {
  if (coefficients.size() == 1) return PiecewiseWideComplex{};
  auto quotient = wide_coefficient(coefficients.back());
  auto result = quotient / PiecewiseWideFloat{coefficients.size() - 1};
  for (std::size_t power = coefficients.size() - 1; power > 1; --power) {
    quotient = wide_coefficient(coefficients[power - 1]) + root * quotient;
    result += quotient / PiecewiseWideFloat{power - 1};
  }
  return result;
}

template <typename Scalar>
inline auto wide_complex_interval_transform(const std::vector<Scalar> &coefficients,
                                            const std::complex<double> subtraction, const double left,
                                            const double right, const std::complex<double> argument) {
  const PiecewiseWideComplex wide_argument{argument.real(), argument.imag()};
  const PiecewiseWideComplex left_delta = wide_argument - PiecewiseWideFloat{left};
  const PiecewiseWideComplex right_delta = wide_argument - PiecewiseWideFloat{right};
  const auto width = PiecewiseWideFloat{right} - PiecewiseWideFloat{left};
  const auto root = boost::multiprecision::abs(left_delta) <= boost::multiprecision::abs(right_delta)
                      ? left_delta / width
                      : PiecewiseWideComplex{1} + right_delta / width;
  const auto logarithm = boost::multiprecision::log(left_delta) - boost::multiprecision::log(right_delta);
  const PiecewiseWideComplex wide_subtraction{subtraction.real(), subtraction.imag()};
  const PiecewiseWideComplex result = (wide_polynomial_value(coefficients, root) - wide_subtraction) * logarithm
                                      - wide_quotient_integral(coefficients, root);
  const PiecewiseWideFloat result_real = result.real();
  const PiecewiseWideFloat result_imaginary = result.imag();
  return std::complex<double>{result_real.convert_to<double>(), result_imaginary.convert_to<double>()};
}

template <typename Scalar>
inline auto wide_real_interval_transform(const std::vector<Scalar> &coefficients,
                                         const std::complex<double> subtraction, const double left,
                                         const double right, const double argument) {
  const auto wide_argument = PiecewiseWideFloat{argument};
  const auto left_delta = wide_argument - PiecewiseWideFloat{left};
  const auto right_delta = wide_argument - PiecewiseWideFloat{right};
  const auto width = PiecewiseWideFloat{right} - PiecewiseWideFloat{left};
  const auto root = boost::multiprecision::abs(left_delta) <= boost::multiprecision::abs(right_delta)
                      ? left_delta / width
                      : PiecewiseWideFloat{1} + right_delta / width;
  const auto logarithm = boost::multiprecision::log(boost::multiprecision::abs(left_delta))
                         - boost::multiprecision::log(boost::multiprecision::abs(right_delta));
  const PiecewiseWideComplex wide_subtraction{subtraction.real(), subtraction.imag()};
  const PiecewiseWideComplex result = (wide_polynomial_value(coefficients, root) - wide_subtraction) * logarithm
                                      - wide_quotient_integral(coefficients, root);
  const PiecewiseWideFloat result_real = result.real();
  const PiecewiseWideFloat result_imaginary = result.imag();
  return std::complex<double>{result_real.convert_to<double>(), result_imaginary.convert_to<double>()};
}

template <typename Moment, typename ExactMoment, typename ExactRoot>
inline auto evaluate_far_series(const std::complex<long double> inverse_root,
                                const std::optional<FarComplex> &inverse_root_extended,
                                const long double real_logarithmic_bound,
                                const long double imaginary_logarithmic_bound,
                                const std::size_t maximum_expansions, Moment moment, ExactMoment exact_moment,
                                ExactRoot exact_root, const char *failure_message) {
  FarSeriesAccumulator sum;
  const auto use_extended_root = inverse_root_extended
                                 && (inverse_root == std::complex<long double>{} || !native_extended_precision);
  if (use_extended_root) sum.force_wide(*inverse_root_extended);

  auto inverse_power = inverse_root;
  const auto inverse_logarithm = use_extended_root ? logarithmic_magnitude(*inverse_root_extended)
                                                    : std::log(std::abs(inverse_root));
  const auto inverse_magnitude = std::exp(inverse_logarithm);
  const auto inverse_is_real = use_extended_root ? inverse_root_extended->imaginary == 0
                                                  : inverse_root.imag() == 0.0L;
  const auto combined_bound = log_add_positive(real_logarithmic_bound, imaginary_logarithmic_bound);
  const auto output_real_bound = inverse_is_real ? real_logarithmic_bound : combined_bound;
  const auto output_imaginary_bound = inverse_is_real ? imaginary_logarithmic_bound : combined_bound;

  for (std::size_t expansion = 0; expansion < maximum_expansions; ++expansion) {
    sum.add_moment(moment(expansion), inverse_root, inverse_root_extended, inverse_power, expansion);
    if (sum.output_is_ill_conditioned())
      return exact_far_series(exact_root(), inverse_logarithm, inverse_magnitude, output_real_bound,
                              output_imaginary_bound, maximum_expansions, exact_moment);
    if (sum.converged(real_logarithmic_bound, imaginary_logarithmic_bound, output_real_bound,
                      output_imaginary_bound, expansion, inverse_logarithm, inverse_magnitude))
      return sum.value();
    const auto next_inverse_power = inverse_power * inverse_root;
    sum.advance_or_promote(inverse_root, next_inverse_power, expansion + 1);
    inverse_power = next_inverse_power;
  }
  throw std::overflow_error(failure_message);
}

template <typename Scalar>
inline std::complex<double> far_interval_transform(const std::vector<Scalar> &coefficients,
                                                   const std::complex<double> subtraction,
                                                   const std::complex<long double> inverse_root,
                                                   const std::optional<FarComplex> &inverse_root_extended,
                                                   const double left, const double right,
                                                   const std::complex<double> argument) {
  auto effective_inverse_root_extended = inverse_root_extended;
  if constexpr (!native_extended_precision) {
    if (!effective_inverse_root_extended)
      effective_inverse_root_extended.emplace(
        far_inverse_root_from_argument(right, left, argument.real(), argument.imag()));
  }
  const std::complex<long double> extended_subtraction{subtraction.real(), subtraction.imag()};
  auto real_logarithmic_bound = -std::numeric_limits<long double>::infinity();
  auto imaginary_logarithmic_bound = -std::numeric_limits<long double>::infinity();
  auto real_absolute_bound = 0.0L;
  auto imaginary_absolute_bound = 0.0L;
  auto bounds_are_finite = true;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
    const std::complex<long double> extended_coefficient{coefficient.real(), coefficient.imag()};
    const auto adjusted = extended_coefficient - (power == 0 ? extended_subtraction : std::complex<long double>{});
    real_absolute_bound += std::abs(adjusted.real()) / static_cast<long double>(power + 1);
    imaginary_absolute_bound += std::abs(adjusted.imag()) / static_cast<long double>(power + 1);
    bounds_are_finite = bounds_are_finite && std::isfinite(real_absolute_bound) && std::isfinite(imaginary_absolute_bound);
  }
  if (bounds_are_finite) {
    real_logarithmic_bound = log_absolute(real_absolute_bound);
    imaginary_logarithmic_bound = log_absolute(imaginary_absolute_bound);
  } else {
    for (std::size_t power = 0; power < coefficients.size(); ++power) {
      const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
      const auto denominator_logarithm = std::log(static_cast<long double>(power + 1));
      const auto real_logarithm = power == 0 ? log_absolute_difference(coefficient.real(), subtraction.real())
                                             : log_absolute(static_cast<long double>(coefficient.real()));
      const auto imaginary_logarithm = power == 0 ? log_absolute_difference(coefficient.imag(), subtraction.imag())
                                                  : log_absolute(static_cast<long double>(coefficient.imag()));
      real_logarithmic_bound = log_add_positive(real_logarithmic_bound, real_logarithm - denominator_logarithm);
      imaginary_logarithmic_bound = log_add_positive(imaginary_logarithmic_bound,
                                                      imaginary_logarithm - denominator_logarithm);
    }
  }
  if (real_logarithmic_bound == -std::numeric_limits<long double>::infinity()
      && imaginary_logarithmic_bound == -std::numeric_limits<long double>::infinity())
    return std::complex<double>{};

  if (coefficients.size() > std::numeric_limits<std::size_t>::max() - double_exponent_span)
    throw std::length_error("Piecewise-polynomial coefficient count is too large.");
  const auto maximum_expansions = coefficients.size() + double_exponent_span;
  auto moment = [&](const std::size_t expansion) {
    CompensatedLongComplexSum moment_sum;
    auto real_moment_absolute_bound = 0.0L;
    auto imaginary_moment_absolute_bound = 0.0L;
    for (std::size_t power = 0; power < coefficients.size(); ++power) {
      const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
      const std::complex<long double> extended_coefficient{coefficient.real(), coefficient.imag()};
      const auto contribution = (extended_coefficient - (power == 0 ? extended_subtraction : std::complex<long double>{}))
                                / static_cast<long double>(power + expansion + 1);
      moment_sum.add(contribution);
      real_moment_absolute_bound += std::abs(contribution.real());
      imaginary_moment_absolute_bound += std::abs(contribution.imag());
    }
    const auto moment_value = moment_sum.value();
    if constexpr (native_extended_precision) {
      return resolve_uncertain_moment(
        moment_value, log_absolute(real_moment_absolute_bound), log_absolute(imaginary_moment_absolute_bound),
        true,
        [&] { return high_precision_local_moment(coefficients, subtraction, expansion); },
        [&] { return exact_local_moment(coefficients, subtraction, expansion); });
    } else {
      return resolve_high_precision_moment(
        log_absolute(real_moment_absolute_bound), log_absolute(imaginary_moment_absolute_bound),
        [&] { return high_precision_local_moment(coefficients, subtraction, expansion); },
        [&] { return exact_local_moment(coefficients, subtraction, expansion); });
    }
  };
  auto exact_moment = [&](const std::size_t expansion) {
    return exact_local_moment(coefficients, subtraction, expansion);
  };
  auto exact_root = [&] {
    return exact_inverse_root_from_argument(right, left, argument.real(), argument.imag());
  };
  return evaluate_far_series(inverse_root, effective_inverse_root_extended, real_logarithmic_bound,
                             imaginary_logarithmic_bound, maximum_expansions, moment, exact_moment, exact_root,
                             "Far-interval Cauchy-transform series did not converge.");
}

template <typename Scalar>
inline std::optional<std::complex<double>> try_fast_real_far_interval_transform(
  const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
  const long double inverse_root) {
  // Returning no value restarts the interval with the hardened logarithmic path.
  if constexpr (!std::is_same_v<Scalar, double> || !native_extended_precision) {
    return std::nullopt;
  } else {
    if (subtraction.imag() != 0.0 || !std::isfinite(inverse_root) || inverse_root == 0.0L
        || !(std::abs(inverse_root) < 0.5L))
      return std::nullopt;
    constexpr auto maximum_coefficients = std::size_t{4};
    if (coefficients.size() > maximum_coefficients) return std::nullopt;
    static const auto reciprocals = [] {
      std::array<long double, double_exponent_span + 2 * maximum_coefficients> values{};
      for (std::size_t denominator = 1; denominator < values.size(); ++denominator)
        values[denominator] = 1.0L / static_cast<long double>(denominator);
      return values;
    }();

    const auto extended_subtraction = static_cast<long double>(subtraction.real());
    std::array<long double, maximum_coefficients> adjusted_coefficients{};
    auto moment_zero_bound = 0.0L;
    for (std::size_t power = 0; power < coefficients.size(); ++power) {
      adjusted_coefficients[power] = static_cast<long double>(coefficients[power])
                                     - (power == 0 ? extended_subtraction : 0.0L);
      moment_zero_bound += std::abs(adjusted_coefficients[power]) * reciprocals[power + 1];
    }
    if (moment_zero_bound == 0.0L) return std::complex<double>{};
    if (!std::isfinite(moment_zero_bound)) return std::nullopt;

    constexpr auto conditioning_threshold = 8.0L * std::numeric_limits<long double>::epsilon()
                                            / std::numeric_limits<double>::epsilon();
    auto well_conditioned = [](const long double value, const long double absolute_bound) {
      if (!std::isfinite(value) || !std::isfinite(absolute_bound)) return false;
      if (absolute_bound == 0.0L) return true;
      const auto ratio = std::abs(value) / absolute_bound;
      return std::isfinite(ratio) && ratio > conditioning_threshold;
    };
    auto add_compensated = [](const long double value, long double &sum, long double &correction) {
      const auto updated = sum + value;
      if (std::abs(sum) >= std::abs(value))
        correction += (sum - updated) + value;
      else
        correction += (value - updated) + sum;
      sum = updated;
    };

    auto sum = 0.0L;
    auto correction = 0.0L;
    auto absolute_sum = 0.0L;
    auto inverse_power = inverse_root;
    const auto inverse_magnitude = std::abs(inverse_root);
    auto inverse_magnitude_power = inverse_magnitude;
    if (coefficients.size() > std::numeric_limits<std::size_t>::max() - double_exponent_span)
      return std::nullopt;
    const auto maximum_expansions = coefficients.size() + double_exponent_span;

    for (std::size_t expansion = 0; expansion < maximum_expansions; ++expansion) {
      auto moment_sum = 0.0L;
      auto moment_correction = 0.0L;
      auto moment_absolute_bound = 0.0L;
      for (std::size_t power = 0; power < coefficients.size(); ++power) {
        const auto contribution = adjusted_coefficients[power] * reciprocals[power + expansion + 1];
        add_compensated(contribution, moment_sum, moment_correction);
        moment_absolute_bound += std::abs(contribution);
      }
      const auto moment = moment_sum + moment_correction;
      if (!well_conditioned(moment, moment_absolute_bound)) return std::nullopt;

      const auto term = moment * inverse_power;
      if (!std::isfinite(term) || (term == 0.0L && moment != 0.0L && inverse_power != 0.0L)) return std::nullopt;
      add_compensated(term, sum, correction);
      absolute_sum += std::abs(term);
      if (!std::isfinite(absolute_sum)) return std::nullopt;
      const auto value = sum + correction;
      if (!well_conditioned(value, absolute_sum)) return std::nullopt;

      const auto next_inverse_magnitude_power = inverse_magnitude_power * inverse_magnitude;
      if (!std::isfinite(next_inverse_magnitude_power)
          || (next_inverse_magnitude_power == 0.0L && inverse_magnitude_power != 0.0L))
        return std::nullopt;
      const auto remaining_numerator = moment_zero_bound * next_inverse_magnitude_power;
      if (!std::isfinite(remaining_numerator)
          || (remaining_numerator == 0.0L && moment_zero_bound != 0.0L
              && next_inverse_magnitude_power != 0.0L))
        return std::nullopt;
      const auto remaining_bound = remaining_numerator / (1.0L - inverse_magnitude);
      const auto target = std::numeric_limits<double>::epsilon() * std::abs(value);
      if (std::isfinite(remaining_bound) && target != 0.0L && remaining_bound <= target) {
        const auto magnitude = std::abs(value);
        if (magnitude > static_cast<long double>(std::numeric_limits<double>::max())
            || magnitude < static_cast<long double>(std::numeric_limits<double>::denorm_min()))
          return std::nullopt;
        const auto result = static_cast<double>(value);
        if (!std::isfinite(result)) return std::nullopt;
        return std::complex<double>{result, 0.0};
      }

      const auto next_inverse_power = inverse_power * inverse_root;
      if (!std::isfinite(next_inverse_power) || (next_inverse_power == 0.0L && inverse_power != 0.0L))
        return std::nullopt;
      inverse_power = next_inverse_power;
      inverse_magnitude_power = next_inverse_magnitude_power;
    }
    return std::nullopt;
  }
}

template <typename Scalar>
inline std::complex<double> complex_interval_transform(
  const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
  const std::complex<long double> root, const std::complex<long double> inverse_root,
  const std::optional<FarComplex> &inverse_root_extended, const std::complex<long double> left_delta,
  const std::complex<long double> right_delta, const double left, const double right,
  const std::complex<double> argument) {
  if constexpr (!native_extended_precision) {
    const auto wide_inverse_root = far_inverse_root_from_argument(right, left, argument.real(), argument.imag());
    if (logarithmic_magnitude(wide_inverse_root) < -std::log(2.0L))
      return far_interval_transform(coefficients, subtraction, inverse_root,
                                    std::optional<FarComplex>{wide_inverse_root}, left, right, argument);
    return wide_complex_interval_transform(coefficients, subtraction, left, right, argument);
  }
  if (std::abs(inverse_root) < 0.5)
    return far_interval_transform(coefficients, subtraction, inverse_root, inverse_root_extended, left, right, argument);
  const auto logarithm = std::log(left_delta) - std::log(right_delta);
  const auto extended_subtraction = std::complex<long double>{subtraction.real(), subtraction.imag()};
  const auto value_at_root = evaluate_polynomial(coefficients, root) - extended_subtraction;
  const auto result = value_at_root * logarithm - quotient_integral(coefficients, root);
  return {static_cast<double>(result.real()), static_cast<double>(result.imag())};
}

template <typename Scalar>
inline std::complex<double> real_interval_transform(
  const std::vector<Scalar> &coefficients, const std::complex<double> subtraction, const long double root,
  const long double inverse_root, const std::optional<FarComplex> &inverse_root_extended,
  const long double left_delta, const long double right_delta, const bool singular, const double left,
  const double right, const double argument) {
  if constexpr (!native_extended_precision) {
    if (singular) {
      const auto root_numerator = PiecewiseWideFloat{argument} - PiecewiseWideFloat{left};
      const auto width = PiecewiseWideFloat{right} - PiecewiseWideFloat{left};
      const PiecewiseWideComplex result = -wide_quotient_integral(coefficients, root_numerator / width);
      const PiecewiseWideFloat result_real = result.real();
      const PiecewiseWideFloat result_imaginary = result.imag();
      return {result_real.convert_to<double>(), result_imaginary.convert_to<double>()};
    }
    const auto wide_inverse_root = far_inverse_root_from_argument(right, left, argument, 0.0);
    if (logarithmic_magnitude(wide_inverse_root) < -std::log(2.0L))
      return far_interval_transform(coefficients, subtraction, std::complex<long double>{inverse_root, 0.0L},
                                    std::optional<FarComplex>{wide_inverse_root}, left, right, {argument, 0.0});
    return wide_real_interval_transform(coefficients, subtraction, left, right, argument);
  } else if (singular) {
    const auto result = -quotient_integral(coefficients, root);
    return {static_cast<double>(result.real()), static_cast<double>(result.imag())};
  }
  if (std::abs(inverse_root) < 0.5) {
    if (const auto fast = try_fast_real_far_interval_transform(coefficients, subtraction, inverse_root)) return *fast;
    return far_interval_transform(coefficients, subtraction, std::complex<long double>{inverse_root, 0.0L},
                                  inverse_root_extended, left, right, {argument, 0.0});
  }
  const auto logarithm = std::log(std::abs(left_delta)) - std::log(std::abs(right_delta));
  const auto extended_subtraction = std::complex<long double>{subtraction.real(), subtraction.imag()};
  const auto value_at_root = evaluate_polynomial(coefficients, root) - extended_subtraction;
  const auto result = value_at_root * logarithm - quotient_integral(coefficients, root);
  return {static_cast<double>(result.real()), static_cast<double>(result.imag())};
}

template <typename Scalar>
inline std::complex<double> global_far_transform(const PiecewisePolynomial<Scalar> &polynomial,
                                                 const std::complex<long double> inverse_root,
                                                 const std::optional<FarComplex> &inverse_root_extended,
                                                 const std::complex<double> argument) {
  const auto global_left = polynomial.lower_bound();
  const auto extended_global_width = static_cast<long double>(polynomial.upper_bound())
                                     - static_cast<long double>(global_left);
  const auto global_width_logarithm = native_extended_precision
                                        ? std::log(extended_global_width)
                                        : logarithmic_absolute(FarReal{polynomial.upper_bound()} - FarReal{global_left});
  auto effective_inverse_root_extended = inverse_root_extended;
  if constexpr (!native_extended_precision) {
    if (!effective_inverse_root_extended)
      effective_inverse_root_extended.emplace(far_inverse_root_from_argument(
        polynomial.upper_bound(), polynomial.lower_bound(), argument.real(), argument.imag()));
  }
  auto real_logarithmic_bound = -std::numeric_limits<long double>::infinity();
  auto imaginary_logarithmic_bound = -std::numeric_limits<long double>::infinity();
  for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
    const auto interval_width_logarithm = native_extended_precision
                                            ? std::log(static_cast<long double>(polynomial.knots()[interval + 1])
                                                       - static_cast<long double>(polynomial.knots()[interval]))
                                            : logarithmic_absolute(FarReal{polynomial.knots()[interval + 1]}
                                                                   - FarReal{polynomial.knots()[interval]});
    const auto width_fraction_logarithm = interval_width_logarithm - global_width_logarithm;
    for (std::size_t power = 0; power < polynomial.coefficients()[interval].size(); ++power) {
      const auto coefficient = polynomial_scalar_as_complex(polynomial.coefficients()[interval][power]);
      const auto factor_logarithm = width_fraction_logarithm - std::log(static_cast<long double>(power + 1));
      real_logarithmic_bound = log_add_positive(real_logarithmic_bound,
                                                log_absolute(static_cast<long double>(coefficient.real())) + factor_logarithm);
      imaginary_logarithmic_bound = log_add_positive(imaginary_logarithmic_bound,
                                                     log_absolute(static_cast<long double>(coefficient.imag())) + factor_logarithm);
    }
  }
  if (real_logarithmic_bound == -std::numeric_limits<long double>::infinity()
      && imaginary_logarithmic_bound == -std::numeric_limits<long double>::infinity())
    return std::complex<double>{};

  std::size_t coefficient_count = 0;
  for (const auto &coefficients : polynomial.coefficients()) {
    if (coefficient_count > std::numeric_limits<std::size_t>::max() - coefficients.size())
      throw std::length_error("Piecewise-polynomial coefficient count is too large.");
    coefficient_count += coefficients.size();
  }
  if (coefficient_count == std::numeric_limits<std::size_t>::max())
    throw std::length_error("Piecewise-polynomial coefficient count is too large.");
  const auto maximum_expansions = std::max<std::size_t>(double_exponent_span, coefficient_count + 1);
  auto moment = [&](const std::size_t expansion) {
    if constexpr (native_extended_precision) {
      const auto native_moment = native_global_moment(polynomial, expansion);
      return resolve_uncertain_moment(
        native_moment.value, log_absolute(native_moment.real_absolute_bound),
        log_absolute(native_moment.imaginary_absolute_bound),
        true,
        [&] { return high_precision_global_moment(polynomial, expansion); },
        [&] { return exact_global_moment(polynomial, expansion); });
    } else {
      return resolve_high_precision_moment(
        real_logarithmic_bound, imaginary_logarithmic_bound,
        [&] { return high_precision_global_moment(polynomial, expansion); },
        [&] { return exact_global_moment(polynomial, expansion); });
    }
  };
  auto exact_moment = [&](const std::size_t expansion) { return exact_global_moment(polynomial, expansion); };
  auto exact_root = [&] {
    return exact_inverse_root_from_argument(polynomial.upper_bound(), polynomial.lower_bound(), argument.real(),
                                            argument.imag());
  };
  return evaluate_far_series(inverse_root, effective_inverse_root_extended, real_logarithmic_bound,
                             imaginary_logarithmic_bound, maximum_expansions, moment, exact_moment, exact_root,
                             "Global far-field Cauchy-transform series did not converge.");
}

inline void require_finite_transform(const std::complex<double> value) {
  if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
    throw std::overflow_error("Piecewise-polynomial Cauchy transform is not finite.");
}

} // namespace detail

template <typename Scalar>
auto cauchy_transform(const PiecewisePolynomial<Scalar> &polynomial, const std::complex<double> argument) {
  if (!std::isfinite(argument.real()) || !std::isfinite(argument.imag()))
    throw std::invalid_argument("Cauchy-transform argument must be finite.");
  if (argument.imag() == 0.0)
    throw std::invalid_argument("Use the principal-value transform for a real Cauchy-transform argument.");

  const auto global_left = polynomial.lower_bound();
  const auto global_width = static_cast<long double>(polynomial.upper_bound()) - static_cast<long double>(global_left);
  const std::complex<long double> extended_argument{argument.real(), argument.imag()};
  const auto global_inverse_root = global_width / (extended_argument - static_cast<long double>(global_left));
  std::optional<detail::FarComplex> global_inverse_root_extended;
  if (global_inverse_root == std::complex<long double>{} || !detail::native_extended_precision)
    global_inverse_root_extended.emplace(detail::far_inverse_root_from_argument(
      polynomial.upper_bound(), global_left, argument.real(), argument.imag()));
  const auto global_far = detail::native_extended_precision
                            ? std::abs(global_inverse_root) < 0.5
                            : detail::logarithmic_magnitude(*global_inverse_root_extended) < -std::log(2.0L);
  if (global_far) {
    const auto result = detail::global_far_transform(polynomial, global_inverse_root, global_inverse_root_extended,
                                                     argument);
    detail::require_finite_transform(result);
    return result;
  }

  std::complex<double> subtraction{};
  if (polynomial.lower_bound() <= argument.real() && argument.real() <= polynomial.upper_bound())
    subtraction = polynomial_scalar_as_complex(polynomial.evaluate(argument.real()));

  CompensatedComplexSum sum;
  for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
    const auto left = polynomial.knots()[interval];
    const auto right = polynomial.knots()[interval + 1];
    const auto delta = extended_argument - static_cast<long double>(left);
    const auto root = detail::normalized_root(extended_argument, left, right);
    const auto interval_width = static_cast<long double>(right) - static_cast<long double>(left);
    const auto inverse_root = interval_width / (extended_argument - static_cast<long double>(left));
    std::optional<detail::FarComplex> inverse_root_extended;
    if (inverse_root == std::complex<long double>{})
      inverse_root_extended.emplace(detail::far_inverse_root_from_argument(right, left, argument.real(),
                                                                            argument.imag()));
    sum.add(detail::complex_interval_transform(polynomial.coefficients()[interval], subtraction, root, inverse_root,
                                               inverse_root_extended, delta,
                                               extended_argument - static_cast<long double>(right), left, right,
                                               argument));
  }

  if (subtraction != std::complex<double>{}) {
    const auto left = polynomial.lower_bound();
    const auto right = polynomial.upper_bound();
    const auto delta = extended_argument - static_cast<long double>(left);
    const auto root = detail::normalized_root(extended_argument, left, right);
    const auto interval_width = static_cast<long double>(right) - static_cast<long double>(left);
    const auto inverse_root = interval_width / (extended_argument - static_cast<long double>(left));
    std::optional<detail::FarComplex> inverse_root_extended;
    if (inverse_root == std::complex<long double>{})
      inverse_root_extended.emplace(detail::far_inverse_root_from_argument(right, left, argument.real(),
                                                                            argument.imag()));
    const std::vector<double> constant{1.0};
    sum.add(subtraction * detail::complex_interval_transform(constant, {}, root, inverse_root, inverse_root_extended,
                                                              delta, extended_argument - static_cast<long double>(right),
                                                              left, right, argument));
  }
  const auto result = sum.value();
  detail::require_finite_transform(result);
  return result;
}

template <typename Scalar>
auto cauchy_principal_value(const PiecewisePolynomial<Scalar> &polynomial, const double argument,
                            const CauchyEndpointPolicy endpoint_policy = CauchyEndpointPolicy::reject) {
  if (!std::isfinite(argument)) throw std::invalid_argument("Principal-value argument must be finite.");
  const auto at_endpoint = argument == polynomial.lower_bound() || argument == polynomial.upper_bound();
  if (at_endpoint && endpoint_policy == CauchyEndpointPolicy::reject)
    throw std::domain_error("The Cauchy transform is singular at the interpolation-domain endpoint.");
  const auto inside = polynomial.lower_bound() <= argument && argument <= polynomial.upper_bound();
  if (!inside) {
    const auto global_width = static_cast<long double>(polynomial.upper_bound())
                              - static_cast<long double>(polynomial.lower_bound());
    const auto global_inverse_root = global_width
                                     / (static_cast<long double>(argument) - static_cast<long double>(polynomial.lower_bound()));
    std::optional<detail::FarComplex> inverse_root_extended;
    if (global_inverse_root == 0.0L || !detail::native_extended_precision)
      inverse_root_extended.emplace(detail::far_inverse_root_from_argument(
        polynomial.upper_bound(), polynomial.lower_bound(), argument, 0.0));
    const auto global_far = detail::native_extended_precision
                              ? std::abs(global_inverse_root) < 0.5
                              : detail::logarithmic_magnitude(*inverse_root_extended) < -std::log(2.0L);
    if (global_far) {
      const auto result = detail::global_far_transform(polynomial,
                                                       std::complex<long double>{global_inverse_root, 0.0L},
                                                       inverse_root_extended, {argument, 0.0});
      detail::require_finite_transform(result);
      return result;
    }
  }
  if (inside && !at_endpoint) {
    const auto knot = std::lower_bound(polynomial.knots().begin(), polynomial.knots().end(), argument);
    if (knot != polynomial.knots().end() && *knot == argument) {
      const auto right_interval = static_cast<std::size_t>(std::distance(polynomial.knots().begin(), knot));
      const auto left_value = detail::evaluate_polynomial(polynomial.coefficients()[right_interval - 1], 1.0);
      const auto right_value = detail::evaluate_polynomial(polynomial.coefficients()[right_interval], 0.0);
      const auto scale = std::max({1.0, std::abs(left_value), std::abs(right_value)});
      if (std::abs(left_value - right_value) > 64.0 * std::numeric_limits<double>::epsilon() * scale)
        throw std::domain_error("The principal value diverges at a discontinuous polynomial knot.");
    }
  }
  const auto subtraction = inside ? polynomial_scalar_as_complex(polynomial.evaluate(argument)) : std::complex<double>{};

  CompensatedComplexSum sum;
  for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
    const auto left = polynomial.knots()[interval];
    const auto right = polynomial.knots()[interval + 1];
    const auto extended_argument = static_cast<long double>(argument);
    const auto delta = extended_argument - static_cast<long double>(left);
    const auto root = detail::normalized_root(extended_argument, left, right);
    const auto interval_width = static_cast<long double>(right) - static_cast<long double>(left);
    const auto inverse_root = interval_width / (static_cast<long double>(argument) - static_cast<long double>(left));
    std::optional<detail::FarComplex> inverse_root_extended;
    if (inverse_root == 0.0L)
      inverse_root_extended.emplace(detail::far_inverse_root_from_argument(right, left, argument, 0.0));
    const auto singular = inside && left <= argument && argument <= right;
    sum.add(detail::real_interval_transform(polynomial.coefficients()[interval], subtraction, root, inverse_root,
                                             inverse_root_extended, delta,
                                             extended_argument - static_cast<long double>(right), singular, left, right,
                                             argument));
  }

  if (inside && !at_endpoint && subtraction != std::complex<double>{}) {
    const auto logarithm = std::log(std::abs(argument - polynomial.lower_bound()))
                           - std::log(std::abs(argument - polynomial.upper_bound()));
    sum.add(subtraction * logarithm);
  }
  const auto result = sum.value();
  detail::require_finite_transform(result);
  return result;
}

} // namespace NRG::Tools

#endif
