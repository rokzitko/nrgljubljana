#ifndef _tools_common_piecewise_polynomial_hpp_
#define _tools_common_piecewise_polynomial_hpp_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace NRG::Tools {

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

template <typename Real>
class CompensatedComplexSumType {
  Real real_sum = 0.0;
  Real real_correction = 0.0;
  Real imaginary_sum = 0.0;
  Real imaginary_correction = 0.0;

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
  }

  auto value() const { return std::complex<Real>{real_sum + real_correction, imaginary_sum + imaginary_correction}; }
};

using CompensatedComplexSum = CompensatedComplexSumType<double>;
using CompensatedLongComplexSum = CompensatedComplexSumType<long double>;

template <typename Scalar>
class PiecewisePolynomial {
  static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, std::complex<double>>);

  std::vector<double> knots_;
  std::vector<std::vector<Scalar>> coefficients_;

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

  auto evaluate(const double x) const {
    const auto index = interval_index(x);
    const auto width = knots_[index + 1] - knots_[index];
    const auto u = (x - knots_[index]) / width;
    Scalar result{};
    for (auto coefficient = coefficients_[index].rbegin(); coefficient != coefficients_[index].rend(); ++coefficient)
      result = result * u + *coefficient;
    return result;
  }

  auto integral() const {
    CompensatedComplexSum sum;
    for (std::size_t interval = 0; interval < coefficients_.size(); ++interval) {
      const auto width = knots_[interval + 1] - knots_[interval];
      std::complex<double> value{};
      for (std::size_t power = 0; power < coefficients_[interval].size(); ++power)
        value += polynomial_scalar_as_complex(coefficients_[interval][power]) / static_cast<double>(power + 1);
      sum.add(width * value);
    }
    if constexpr (std::is_same_v<Scalar, double>)
      return sum.value().real();
    else
      return sum.value();
  }

  auto multiply_by_monomial(const std::size_t exponent) const {
    if (exponent == 0) return *this;
    std::vector<std::vector<Scalar>> weighted;
    weighted.reserve(coefficients_.size());
    for (std::size_t interval = 0; interval < coefficients_.size(); ++interval) {
      const auto left = knots_[interval];
      const auto width = knots_[interval + 1] - left;
      std::vector<double> monomial(exponent + 1, 0.0);
      if (left == 0.0) {
        monomial[exponent] = std::pow(width, static_cast<double>(exponent));
      } else if (std::abs(left) >= width) {
        monomial[0] = std::pow(left, static_cast<double>(exponent));
        for (std::size_t power = 0; power < exponent; ++power)
          monomial[power + 1] = monomial[power] * static_cast<double>(exponent - power)
                                / static_cast<double>(power + 1) * width / left;
      } else {
        monomial[exponent] = std::pow(width, static_cast<double>(exponent));
        for (std::size_t power = exponent; power > 0; --power)
          monomial[power - 1] = monomial[power] * static_cast<double>(power)
                                / static_cast<double>(exponent - power + 1) * left / width;
      }

      std::vector<Scalar> product(coefficients_[interval].size() + exponent, Scalar{});
      for (std::size_t first = 0; first < coefficients_[interval].size(); ++first)
        for (std::size_t second = 0; second < monomial.size(); ++second)
          product[first + second] += coefficients_[interval][first] * monomial[second];
      for (const auto coefficient : product)
        if (!polynomial_scalar_is_finite(coefficient))
          throw std::overflow_error("Energy weighting produced a nonfinite polynomial coefficient.");
      weighted.push_back(std::move(product));
    }
    return PiecewisePolynomial<Scalar>{knots_, std::move(weighted)};
  }
};

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

template <typename Argument>
inline auto normalized_root(const Argument argument, const double left, const double right) {
  const auto width = right - left;
  if (std::abs(argument - left) <= std::abs(argument - right)) return (argument - left) / width;
  return Argument{1.0} + (argument - right) / width;
}

template <typename Scalar, typename Argument>
inline auto evaluate_polynomial(const std::vector<Scalar> &coefficients, const Argument argument) {
  using Result = decltype(polynomial_scalar_as_complex(Scalar{}) * Argument{});
  Result result{};
  for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
    result = result * argument + polynomial_scalar_as_complex(*coefficient);
  return result;
}

template <typename Scalar, typename Argument>
inline auto quotient_integral(const std::vector<Scalar> &coefficients, const Argument root) {
  if (coefficients.size() == 1) return std::complex<double>{};
  auto quotient = polynomial_scalar_as_complex(coefficients.back());
  auto result = quotient / static_cast<double>(coefficients.size() - 1);
  for (std::size_t power = coefficients.size() - 1; power > 1; --power) {
    quotient = polynomial_scalar_as_complex(coefficients[power - 1]) + root * quotient;
    result += quotient / static_cast<double>(power - 1);
  }
  return result;
}

template <typename Scalar>
inline auto far_interval_transform(const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
                                   const std::complex<long double> inverse_root) {
  const std::complex<long double> extended_subtraction{subtraction.real(), subtraction.imag()};
  long double absolute_moment_bound = 0.0L;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
    const std::complex<long double> extended_coefficient{coefficient.real(), coefficient.imag()};
    const auto adjusted = extended_coefficient - (power == 0 ? extended_subtraction : std::complex<long double>{});
    absolute_moment_bound += std::abs(adjusted) / static_cast<long double>(power + 1);
  }
  if (absolute_moment_bound == 0.0L) return std::complex<double>{};

  CompensatedLongComplexSum sum;
  auto inverse_power = inverse_root;
  auto inverse_magnitude_power = std::abs(inverse_root);
  const auto inverse_magnitude = inverse_magnitude_power;
  const auto maximum_expansions = coefficients.size() + static_cast<std::size_t>(std::numeric_limits<long double>::max_exponent);
  for (std::size_t expansion = 0; expansion < maximum_expansions; ++expansion) {
    std::complex<long double> moment{};
    for (std::size_t power = 0; power < coefficients.size(); ++power) {
      const auto coefficient = polynomial_scalar_as_complex(coefficients[power]);
      const std::complex<long double> extended_coefficient{coefficient.real(), coefficient.imag()};
      moment += (extended_coefficient - (power == 0 ? extended_subtraction : std::complex<long double>{}))
                / static_cast<long double>(power + expansion + 1);
    }
    const auto term = moment * inverse_power;
    sum.add(term);
    const auto next_inverse_magnitude_power = inverse_magnitude_power * inverse_magnitude;
    const auto remaining_bound = absolute_moment_bound * next_inverse_magnitude_power / (1.0L - inverse_magnitude);
    if (std::abs(sum.value()) != 0.0
        && remaining_bound <= static_cast<long double>(std::numeric_limits<double>::epsilon())
                              * std::abs(sum.value()))
      return std::complex<double>{static_cast<double>(sum.value().real()), static_cast<double>(sum.value().imag())};
    inverse_power *= inverse_root;
    inverse_magnitude_power = next_inverse_magnitude_power;
  }
  throw std::overflow_error("Far-interval Cauchy-transform series did not converge.");
}

template <typename Scalar>
inline auto complex_interval_transform(const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
                                       const std::complex<double> root, const std::complex<long double> inverse_root,
                                       const std::complex<double> left_delta, const std::complex<double> right_delta) {
  if (std::abs(inverse_root) < 0.5) return far_interval_transform(coefficients, subtraction, inverse_root);
  const auto logarithm = std::log(left_delta) - std::log(right_delta);
  const auto value_at_root = evaluate_polynomial(coefficients, root) - subtraction;
  return value_at_root * logarithm - quotient_integral(coefficients, root);
}

template <typename Scalar>
inline auto real_interval_transform(const std::vector<Scalar> &coefficients, const std::complex<double> subtraction,
                                    const double root, const long double inverse_root, const double left_delta,
                                    const double right_delta, const bool singular) {
  if (singular) return -quotient_integral(coefficients, root);
  if (std::abs(inverse_root) < 0.5)
    return far_interval_transform(coefficients, subtraction, std::complex<long double>{inverse_root, 0.0L});
  const auto logarithm = std::log(std::abs(left_delta)) - std::log(std::abs(right_delta));
  const auto value_at_root = evaluate_polynomial(coefficients, root) - subtraction;
  return value_at_root * logarithm - quotient_integral(coefficients, root);
}

template <typename Scalar>
inline auto global_far_transform(const PiecewisePolynomial<Scalar> &polynomial,
                                 const std::complex<long double> inverse_root) {
  const auto global_left = polynomial.lower_bound();
  const auto global_width = polynomial.upper_bound() - global_left;
  const auto extended_global_width = static_cast<long double>(polynomial.upper_bound()) - static_cast<long double>(global_left);
  std::vector<std::vector<double>> global_powers(polynomial.interval_count(), std::vector<double>{1.0});
  long double absolute_moment_bound = 0.0L;
  for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
    const auto width_fraction = (static_cast<long double>(polynomial.knots()[interval + 1])
                                 - static_cast<long double>(polynomial.knots()[interval])) / extended_global_width;
    for (std::size_t power = 0; power < polynomial.coefficients()[interval].size(); ++power)
      absolute_moment_bound += width_fraction
                               * static_cast<long double>(std::abs(polynomial_scalar_as_complex(polynomial.coefficients()[interval][power])))
                               / static_cast<long double>(power + 1);
  }
  if (absolute_moment_bound == 0.0L) return std::complex<double>{};

  CompensatedLongComplexSum sum;
  auto inverse_power = inverse_root;
  auto inverse_magnitude_power = static_cast<long double>(std::abs(inverse_root));
  const auto inverse_magnitude = inverse_magnitude_power;
  for (std::size_t expansion = 0; expansion < 128; ++expansion) {
    CompensatedComplexSum moment;
    for (std::size_t interval = 0; interval < polynomial.interval_count(); ++interval) {
      const auto left_fraction = (polynomial.knots()[interval] - global_left) / global_width;
      const auto width_fraction = (polynomial.knots()[interval + 1] - polynomial.knots()[interval]) / global_width;
      const auto extended_width_fraction = (static_cast<long double>(polynomial.knots()[interval + 1])
                                            - static_cast<long double>(polynomial.knots()[interval])) / extended_global_width;
      const auto &power_coefficients = global_powers[interval];
      for (std::size_t local_power = 0; local_power < polynomial.coefficients()[interval].size(); ++local_power) {
        const auto coefficient = polynomial_scalar_as_complex(polynomial.coefficients()[interval][local_power]);
        for (std::size_t global_power = 0; global_power < power_coefficients.size(); ++global_power) {
          if (width_fraction != 0.0) {
            moment.add(width_fraction * coefficient * power_coefficients[global_power]
                       / static_cast<double>(local_power + global_power + 1));
          } else {
            const std::complex<long double> extended_coefficient{coefficient.real(), coefficient.imag()};
            const auto contribution = extended_width_fraction * extended_coefficient
                                      * static_cast<long double>(power_coefficients[global_power])
                                      / static_cast<long double>(local_power + global_power + 1);
            moment.add({static_cast<double>(contribution.real()), static_cast<double>(contribution.imag())});
          }
        }
      }

      std::vector<double> next_power(power_coefficients.size() + 1, 0.0);
      for (std::size_t power = 0; power < power_coefficients.size(); ++power) {
        next_power[power] += left_fraction * power_coefficients[power];
        next_power[power + 1] += width_fraction * power_coefficients[power];
      }
      global_powers[interval] = std::move(next_power);
    }

    const auto moment_value = moment.value();
    sum.add(std::complex<long double>{moment_value.real(), moment_value.imag()} * inverse_power);
    const auto next_inverse_magnitude_power = inverse_magnitude_power * inverse_magnitude;
    const auto remaining_bound = absolute_moment_bound * next_inverse_magnitude_power / (1.0L - inverse_magnitude);
    if (std::abs(sum.value()) != 0.0
        && remaining_bound <= static_cast<long double>(std::numeric_limits<double>::epsilon())
                              * std::abs(sum.value()))
      return std::complex<double>{static_cast<double>(sum.value().real()), static_cast<double>(sum.value().imag())};
    inverse_power *= inverse_root;
    inverse_magnitude_power = next_inverse_magnitude_power;
  }
  throw std::overflow_error("Global far-field Cauchy-transform series did not converge.");
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
  const auto global_width = polynomial.upper_bound() - global_left;
  const std::complex<long double> extended_argument{argument.real(), argument.imag()};
  const auto global_inverse_root = static_cast<long double>(global_width)
                                   / (extended_argument - static_cast<long double>(global_left));
  if (std::abs(global_inverse_root) < 0.5) {
    const auto result = detail::global_far_transform(polynomial, global_inverse_root);
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
    const auto delta = argument - left;
    const auto root = detail::normalized_root(argument, left, right);
    const auto inverse_root = (static_cast<long double>(right) - static_cast<long double>(left))
                              / (extended_argument - static_cast<long double>(left));
    sum.add(detail::complex_interval_transform(polynomial.coefficients()[interval], subtraction, root, inverse_root,
                                               delta, argument - right));
  }

  if (subtraction != std::complex<double>{}) {
    const auto left = polynomial.lower_bound();
    const auto right = polynomial.upper_bound();
    const auto delta = argument - left;
    const auto root = detail::normalized_root(argument, left, right);
    const auto inverse_root = (static_cast<long double>(right) - static_cast<long double>(left))
                              / (extended_argument - static_cast<long double>(left));
    const std::vector<double> constant{1.0};
    sum.add(subtraction * detail::complex_interval_transform(constant, {}, root, inverse_root, delta, argument - right));
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
    if (std::abs(global_inverse_root) < 0.5) {
      const auto result = detail::global_far_transform(polynomial, std::complex<long double>{global_inverse_root, 0.0L});
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
    const auto delta = argument - left;
    const auto root = detail::normalized_root(argument, left, right);
    const auto inverse_root = (static_cast<long double>(right) - static_cast<long double>(left))
                              / (static_cast<long double>(argument) - static_cast<long double>(left));
    const auto singular = inside && left <= argument && argument <= right;
    sum.add(detail::real_interval_transform(polynomial.coefficients()[interval], subtraction, root, inverse_root,
                                            delta, argument - right, singular));
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
