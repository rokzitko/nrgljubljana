#ifndef _tools_common_tabulated_density_hpp_
#define _tools_common_tabulated_density_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gsl_config.hpp"
#include "gsl_piecewise_polynomial.hpp"
#include "linint.hpp"

namespace NRG::Tools {

inline auto parse_density_interpolation_method(const std::string_view value) {
  if (value == "linear") return InterpolationMethod::linear;
  if (value == "steffen") return InterpolationMethod::steffen;
  throw std::invalid_argument("density_interpolation must be either 'linear' or 'steffen'.");
}

using DensityTable = std::vector<std::pair<double, double>>;

struct DensityInterpolationError {
  [[noreturn]] void operator()(const std::string &message) const { throw std::runtime_error(message); }
};

class TabulatedDensity {
  using LinearValue = LinIntBase<DensityTable, DensityInterpolationError>;

  bool initialized_{};
  InterpolationMethod method_ = InterpolationMethod::linear;
  DensityTable samples_;
  LinearValue linear_value_;
  std::vector<double> linear_prefix_integrals_;
  std::size_t linear_cumulative_lower_interval_ = std::numeric_limits<std::size_t>::max();
  std::size_t linear_cumulative_upper_interval_ = std::numeric_limits<std::size_t>::max();
  std::optional<PiecewisePolynomial<double>> polynomial_;
  std::vector<long double> prefix_integrals_;
  std::vector<long double> interval_mass_tree_;
  std::size_t interval_mass_tree_base_{};
  std::size_t integral_lower_interval_ = std::numeric_limits<std::size_t>::max();
  std::size_t integral_upper_interval_ = std::numeric_limits<std::size_t>::max();
  std::size_t tracked_interval_ = std::numeric_limits<std::size_t>::max();
  bool interval_changed_{};

  void require_initialized() const {
    if (!initialized_) throw std::logic_error("Tabulated density has not been initialized.");
  }

  static void validate_method(const InterpolationMethod method) {
    if (method != InterpolationMethod::linear && method != InterpolationMethod::steffen)
      throw std::invalid_argument("Density interpolation supports only linear and steffen methods.");
  }

  static void validate_samples(const DensityTable &samples, const InterpolationMethod method) {
    const auto minimum_size = interpolation_minimum_size(method);
    if (samples.size() < minimum_size)
      throw std::invalid_argument(std::string(interpolation_method_name(method)) + " density interpolation requires at least "
                                  + std::to_string(minimum_size) + " points.");
    for (std::size_t index = 0; index < samples.size(); ++index) {
      const auto [energy, value] = samples[index];
      if (!std::isfinite(energy) || !std::isfinite(value))
        throw std::invalid_argument("Density interpolation data must be finite.");
      if (value < 0.0) throw std::invalid_argument("Density interpolation values must be nonnegative.");
      if (index != 0 && !(samples[index - 1].first < energy))
        throw std::invalid_argument("Density interpolation energies must be strictly increasing.");
    }
    if (!std::isfinite(samples.back().first - samples.front().first))
      throw std::invalid_argument("Density interpolation domain width must be finite.");
  }

  auto interval_index(const double x) const {
    if (x == samples_.back().first) return samples_.size() - 2;
    const auto upper = std::upper_bound(samples_.begin(), samples_.end(), x,
                                        [](const double value, const auto &sample) { return value < sample.first; });
    return static_cast<std::size_t>(std::distance(samples_.begin(), upper) - 1);
  }

  auto cached_interval_index(const double x, std::size_t &cached) const {
    const auto interval_count = samples_.size() - 1;
    if (x == samples_.back().first) return interval_count - 1;
    if (cached >= interval_count) {
      cached = interval_index(x);
      return cached;
    }
    while (cached + 1 < interval_count && x >= samples_[cached + 1].first) ++cached;
    while (cached > 0 && x < samples_[cached].first) --cached;
    return cached;
  }

  auto legacy_linear_interval(const double x, std::size_t &cached) const {
    const auto interval_count = samples_.size() - 1;
    if (cached < interval_count && samples_[cached].first <= x && x < samples_[cached + 1].first) return cached;
    if (cached >= interval_count) {
      for (std::size_t interval = 0; interval < interval_count; ++interval) {
        if (samples_[interval].first <= x && x <= samples_[interval + 1].first) {
          cached = interval;
          return cached;
        }
      }
    } else if (x >= samples_[cached + 1].first) {
      for (auto interval = cached + 1; interval < interval_count; ++interval) {
        if (samples_[interval].first <= x && x <= samples_[interval + 1].first) {
          cached = interval;
          return cached;
        }
      }
    } else {
      for (auto interval = cached; interval-- > 0;) {
        if (samples_[interval].first <= x && x <= samples_[interval + 1].first) {
          cached = interval;
          return cached;
        }
      }
    }
    throw std::logic_error("Failed to locate a linear density interval.");
  }

  void initialize_legacy_linear_cumulative() {
    linear_prefix_integrals_.assign(samples_.size(), 0.0);
    double sum = 0.0;
    for (std::size_t index = 1; index < samples_.size(); ++index) {
      const auto [x0, y0] = samples_[index - 1];
      const auto [x1, y1] = samples_[index];
      const auto dx = x1 - x0;
      const auto yavg = (y1 + y0) / 2.0;
      sum += dx * yavg;
      linear_prefix_integrals_[index] = sum;
    }
  }

  auto legacy_linear_cumulative(const double x, std::size_t &cached) const {
    if (x <= samples_.front().first)
      return linear_prefix_integrals_.front() + samples_.front().second * (x - samples_.front().first);
    if (x >= samples_.back().first)
      return linear_prefix_integrals_.back() + samples_.back().second * (x - samples_.back().first);
    const auto interval = legacy_linear_interval(x, cached);
    const auto x0 = samples_[interval].first;
    const auto f0 = samples_[interval].second;
    const auto deriv = (samples_[interval + 1].second - f0) / (samples_[interval + 1].first - x0);
    const auto dx = x - x0;
    return linear_prefix_integrals_[interval] + dx * f0 + dx * dx * deriv / 2.0;
  }

  auto local_integral(const std::size_t interval, const double lower, const double upper) const {
    const auto left = samples_[interval].first;
    const auto right = samples_[interval + 1].first;
    if (method_ == InterpolationMethod::linear) {
      const auto lower_u = detail::normalized_interval_coordinate(lower, left, right);
      const auto upper_u = detail::normalized_interval_coordinate(upper, left, right);
      const auto value_left = static_cast<long double>(samples_[interval].second);
      const auto delta = static_cast<long double>(samples_[interval + 1].second) - value_left;
      const auto lower_value = value_left + delta * lower_u;
      const auto upper_value = value_left + delta * upper_u;
      const auto span = static_cast<long double>(upper) - static_cast<long double>(lower);
      return span * (lower_value + upper_value) / 2.0L;
    }

    const auto lower_u = detail::normalized_interval_coordinate(lower, left, right);
    const auto upper_u = detail::normalized_interval_coordinate(upper, left, right);
    const auto span = static_cast<long double>(upper) - static_cast<long double>(lower);
    return detail::normalized_polynomial_integral(polynomial_->coefficients()[interval], lower_u, upper_u, span)
      .real();
  }

  void initialize_polynomial() {
    std::vector<double> knots;
    std::vector<double> values;
    knots.reserve(samples_.size());
    values.reserve(samples_.size());
    for (const auto &[energy, value] : samples_) {
      knots.push_back(energy);
      values.push_back(value);
    }
    polynomial_.emplace(make_gsl_piecewise_polynomial(knots, values, method_));
  }

  void initialize_prefix_integrals() {
    prefix_integrals_.assign(samples_.size(), 0.0L);
    interval_mass_tree_base_ = 1;
    const auto interval_count = samples_.size() - 1;
    while (interval_mass_tree_base_ < interval_count) {
      if (interval_mass_tree_base_ > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("Density interpolation grid is too large.");
      interval_mass_tree_base_ *= 2;
    }
    if (interval_mass_tree_base_ > std::numeric_limits<std::size_t>::max() / 2)
      throw std::length_error("Density interpolation grid is too large.");
    interval_mass_tree_.assign(2 * interval_mass_tree_base_, 0.0L);
    long double sum = 0.0L;
    long double correction = 0.0L;
    for (std::size_t interval = 0; interval + 1 < samples_.size(); ++interval) {
      auto mass = local_integral(interval, samples_[interval].first, samples_[interval + 1].first);
      const auto width = static_cast<long double>(samples_[interval + 1].first)
                         - static_cast<long double>(samples_[interval].first);
      long double scale = 0.0L;
      if (method_ == InterpolationMethod::linear) {
        scale = (std::abs(static_cast<long double>(samples_[interval].second))
                 + std::abs(static_cast<long double>(samples_[interval + 1].second)))
                / 2.0L;
      } else {
        const auto &coefficients = polynomial_->coefficients()[interval];
        for (std::size_t degree = 0; degree < coefficients.size(); ++degree)
          scale += std::abs(static_cast<long double>(coefficients[degree])) / static_cast<long double>(degree + 1);
      }
      scale *= width;
      const auto tolerance = 64.0L * std::numeric_limits<long double>::epsilon() * scale;
      if (!std::isfinite(mass) || mass < -tolerance)
        throw std::runtime_error("Density interpolation produced an invalid interval weight.");
      if (mass < 0.0L) mass = 0.0L;
      const auto corrected = mass - correction;
      const auto next = sum + corrected;
      correction = (next - sum) - corrected;
      sum = next;
      prefix_integrals_[interval + 1] = sum;
      interval_mass_tree_[interval_mass_tree_base_ + interval] = mass;
    }
    for (auto node = interval_mass_tree_base_; node-- > 1;)
      interval_mass_tree_[node] = interval_mass_tree_[2 * node] + interval_mass_tree_[2 * node + 1];
    if (!std::isfinite(static_cast<double>(sum)))
      throw std::runtime_error("Integrated density is not finite.");
  }

  auto full_interval_mass(const std::size_t begin, const std::size_t end) const {
    auto left = interval_mass_tree_base_ + begin;
    auto right = interval_mass_tree_base_ + end;
    long double sum = 0.0L;
    long double correction = 0.0L;
    const auto add = [&sum, &correction](const long double value) {
      const auto corrected = value - correction;
      const auto next = sum + corrected;
      correction = (next - sum) - corrected;
      sum = next;
    };
    while (left < right) {
      if (left % 2 != 0) add(interval_mass_tree_[left++]);
      if (right % 2 != 0) add(interval_mass_tree_[--right]);
      left /= 2;
      right /= 2;
    }
    return sum;
  }

  auto cumulative_wide(const double x) const {
    if (x <= samples_.front().first)
      return static_cast<long double>(samples_.front().second)
             * (static_cast<long double>(x) - static_cast<long double>(samples_.front().first));
    if (x >= samples_.back().first)
      return prefix_integrals_.back() + static_cast<long double>(samples_.back().second)
                                              * (static_cast<long double>(x)
                                                 - static_cast<long double>(samples_.back().first));
    const auto interval = interval_index(x);
    return prefix_integrals_[interval] + local_integral(interval, samples_[interval].first, x);
  }

  auto integral_wide(const double lower, const double upper) {
    long double sum = 0.0L;
    long double correction = 0.0L;
    const auto add = [&sum, &correction](const long double value) {
      const auto corrected = value - correction;
      const auto next = sum + corrected;
      correction = (next - sum) - corrected;
      sum = next;
    };

    const auto domain_lower = samples_.front().first;
    const auto domain_upper = samples_.back().first;
    if (lower < domain_lower) {
      const auto tail_upper = std::min(upper, domain_lower);
      add(static_cast<long double>(samples_.front().second)
          * (static_cast<long double>(tail_upper) - static_cast<long double>(lower)));
    }

    const auto inside_lower = std::max(lower, domain_lower);
    const auto inside_upper = std::min(upper, domain_upper);
    if (inside_lower < inside_upper) {
      const auto first = cached_interval_index(inside_lower, integral_lower_interval_);
      const auto last = cached_interval_index(std::nextafter(inside_upper, inside_lower), integral_upper_interval_);
      if (first == last) {
        add(local_integral(first, inside_lower, inside_upper));
      } else {
        add(local_integral(first, inside_lower, samples_[first + 1].first));
        if (first + 1 < last) {
          if (method_ == InterpolationMethod::linear) {
            for (auto interval = first + 1; interval < last; ++interval)
              add(local_integral(interval, samples_[interval].first, samples_[interval + 1].first));
          } else {
            add(full_interval_mass(first + 1, last));
          }
        }
        add(local_integral(last, samples_[last].first, inside_upper));
      }
    }

    if (upper > domain_upper) {
      const auto tail_lower = std::max(lower, domain_upper);
      add(static_cast<long double>(samples_.back().second)
          * (static_cast<long double>(upper) - static_cast<long double>(tail_lower)));
    }
    return sum;
  }

  static auto checked_weight(const long double value) {
    const auto result = static_cast<double>(value);
    if (!std::isfinite(result)) throw std::runtime_error("Integrated density weight is not finite.");
    const auto tolerance = 64.0L * std::numeric_limits<long double>::epsilon()
                           * std::max(1.0L, std::abs(value));
    if (value < -tolerance) throw std::runtime_error("Integrated density weight is negative.");
    return value < 0.0L ? 0.0 : result;
  }

 public:
  TabulatedDensity() = default;

  explicit TabulatedDensity(const DensityTable &samples,
                            const InterpolationMethod method = InterpolationMethod::linear)
      : method_{method}, samples_{samples} {
    validate_method(method_);
    validate_samples(samples_, method_);
    if (method_ == InterpolationMethod::linear) {
      linear_value_ = LinearValue(samples_);
      initialize_legacy_linear_cumulative();
    } else {
      initialize_polynomial();
      initialize_prefix_integrals();
    }
    initialized_ = true;
  }

  auto interpolation_method() const {
    require_initialized();
    return method_;
  }

  double operator()(const double x) {
    require_initialized();
    if (!std::isfinite(x)) throw std::invalid_argument("Density interpolation argument must be finite.");
    if (method_ == InterpolationMethod::linear) return linear_value_(x);

    if (x <= samples_.front().first) return samples_.front().second;
    if (x >= samples_.back().first) return samples_.back().second;
    const auto interval = interval_index(x);
    if (tracked_interval_ != interval) {
      tracked_interval_ = interval;
      interval_changed_ = true;
    }
    const auto left = samples_[interval].first;
    const auto right = samples_[interval + 1].first;
    const auto u = detail::normalized_interval_coordinate(x, left, right);
    long double result = 0.0L;
    long double absolute_bound = 0.0L;
    for (auto coefficient = polynomial_->coefficients()[interval].rbegin();
         coefficient != polynomial_->coefficients()[interval].rend(); ++coefficient) {
      result = result * u + static_cast<long double>(*coefficient);
      absolute_bound = absolute_bound * std::abs(u) + std::abs(static_cast<long double>(*coefficient));
    }
    const auto negative_tolerance = 64.0L * std::numeric_limits<double>::epsilon() * absolute_bound;
    if (!std::isfinite(result) || result < -negative_tolerance)
      throw std::runtime_error("Density interpolation produced a negative or non-finite value.");
    if (result < 0.0L) return 0.0;
    const auto narrowed = static_cast<double>(result);
    if (!std::isfinite(narrowed)) throw std::runtime_error("Density interpolation value is not finite.");
    return narrowed;
  }

  double cumulative(const double x) {
    require_initialized();
    if (!std::isfinite(x)) throw std::invalid_argument("Cumulative-density argument must be finite.");
    if (method_ == InterpolationMethod::linear)
      return legacy_linear_cumulative(x, linear_cumulative_lower_interval_);
    const auto result = static_cast<double>(cumulative_wide(x));
    if (!std::isfinite(result)) throw std::runtime_error("Cumulative density is not finite.");
    return result;
  }

  double integral(const double lower, const double upper) {
    require_initialized();
    if (!std::isfinite(lower) || !std::isfinite(upper))
      throw std::invalid_argument("Density integration bounds must be finite.");
    if (lower > upper) throw std::invalid_argument("Density integration bounds must be ordered.");
    if (lower == upper) return 0.0;
    if (method_ != InterpolationMethod::linear) return checked_weight(integral_wide(lower, upper));

    const auto upper_cumulative = legacy_linear_cumulative(upper, linear_cumulative_upper_interval_);
    const auto lower_cumulative = legacy_linear_cumulative(lower, linear_cumulative_lower_interval_);
    const auto legacy = upper_cumulative - lower_cumulative;
    const auto subtraction_floor = 64.0 * std::numeric_limits<double>::epsilon()
                                   * std::max(std::abs(upper_cumulative), std::abs(lower_cumulative));
    constexpr double well_conditioned_margin = 1e10;
    if (std::isfinite(legacy) && legacy >= 0.0 && legacy > well_conditioned_margin * subtraction_floor)
      return legacy;

    const auto direct = checked_weight(integral_wide(lower, upper));
    const auto compatibility_tolerance = 1024.0 * std::numeric_limits<double>::epsilon()
                                         * std::max(direct, std::abs(legacy));
    if (std::isfinite(legacy) && legacy >= 0.0 && std::abs(direct - legacy) <= compatibility_tolerance)
      return legacy;
    return direct;
  }

  bool flag() const {
    require_initialized();
    return method_ == InterpolationMethod::linear ? linear_value_.flag() : interval_changed_;
  }

  void clear_flag() {
    require_initialized();
    if (method_ == InterpolationMethod::linear)
      linear_value_.clear_flag();
    else
      interval_changed_ = false;
  }
};

} // namespace NRG::Tools

#endif
