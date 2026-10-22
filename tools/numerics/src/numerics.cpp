#include <nrgljubljana/numerics.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

#include "../../common/gsl_piecewise_polynomial.hpp"
#include "../../hilb/hilb.hpp"
#include "../../integ/integ.hpp"
#include "../../kk/kk.hpp"

namespace NRG::Numerics {
namespace {

void validate_grid(const std::vector<double> &grid) {
  if (grid.empty()) throw std::invalid_argument("Grid must not be empty.");
  for (std::size_t i = 0; i < grid.size(); ++i) {
    if (!std::isfinite(grid[i])) throw std::invalid_argument("Grid coordinates must be finite.");
    if (i != 0 && !(grid[i - 1] < grid[i])) throw std::invalid_argument("Grid must be strictly increasing.");
  }
  if (!std::isfinite(grid.back() - grid.front())) throw std::invalid_argument("Grid width must be finite.");
}

void validate_table(const Table &table) {
  if (table.x.size() != table.y.size()) throw std::invalid_argument("Table x/y sizes must match.");
  validate_grid(table.x);
  for (const auto y : table.y)
    if (!std::isfinite(y)) throw std::invalid_argument("Table values must be finite.");
}

auto steffen(const Table &table) {
  return Tools::make_gsl_piecewise_polynomial(table.x, table.y, Tools::InterpolationMethod::steffen);
}

double finite_result(const double value) {
  if (!std::isfinite(value)) throw std::overflow_error("Numerical result is not finite.");
  return value;
}

} // namespace

double integrate(const Table &table) {
  return finite_result(steffen(table).integral());
}

double integrate_moment(const Table &table, const unsigned moment) {
  if (moment > static_cast<unsigned>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Moment must fit in a nonnegative int.");
  return finite_result(steffen(table).multiply_by_monomial(moment).integral());
}

double integrate_abs(const Table &table) {
  return finite_result(Tools::absolute_integral(steffen(table)));
}

double integrate_fermi(const Table &table, const double temperature, const double epsabs, const double epsrel) {
  if (!std::isfinite(temperature) || temperature < 0.0)
    throw std::invalid_argument("Temperature must be finite and nonnegative.");
  Tools::validate_tolerances(epsabs, epsrel);
  Integ::Options options;
  options.temperature = temperature;
  options.epsabs = epsabs;
  options.epsrel = epsrel;
  options.gsl_error_policy = Tools::GslErrorPolicy::fail;
  return finite_result(Integ::integrate_fermi(steffen(table), options));
}

Table resample(const Table &table, const std::vector<double> &grid) {
  const auto polynomial = steffen(table);
  validate_grid(grid);
  const auto lower = polynomial.lower_bound();
  const auto upper = polynomial.upper_bound();
  const auto scale = std::max({std::abs(lower), std::abs(upper), upper - lower});
  const auto tolerance = 32.0L * std::numeric_limits<double>::epsilon() * static_cast<long double>(scale);
  Table result{grid, {}};
  result.y.reserve(grid.size());
  for (const auto x : grid) {
    if (static_cast<long double>(x) < static_cast<long double>(lower) - tolerance
        || static_cast<long double>(x) > static_cast<long double>(upper) + tolerance)
      throw std::domain_error("Resampling would extrapolate outside the tabulated domain.");
    result.y.push_back(finite_result(polynomial.evaluate(std::clamp(x, lower, upper))));
  }
  return result;
}

Table kramers_kronig(const Table &table) {
  // The CLI library constructor sorts its input. Reject unsorted data first.
  validate_table(table);
  KK::XYFUNC input;
  input.reserve(table.x.size());
  for (std::size_t i = 0; i < table.x.size(); ++i) input.emplace_back(table.x[i], table.y[i]);
  KK::NumericalOptions options;
  options.interpolation = Tools::InterpolationMethod::steffen;
  options.algorithm = KK::Algorithm::analytic;
  options.jobs = 1;
  const KK::KK transform(std::move(input), options);
  Table result{table.x, {}};
  result.y.reserve(table.x.size());
  for (const auto x : table.x) result.y.push_back(finite_result(transform.calc(x)));
  return result;
}

Lattice lattice(const Table &dos, const Table &re_sigma, const Table &im_sigma, const double mu, const double clip_sigma,
                const unsigned moment) {
  validate_table(re_sigma);
  validate_table(im_sigma);
  if (re_sigma.x != im_sigma.x) throw std::invalid_argument("Self-energy frequency grids must match exactly.");
  if (!std::isfinite(mu)) throw std::invalid_argument("Chemical potential must be finite.");
  if (!std::isfinite(clip_sigma) || clip_sigma < Hilb::minimum_safe_imaginary_part())
    throw std::invalid_argument("Self-energy clipping must be finite and at least sqrt(min normal double).");
  if (moment > static_cast<unsigned>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Moment must fit in a nonnegative int.");
  const auto density = steffen(dos);
  for (const auto y : dos.y)
    if (y < 0.0) throw std::invalid_argument("DOS must be nonnegative.");
  // Prepare once, exactly as in hilb's analytic tabulated path. In particular,
  // n=1 must retain the far-field accuracy lost by subtracting z*G_0 - mass.
  const auto weighted_density = density.multiply_by_monomial(moment);
  const auto bandwidth = Hilb::bandwidth(dos.x);
  Lattice result{{re_sigma.x, {}}, {re_sigma.x, {}}};
  result.real.y.reserve(re_sigma.x.size());
  result.spectral.y.reserve(re_sigma.x.size());
  for (std::size_t i = 0; i < re_sigma.x.size(); ++i) {
    const std::complex<double> z{re_sigma.x[i] + mu - re_sigma.y[i], std::max(-im_sigma.y[i], clip_sigma)};
    Hilb::validate_hilbert_transform_argument(bandwidth, z, static_cast<int>(moment));
    const auto green = Tools::cauchy_transform(weighted_density, z);
    result.real.y.push_back(finite_result(green.real()));
    result.spectral.y.push_back(finite_result(-green.imag() / std::numbers::pi));
  }
  return result;
}

} // namespace NRG::Numerics
