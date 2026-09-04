#ifndef _floquet_hpp_
#define _floquet_hpp_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Eigenvalues>
#include <fmt/format.h>

#include "eigen.hpp"
#include "operators.hpp"
#include "portabil.hpp"
#include "traits.hpp"

namespace NRG {

struct FloquetCriteria {
  std::map<Invar, std::vector<double>> values;
  double minimum_energy = std::numeric_limits<double>::infinity();
};

inline double floquet_numerical_tolerance(const double scale, const size_t dimension) {
  return 256.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, scale) * std::max<size_t>(1, dimension);
}

template<scalar S>
auto floquet_conjugate(const S value) {
  if constexpr (is_complex<S>::value)
    return std::conj(value);
  else
    return value;
}

template<scalar S>
double floquet_imaginary_part(const S value) {
  if constexpr (is_complex<S>::value)
    return value.imag();
  else
    return 0.0;
}

template<scalar S>
double validate_floquet_matrix(const EigenMatrix<S> &matrix, const std::string &description) {
  if (matrix.rows() == 0 || matrix.cols() == 0)
    throw std::invalid_argument(fmt::format("{} is empty", description));

  double scale = 0.0;
  for (::Eigen::Index row = 0; row < matrix.rows(); row++)
    for (::Eigen::Index column = 0; column < matrix.cols(); column++) {
      const auto value = matrix(row, column);
      if (!my_isfinite(value))
        throw std::invalid_argument(fmt::format("{} contains a non-finite element", description));
      scale = std::max(scale, std::abs(value));
    }

  if (matrix.rows() != matrix.cols())
    throw std::invalid_argument(fmt::format("{} must be square", description));
  const auto tolerance = floquet_numerical_tolerance(scale, static_cast<size_t>(matrix.rows()));
  for (::Eigen::Index row = 0; row < matrix.rows(); row++)
    for (::Eigen::Index column = 0; column <= row; column++)
      if (std::abs(matrix(row, column) - floquet_conjugate(matrix(column, row))) > tolerance)
        throw std::invalid_argument(fmt::format("{} is not Hermitian", description));
  return scale;
}

template<scalar S>
std::pair<double, double> validate_floquet_mode_operator(const DiagInfo<S> &diag,
                                                         const Operators<S> &operators) {
  const size_t container_count = operators.ops.count("m") + operators.opsp.count("m") + operators.opsg.count("m") +
                                 operators.opd.count("m") + operators.opt.count("m") + operators.opq.count("m") +
                                 operators.opot.count("m");
  if (container_count != 1 || operators.ops.count("m") != 1)
    throw std::invalid_argument("Floquet operator m must occur only as an even-parity singlet");
  if (diag.empty()) throw std::invalid_argument("Floquet seed eigenspectrum has no invariant sectors");

  const auto &mode = operators.ops.at("m");
  if (mode.size() != diag.size())
    throw std::invalid_argument("Floquet operator m must have one diagonal block for every seed invariant sector");

  double lower = std::numeric_limits<double>::infinity();
  double upper = -std::numeric_limits<double>::infinity();
  size_t states = 0;
  for (const auto &[I, eig] : diag) {
    const auto count = eig.getnrstored();
    if (count == 0)
      throw std::invalid_argument(fmt::format("Floquet seed sector {} has no states", I.str()));
    states += count;
    const auto block = mode.find(Twoinvar(I, I));
    if (block == mode.end())
      throw std::invalid_argument(fmt::format("Floquet operator m is missing seed-sector block ({},{})", I.str(), I.str()));
    const auto &matrix = block->second;
    if (size1(matrix) != count || size2(matrix) != count)
      throw std::invalid_argument(fmt::format("Floquet operator m block ({},{}) has the wrong dimensions", I.str(), I.str()));
    validate_floquet_matrix(matrix, fmt::format("Floquet operator m block ({},{})", I.str(), I.str()));

    ::Eigen::SelfAdjointEigenSolver<EigenMatrix<S>> solver(matrix, ::Eigen::EigenvaluesOnly);
    if (solver.info() != ::Eigen::Success)
      throw std::invalid_argument(fmt::format("Could not determine the spectrum of Floquet operator m in sector {}", I.str()));
    lower = std::min(lower, solver.eigenvalues().minCoeff());
    upper = std::max(upper, solver.eigenvalues().maxCoeff());
  }

  for (const auto &[sectors, matrix] : mode) {
    const auto &[I1, I2] = sectors;
    if (I1 != I2 || diag.count(I1) == 0)
      throw std::invalid_argument(fmt::format("Floquet operator m contains unexpected block ({},{})", I1.str(), I2.str()));
  }
  if (states == 0 || !std::isfinite(lower) || !std::isfinite(upper) || lower > upper)
    throw std::invalid_argument("Floquet operator m has no finite spectral range");
  return {lower, upper};
}

template<scalar S>
FloquetCriteria calculate_floquet_criteria(const DiagInfo<S> &diag,
                                           const MatrixElements<S> &mode,
                                           const double Omega,
                                           const std::pair<double, double> mode_bounds) {
  if (!std::isfinite(Omega) || Omega <= 0.0)
    throw std::invalid_argument("Scaled Floquet Omega must be finite and positive");
  if (!std::isfinite(mode_bounds.first) || !std::isfinite(mode_bounds.second) || mode_bounds.first > mode_bounds.second)
    throw std::invalid_argument("Floquet operator m has invalid spectral bounds");
  if (diag.empty()) throw std::invalid_argument("Floquet eigenspectrum has no invariant sectors");
  if (mode.size() != diag.size())
    throw std::invalid_argument("Recalculated Floquet operator m must have one block for every invariant sector");

  FloquetCriteria result;
  size_t states = 0;
  for (const auto &[I, eig] : diag) {
    const auto count = eig.getnrcomputed();
    if (count == 0)
      throw std::invalid_argument(fmt::format("Floquet sector {} has no computed states", I.str()));
    states += count;
    const auto block = mode.find(Twoinvar(I, I));
    if (block == mode.end())
      throw std::invalid_argument(fmt::format("Recalculated Floquet operator m is missing block ({},{})", I.str(), I.str()));
    const auto &matrix = block->second;
    if (size1(matrix) != count || size2(matrix) != count)
      throw std::invalid_argument(fmt::format("Recalculated Floquet operator m block ({},{}) has the wrong dimensions", I.str(), I.str()));
    const auto matrix_scale = validate_floquet_matrix(
      matrix, fmt::format("Recalculated Floquet operator m block ({},{})", I.str(), I.str()));
    const auto bound_scale = std::max(std::abs(mode_bounds.first), std::abs(mode_bounds.second));
    const auto tolerance = floquet_numerical_tolerance(std::max(matrix_scale, bound_scale), count);

    auto &criteria = result.values[I];
    criteria.reserve(count);
    for (size_t index = 0; index < count; index++) {
      const auto energy = eig.values.raw(index);
      if (!std::isfinite(energy))
        throw std::invalid_argument(fmt::format("Floquet sector {} contains a non-finite energy", I.str()));
      const auto expectation_value = matrix(static_cast<::Eigen::Index>(index), static_cast<::Eigen::Index>(index));
      const auto expectation = std::real(expectation_value);
      if (std::abs(floquet_imaginary_part(expectation_value)) > tolerance)
        throw std::invalid_argument(fmt::format("Floquet mode expectation in sector {} is not real", I.str()));
      if (expectation < mode_bounds.first - tolerance || expectation > mode_bounds.second + tolerance)
        throw std::invalid_argument(fmt::format("Floquet mode expectation {} in sector {} is outside [{},{}]",
                                                expectation, I.str(), mode_bounds.first, mode_bounds.second));

      const auto mode_energy = expectation * Omega;
      const auto mean_energy = energy - mode_energy;
      const auto criterion = mean_energy + std::abs(mode_energy);
      if (!std::isfinite(mode_energy) || !std::isfinite(mean_energy) || !std::isfinite(criterion))
        throw std::invalid_argument(fmt::format("Floquet truncation criterion in sector {} is non-finite", I.str()));
      result.minimum_energy = std::min(result.minimum_energy, energy);
      criteria.push_back(criterion);
    }
  }

  for (const auto &[sectors, matrix] : mode) {
    const auto &[I1, I2] = sectors;
    if (I1 != I2 || diag.count(I1) == 0)
      throw std::invalid_argument(fmt::format("Recalculated Floquet operator m contains unexpected block ({},{})",
                                              I1.str(), I2.str()));
  }
  if (states == 0 || !std::isfinite(result.minimum_energy))
    throw std::invalid_argument("Floquet eigenspectrum contains no finite states");
  return result;
}

} // namespace NRG

#endif
