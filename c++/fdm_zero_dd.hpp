#ifndef NRG_FDM_ZERO_DD_HPP
#define NRG_FDM_ZERO_DD_HPP

#include <cmath>
#include <complex>
#include <type_traits>

namespace NRG {

template<typename Range, typename Values, typename Matrix, typename Factor>
[[nodiscard]] bool fdm_zero_dd_block(const Range &di, const Range &dj, const Values &boltz_i, const Values &boltz_j,
                                     const Values &energy_i, const Values &energy_j, const Matrix &op1, const Matrix &op2,
                                     double temperature, double wnf, Factor factor) {
  // op1 * op2 * wnf is evaluated before the zero Boltzmann factor. This
  // conservative component bound keeps even the complex intermediate products finite.
  constexpr double bound = 0x1p255;
  if (!(temperature > 0 && std::isfinite(temperature) && std::abs(wnf) <= bound &&
        std::isfinite(std::real(factor)) && std::isfinite(std::imag(factor)))) return false;

  const auto all_zero = [](const Range &d, const Values &boltz, const Values &energy) {
    if (d.empty()) return false;
    for (const auto k : d)
      if (boltz[k] != 0.0 || !std::isfinite(energy[k])) return false;
    return true;
  };
  if (!all_zero(di, boltz_i, energy_i) || !all_zero(dj, boltz_j, energy_j)) return false;

  const auto safe = [&](const Matrix &op) {
    const auto block = op.block(*dj.begin(), *di.begin(), dj.size(), di.size());
    if (!(block.real().array().abs() <= bound).all()) return false;
    if constexpr (!std::is_arithmetic_v<typename Matrix::Scalar>) return (block.imag().array().abs() <= bound).all();
    return true;
  };
  return safe(op1) && (&op1 == &op2 || safe(op2));
}

} // namespace NRG
#endif
