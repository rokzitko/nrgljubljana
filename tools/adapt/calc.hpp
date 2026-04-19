// Discretization ODE solver for NRG
// ** Integration code

#ifndef _adapt_calc_hpp_
#define _adapt_calc_hpp_

#include "../common/calc.hpp"

namespace NRG::Adapt {

inline auto integrate(Vec &vec) {
  return NRG::Tools::integrate_trapezoidal(vec);
}

inline auto integrate_ab(const Vec &vec, [[maybe_unused]] const double a, const double b) {
  return NRG::Tools::integrate_trapezoidal_ab(vec, a, b);
}

} // namespace

#endif
