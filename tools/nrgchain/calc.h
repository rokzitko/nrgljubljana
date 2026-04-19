// Discretization ODE solver for NRG
//
// ** Integration code

#include "../common/calc.hpp"

inline double integrate(Vec &vec) {
  return NRG::Tools::integrate_trapezoidal(vec);
}

inline double integrate_ab(const Vec &vec, [[maybe_unused]] double a, double b) {
  return NRG::Tools::integrate_trapezoidal_ab(vec, a, b);
}
