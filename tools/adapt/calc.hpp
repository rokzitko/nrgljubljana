// Discretization ODE solver for NRG
// ** Integration code

#ifndef _adapt_calc_hpp_
#define _adapt_calc_hpp_

#include <vector>
#include <algorithm>
#include <cassert>

namespace NRG::Adapt {

// Integrate a vector using the trapezoidal rule (this is consistent with
// linear interpolation between the tabulated function values). Returns the
// value of the integral over all data points.
auto integrate(Vec &vec) {
  Vec temp(vec);
  const int len = vec.size();
  assert(len >= 2);
  auto sum = 0.0;
  for (int i = 1; i < len; i++) {
    const auto [x0, y0] = temp[i-1];
    const auto [x1, y1] = temp[i];
    const auto dx   = x1-x0;
    const auto yavg = (y1+y0) / 2.0;
    sum += dx * yavg;
    vec[i].second = sum; // total so far
  }
  vec[0].second = 0.0;
  return sum;
}

// Integrate a tabulated function using the trapezoidal rule on the
// interval [a:b]. It is assued that the upper boundary 'b' is within the
// interval of tabulation, while the lower boundary 'a' is below it
// (typically a=0); to obtain the contribution from 'a' to the beginning of
// the tabulation interval, an extrapolation is performed.
auto integrate_ab(const Vec &vec, const double a, const double b) {
  assert(a < b);
  Vec temp(vec); // We need to make a copy to sort the table.
  std::sort(temp.begin(), temp.end());
  const int len = temp.size();
  auto sum = 0.0;
  for (int i = 1; i < len; i++) {
    const auto [x0, y0] = temp[i-1];
    const auto [x1, y1] = temp[i];
    const auto yavg = (y1+y0) / 2.0;
    double dx;
    if (x0 < b && x1 <= b) {
      dx = x1 - x0;
    } else if (x0 < b && x1 > b) {
      dx = b - x0;
    } else {
      dx = 0.0;
    }
    sum += dx * yavg;
  }
  return sum;
}

} // namespace
   
#endif
