#ifndef _tools_common_calc_hpp_
#define _tools_common_calc_hpp_

#include <algorithm>
#include <cassert>

namespace NRG::Tools {

template<typename Vec>
double integrate_trapezoidal(Vec &vec) {
  Vec temp(vec);
  const int len = vec.size();
  assert(len >= 2);
  double sum = 0.0;
  for (int i = 1; i < len; i++) {
    const auto [x0, y0] = temp[i - 1];
    const auto [x1, y1] = temp[i];
    const auto dx = x1 - x0;
    const auto yavg = (y1 + y0) / 2.0;
    sum += dx * yavg;
    vec[i].second = sum;
  }
  vec[0].second = 0.0;
  return sum;
}

template<typename Vec>
double integrate_trapezoidal_ab(const Vec &vec, [[maybe_unused]] const double a, const double b) {
  assert(a < b);
  Vec temp(vec);
  std::sort(temp.begin(), temp.end());
  const int len = temp.size();
  double sum = 0.0;
  for (int i = 1; i < len; i++) {
    const auto [x0, y0] = temp[i - 1];
    const auto [x1, y1] = temp[i];
    const auto yavg = (y1 + y0) / 2.0;
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

} // namespace NRG::Tools

#endif
