// Discretization ODE solver for NRG
//
// ** Integration code
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: calc.h,v 1.1 2009/03/20 09:53:41 rok Exp $

#include <algorithm>

// Integrate a vector using the trapezoidal rule (this is consistent with
// linear interpolation between the tabulated function values). Returns the
// value of the integral over all data points.
double integrate(Vec &vec)
{
  Vec temp(vec);
  int len = vec.size();
  assert(len >= 2);
  double sum = 0.0;
  for (int i = 1; i < len; i++) {
    double dx = vec[i].first-vec[i-1].first;
    double yavg = (temp[i].second+temp[i-1].second)/2.0;
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
double integrate_ab(const Vec &vec, double a, double b)
{
  assert(a < b);

  Vec temp(vec); // We need to make a copy to sort the table.
  sort(temp.begin(), temp.end());

  int len = temp.size();
  double sum = 0.0;
  for (int i = 1; i < len; i++) {
    double x0 = temp[i-1].first;
    double x1 = temp[i].first;
    double yavg = (temp[i].second+temp[i-1].second)/2.0;
    double dx;
    if (x0 < b && x1 <= b) {
      dx = x1-x0;
    } else if (x0 < b && x1 > b) {
      dx = b - x0;
    } else {
      dx = 0.0;
    }
    sum += dx * yavg;
  }
  return sum;
}
