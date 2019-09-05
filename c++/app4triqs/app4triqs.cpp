#include <cmath>
#include "./app4triqs.hpp"

namespace app4triqs {

  toto &toto::operator+=(toto const &b) {
    this->i += b.i;
    return *this;
  }

  toto toto::operator+(toto const &b) const {
    auto res = *this;
    res += b;
    return res;
  }

  bool toto::operator==(toto const &b) const { return (this->i == b.i); }

  int chain(int i, int j) {
    int n_digits_j = j > 0 ? (int)log10(j) + 1 : 1;
    return i * int(pow(10, n_digits_j)) + j;
  }

} // namespace app4triqs
