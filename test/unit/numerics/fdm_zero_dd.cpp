#include <complex>
#include <limits>
#include <vector>
#include <boost/range/irange.hpp>
#include <gtest/gtest.h>
#include "traits.hpp"
#include "fdm_zero_dd.hpp"

using namespace NRG;

TEST(FdmZeroDD, RequiresExactZerosInBothDiscardedRanges) {
  const auto d = boost::irange(1, 3), all = boost::irange(0, 3), empty = boost::irange(3, 3);
  const std::vector<double> boltz{1, 0, -0.0}, energy{0, 1, 2};
  const EigenMatrix<double> op = EigenMatrix<double>::Ones(3, 3);
  const auto skip = [&](auto di, auto dj, const auto &bi, const auto &bj) {
    return fdm_zero_dd_block(di, dj, bi, bj, energy, energy, op, op, 1e-10, 1.0, 1.0);
  };
  EXPECT_TRUE(skip(d, d, boltz, boltz));
  EXPECT_FALSE(skip(empty, d, boltz, boltz));
  EXPECT_FALSE(skip(d, empty, boltz, boltz));
  EXPECT_FALSE(skip(all, all, boltz, boltz)); // FDM's last-shell Drange includes the ground state
  for (double live : {1.0, std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::quiet_NaN()}) {
    const std::vector<double> partial{1, 0, live};
    EXPECT_FALSE(skip(d, d, partial, boltz));
    EXPECT_FALSE(skip(d, d, boltz, partial));
  }
}

TEST(FdmZeroDD, UnsafePreproductsRetainOriginalLoop) {
  using Z = std::complex<double>;
  const auto di = boost::irange(1, 3), dj = boost::irange(2, 4);
  const std::vector<double> bi{1, 0, 0}, bj{1, 1, 0, 0}, ej{0, 0, 1, 2};
  std::vector<double> ei{0, 1, 2};
  const auto inf = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  EigenMatrix<Z> a = EigenMatrix<Z>::Constant(4, 3, Z(2, -3)), b = a;
  const auto skip = [&](double T, double wnf, Z factor) {
    return fdm_zero_dd_block(di, dj, bi, bj, ei, ej, a, b, T, wnf, factor);
  };
  a(0, 0) = Z(nan, 0); // outside the DD rectangle: not part of this decision
  EXPECT_TRUE(skip(1e-10, 1, Z(1, -1)));
  EXPECT_FALSE(skip(0, 1, 1));
  EXPECT_FALSE(skip(inf, 1, 1));
  EXPECT_FALSE(skip(1e-10, nan, 1));
  EXPECT_FALSE(skip(1e-10, 0x1p256, 1));
  EXPECT_FALSE(skip(1e-10, 1, Z(0, inf)));
  ei.back() = inf;
  EXPECT_FALSE(skip(1e-10, 1, 1));
  ei.back() = 2;
  for (auto *op : {&a, &b}) {
    for (Z unsafe : {Z(nan, 0), Z(0, inf), Z(0x1p256, 0), Z(0, -0x1p256)}) {
      (*op)(3, 2) = unsafe;
      EXPECT_FALSE(skip(1e-10, 1, 1));
    }
    (*op)(3, 2) = 1;
  }
}
