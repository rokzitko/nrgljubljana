#include <string>
#include <sstream>
#include <cmath>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <params.hpp>
#include <diag.hpp>

using namespace NRG;

TEST(Diag, diagonalise) { // NOLINT
  Params P;
  auto DP = DiagParams(P, 0.5);
}

TEST(Diag, to_matel) {
  lapack_complex_double z1 {2.0, 1.0};
  const auto z2 = to_matel(z1);
  EXPECT_DOUBLE_EQ(z2.real(), 2.0);
  EXPECT_DOUBLE_EQ(z2.imag(), 1.0);
}

TEST(Diag, copy_val) {
  const std::vector<double> s { 1.0, 2.0, 3.0 };
  std::vector<double> d;
  copy_val(s, d, s.size());
  EXPECT_EQ(d.size(), 3);
  EXPECT_DOUBLE_EQ(d[0], 1.0);
  EXPECT_DOUBLE_EQ(d[1], 2.0);
  EXPECT_DOUBLE_EQ(d[2], 3.0);
}

TEST(Diag, copy_val_M) {
  const std::vector<double> s { 1.0, 2.0, 3.0 };
  std::vector<double> d;
  copy_val(s, d, 2);
  EXPECT_EQ(d.size(), 2);
  EXPECT_DOUBLE_EQ(d[0], 1.0);
  EXPECT_DOUBLE_EQ(d[1], 2.0);
}

TEST(Diag, copy_vec) {
  const std::vector<double> m { 1.0, 2.0, 3.0, 4.0 };
  Matrix_traits<double> vec;
  copy_vec(data(m), vec, 2, 2);
  EXPECT_EQ(size1(vec), 2);
  EXPECT_EQ(size2(vec), 2);
  EXPECT_DOUBLE_EQ(vec(0,0), 1.0);
  EXPECT_DOUBLE_EQ(vec(0,1), 2.0);
  EXPECT_DOUBLE_EQ(vec(1,0), 3.0);
  EXPECT_DOUBLE_EQ(vec(1,1), 4.0);
}

TEST(Diag, copy_results) {
  const std::vector<double> val { 1.0, 2.0 };
  const std::vector<double> vec { 1.0, 0.0, 0.0, 1.0 };
  const auto res = copy_results<double>(val, data(vec), 'V', 2, 2);
  res.check_diag();
  EXPECT_EQ(res.getnrcomputed(), 2);
  EXPECT_EQ(res.getdim(), 2);
}

TEST(Diag, check_is_matrix_upper) {
  Matrix_traits<double> m(2,2);
  m(0,0) = m(0,1) = m(1,1) = 1.0;
  m(1,0) = 0.0;
  EXPECT_TRUE(is_square(m));
  EXPECT_TRUE(is_matrix_upper(m));
}

constexpr double sq2 = 1.0/sqrt(2.0);

TEST(Diag, dsyev) {
  {
    Matrix_traits<double> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise_dsyev(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(res.vec(0,1), 0.0);
    EXPECT_DOUBLE_EQ(res.vec(1,0), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<double> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise_dsyev(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(0,0)*res.vec(0,1), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(1,0)*res.vec(1,1), +0.5);
  }
}

TEST(Diag, dsyevr) {
  {
    Matrix_traits<double> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise_dsyevr(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(res.vec(0,1), 0.0);
    EXPECT_DOUBLE_EQ(res.vec(1,0), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<double> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise_dsyevr(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(0,0)*res.vec(0,1), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(1,0)*res.vec(1,1), +0.5);
  }
}

TEST(Diag, dsyevd) {
  {
    Matrix_traits<double> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise_dsyevd(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(res.vec(0,1), 0.0);
    EXPECT_DOUBLE_EQ(res.vec(1,0), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<double> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise_dsyevd(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(0,0)*res.vec(0,1), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ(res.vec(1,0)*res.vec(1,1), +0.5);
  }
}

TEST(Diag, zheev) {
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise_zheev(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise_zheev(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).real(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).real(), +0.5);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Y
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = std::complex<double>(0.0,1.0);
    const auto res = diagonalise_zheev(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).imag(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).imag(), +0.5);
  }
}

TEST(Diag, zheevr) {
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise_zheevr(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise_zheevr(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).real(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).real(), +0.5);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Y
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = std::complex<double>(0.0,1.0);
    const auto res = diagonalise_zheevr(m);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).imag(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).imag(), +0.5);
  }
}

TEST(Diag, diagonalise_complex) {
  Params P;
  auto DP = DiagParams(P, 0.5);
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Z
    m(0,0) = -1.0; m(1,1) = +1.0;
    m(1,0) = m(0,1) = 0.0;
    const auto res = diagonalise(m, DP, -1);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), 1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), 0.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), 1.0);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_X
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = 1.0;
    const auto res = diagonalise(m, DP, -1);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).real(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).real(), +0.5);
  }
  {
    Matrix_traits<std::complex<double>> m(2,2); // Sigma_Y
    m(0,0) = m(1,1) = m(1,0) = 0.0;
    m(0,1) = std::complex<double>(0.0,1.0);
    const auto res = diagonalise(m, DP, -1);
    EXPECT_EQ(res.getnrcomputed(), 2);
    EXPECT_EQ(res.getdim(), 2);
    EXPECT_DOUBLE_EQ(res.val[0], -1.0);
    EXPECT_DOUBLE_EQ(res.val[1], +1.0);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(0,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(0,0)*res.vec(0,1)).imag(), -0.5);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,0)), sq2);
    EXPECT_DOUBLE_EQ(abs(res.vec(1,1)), sq2);
    EXPECT_DOUBLE_EQ((res.vec(1,0)*res.vec(1,1)).imag(), +0.5);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}

