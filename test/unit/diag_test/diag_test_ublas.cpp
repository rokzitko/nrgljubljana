#include <string>
#include <sstream>
#include <cmath>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <params.hpp>
#include <diag.hpp>

using namespace NRG;

const double sq2 = 1.0/sqrt(2.0);

TEST(Diag, dsyev) {
  {
    ublasMatrix<double> m(2,2); // Sigma_Z
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
    ublasMatrix<double> m(2,2); // Sigma_X
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
    ublasMatrix<double> m(2,2); // Sigma_Z
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
    ublasMatrix<double> m(2,2); // Sigma_X
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
    ublasMatrix<double> m(2,2); // Sigma_Z
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
    ublasMatrix<double> m(2,2); // Sigma_X
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Z
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_X
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Y
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Z
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_X
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Y
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Z
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_X
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
    ublasMatrix<std::complex<double>> m(2,2); // Sigma_Y
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

