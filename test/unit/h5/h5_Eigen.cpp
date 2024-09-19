#include <gtest/gtest.h>

#include <type_traits>
#include <complex>

#include <Eigen/Dense>

#define H5_USE_EIGEN
#include <highfive/H5File.hpp>
#include <h5.hpp>

#include "compare.hpp"

TEST(h5dump, eigen_dump_matrix) { // NOLINT
  H5Easy::File file("eigen_matrix.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  EigenMatrix<double> w(nx,ny);
  double cnt = 1.0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  NRG::h5_dump_matrix(file, "/path", w);

  auto r = H5Easy::load<EigenMatrix<double>>(file, "/path");
  compare(w,r);
}

TEST(h5dump, eigen_matrix_real_part_rect) { // NOLINT
  H5Easy::File file("eigen_matrix_real_part_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  EigenMatrix<std::complex<double>> w(nx,ny);
  auto cnt = std::complex<double>(1.0,1.0);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  NRG::h5_dump_matrix(file, "/path", w);

  auto r = H5Easy::load<EigenMatrix<double>>(file, "/path");
  for (int x = 0; x < nx; x++)
    for (int y = 0; y < ny; y++)
      EXPECT_EQ(w(x,y).real(), r(x,y));
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
