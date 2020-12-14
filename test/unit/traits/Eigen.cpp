#include <gtest/gtest.h>
#include <type_traits>
#include <complex>

#include <traits.hpp>
#include <Eigen/Dense>

using namespace NRG;

TEST(Eigen, default_constructor) {
  Eigen::MatrixX<double> md;
  Eigen::MatrixX<std::complex<double>> mc;
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
