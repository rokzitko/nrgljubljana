#include <gtest/gtest.h>
#include <type_traits>
#include <complex>

#include <traits.hpp>
#include <numerics.hpp>

using namespace NRG;

TEST(generate, Eigen_d) {
  const auto m = generate_Eigen<double>(5, 10);
  EXPECT_EQ(nrvec(m), 5);
  EXPECT_EQ(dim(m), 10);
  EXPECT_EQ(size1(m), 5);
  EXPECT_EQ(size2(m), 10);
}

TEST(generate, ublas_d) {
  const auto m = generate_ublas<double>(5, 10);
  EXPECT_EQ(nrvec(m), 5);
  EXPECT_EQ(dim(m), 10);
  EXPECT_EQ(size1(m), 5);
  EXPECT_EQ(size2(m), 10);
}

TEST(generate, Eigen_c) {
  const auto m = generate_Eigen<std::complex<double>>(5, 10);
  EXPECT_EQ(nrvec(m), 5);
  EXPECT_EQ(dim(m), 10);
  EXPECT_EQ(size1(m), 5);
  EXPECT_EQ(size2(m), 10);
}

TEST(generate, ublas_c) {
  const auto m = generate_ublas<std::complex<double>>(5, 10);
  EXPECT_EQ(nrvec(m), 5);
  EXPECT_EQ(dim(m), 10);
  EXPECT_EQ(size1(m), 5);
  EXPECT_EQ(size2(m), 10);
}


int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
