#include <gtest/gtest.h>
#include <type_traits>
#include <traits.hpp>

using namespace NRG;

TEST(traits, double) {
  const bool b = std::is_same_v<typename traits<double>::t_matel, double>; // cannot be embedded in EXPECT_EQ
  EXPECT_EQ(b, true);
}

TEST(traits, complex) {
  const bool b = std::is_same_v<typename traits<std::complex<double>>::t_matel, std::complex<double>>;
  EXPECT_EQ(b, true);
}

template<matrix M> void foo(const M &m) {
  EXPECT_EQ(size1(m), 0);
  EXPECT_EQ(size2(m), 0);
}

TEST(traits, matrix) {
  EigenMatrix<double> er;
  foo(er);
  EigenMatrix<std::complex<double>> ec;
  foo(ec);
}

template<real_matrix RM> void foo_r(const RM &m) {
  EXPECT_EQ(size1(m), 0);
  EXPECT_EQ(size2(m), 0);
}

TEST(traits, real_matrix) {
  EigenMatrix<double> er;
  foo_r(er);
}

template<complex_matrix CM> void foo_c(const CM &m) {
  EXPECT_EQ(size1(m), 0);
  EXPECT_EQ(size2(m), 0);
}

TEST(traits, complex_matrix) {
  EigenMatrix<std::complex<double>> ec;
  foo_c(ec);
}

template<Eigen_matrix EM> void foo_e(const EM &m) {
  EXPECT_EQ(size1(m), 0);
  EXPECT_EQ(size2(m), 0);
}

TEST(traits, Eigen_matrix) {
  EigenMatrix<double> er;
  foo_e(er);
  EigenMatrix<std::complex<double>> ec;
  foo_e(ec);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
