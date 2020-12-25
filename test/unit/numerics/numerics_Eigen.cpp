#include <gtest/gtest.h>
#include "compare.hpp"

#define INCL_EIGEN
#define USE_EIGEN
#undef INCL_UBLAS
#undef USE_UBLAS

#include "numerics.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics_Eigen, zero_matrix) {
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  Eigen::MatrixXd zero_m1 = NRG::zero_matrix<double>(dim1,dim2);
  Eigen::MatrixXd zero_m2 = NRG::zero_matrix<double>(dim3);
  ASSERT_EQ(zero_m1.rows(), dim1);
  ASSERT_EQ(zero_m1.cols(), dim2);
  ASSERT_EQ(zero_m2.rows(), dim3);
  ASSERT_EQ(zero_m2.cols(), dim3);
  for(size_t i = 0; i < dim1; i++)
    for(size_t j = 0; j < dim2; j++)
      EXPECT_EQ(zero_m1(i,j), 0);
  for(size_t i = 0; i < dim3; i++)
    for(size_t j = 0; j < dim3; j++)
      EXPECT_EQ(zero_m2(i,j), 0);
  // Test modification
  zero_m1(0,0) = 1;
  EXPECT_EQ(zero_m1(0,0), 1);
}

TEST(numerics_Eigen, trace_exp_real) {
  const size_t N = 3;
  Eigen::VectorXd v(N);
  for(size_t i = 0; i < N ; i++)
    v(i) = i + 1;
  Eigen::MatrixXd m(N, N);
  for(size_t i = 0; i < N*N; i++)
    m(i) = i + 1;
  double expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(-2.5 * v[i]) * m(i, i);
  compare(expected, trace_exp(v, m, 2.5));
}

TEST(numerics_Eigen, trace_exp_complex) {
  const size_t N = 3;
  Eigen::VectorXcd v(N);
  for(size_t j = 0; j < N ; j++)
    v(j) = j + 1.0 + 3i * (double)j;
  Eigen::MatrixXcd m(N, N);
  for(size_t j = 0; j < N*N; j++)
    m(j) = j + 3.0 + 2i * (double)j;
  std::complex<double> expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(- 2.5 * v(i)) * m(i, i);
  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics_Eigen, submatrix1) {
  const size_t N = 3; 
  Eigen::MatrixXd a(N,N);
  for(size_t i = 0; i < N*N; i++) a(i) = i;
  Eigen::MatrixXd dg(N - 1, N - 1);
  dg << 0, 3, 1, 4;
  const Eigen::MatrixXd dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics_Eigen, submatrix2) {
  Eigen::MatrixX<double> m(2,2);
  m(0,0) = m(0,1) = m(1,0) = m(1,1) = 0;
  auto sub = submatrix(m, {1,1}, {1,1});
  sub(0,0) = 1;
  EXPECT_DOUBLE_EQ(m(1,1), 1);
}

