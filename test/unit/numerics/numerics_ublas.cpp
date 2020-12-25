#include <gtest/gtest.h>

#define USE_UBLAS
#include <traits.hpp>
#include <numerics.hpp>

#include "compare.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics_ublas, generate) {
  const auto m1 = generate_ublas<double>(2,3);
  EXPECT_EQ(size1(m1), 2);
  EXPECT_EQ(size2(m1), 3);
}

TEST(numerics_ublas, zero_matrix) {
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  auto zero_m1 = NRG::zero_matrix<double>(dim1,dim2); // non const
  const auto zero_m2 = NRG::zero_matrix<double>(dim3);
  ASSERT_EQ(zero_m1.size1(), dim1);
  ASSERT_EQ(zero_m1.size2(), dim2);
  ASSERT_EQ(zero_m2.size1(), dim3);
  ASSERT_EQ(zero_m2.size2(), dim3);
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

TEST(numerics_ublas, std_trace_exp_real) {
  const size_t N = 3;
  std::vector<double> v(N);
  for(size_t i = 0; i < N; i++)
    v[i] = i + 1;
  ublas::matrix<double> m(N, N);
  for(size_t i = 0; i < N*N; i++)
    m(i) = i + 1;
  double expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(-2.5 * v[i]) * m(i, i);
  EXPECT_DOUBLE_EQ(expected, trace_exp(v, m, 2.5));
}

TEST(numerics_ublas, std_trace_exp_complex){
  const size_t N = 3;
  std::vector<std::complex<double>> v(N);
  for(size_t j = 0; j < N ; j++)
    v[j] = j + 1 + 3i * j;
  ublas::matrix<std::complex<double>> m(N, N);
  for(size_t j = 0; j < N*N; j++)
    m(j) = j + 3 + 2i * j;
  std::complex<double> expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(-2.5 * v[i]) * m(i, i);
  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics_ublas, sum_of_exp) {
  const size_t N = 3;
  const double factor = 1.8;
  ublas::matrix<double> a(N, N);
  for(size_t i = 0; i < N*N; i++) a(i) = i;
  const double result_ublas = sum_of_exp(a.data(), factor);
  double sum = 0;
  for(size_t i = 0; i < N*N; i++) sum += exp(-factor*a(i));
  compare(result_ublas, sum);
}

TEST(numerics_ublas, submatrix1) {
  const size_t N = 3; 
  ublas::matrix<double> a(N, N);
  for(size_t i = 0; i < N*N; i++) a(i) = i;
  ublas::matrix<double> dg(N - 1, N - 1);
  dg(0, 0) = 0;
  dg(0, 1) = 1;
  dg(1, 0) = 3;
  dg(1, 1) = 4;
  ublas::matrix<double> dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics_ublas, submatrix2) {
  ublas::matrix<double> m(2,2);
  m(0,0) = m(0,1) = m(1,0) = m(1,1) = 0;
  auto sub = submatrix(m, {1,1}, {1,1});
  sub(0,0) = 1;
  EXPECT_DOUBLE_EQ(m(1,1), 1);
}

TEST(numerics_Eigen, data_r) {
  ublas::matrix<double> m(2, 3);
  m(0,0) = 1.;
  m(0,1) = 2.;
  m(0,2) = 3.;
  m(1,0) = 4.;
  m(1,1) = 5.;
  m(1,2) = 6.;
  double * d = data(m);
  EXPECT_EQ(*(d+0), 1.);
  EXPECT_EQ(*(d+1), 2.);
  EXPECT_EQ(*(d+2), 3.);
  EXPECT_EQ(*(d+3), 4.);
  EXPECT_EQ(*(d+4), 5.);
  EXPECT_EQ(*(d+5), 6.);
}

TEST(numerics_Eigen, data_c) {
  ublas::matrix<std::complex<double>> m(2, 3);
  m(0,0) = 1.*(1.+1.i);
  m(0,1) = 2.*(1.+1.i);
  m(0,2) = 3.*(1.+1.i);
  m(1,0) = 4.*(1.+1.i);
  m(1,1) = 5.*(1.+1.i);
  m(1,2) = 6.*(1.+1.i);
  std::complex<double> * d = data(m);
  EXPECT_EQ(*(d+0), 1.*(1.+1.i));
  EXPECT_EQ(*(d+1), 2.*(1.+1.i));
  EXPECT_EQ(*(d+2), 3.*(1.+1.i));
  EXPECT_EQ(*(d+3), 4.*(1.+1.i));
  EXPECT_EQ(*(d+4), 5.*(1.+1.i));
  EXPECT_EQ(*(d+5), 6.*(1.+1.i));
}
