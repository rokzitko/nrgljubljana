#include <gtest/gtest.h>
#include "compare.hpp"

#ifdef USE_UBLAS

#include "numerics.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics, zero_matrix){
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  auto zero_m1 = NRG::zero_matrix<double>(dim1,dim2);
  auto zero_m2 = NRG::zero_matrix<double>(dim3);

  ASSERT_EQ(zero_m1.size1(), dim1);
  ASSERT_EQ(zero_m1.size2(), dim2);
  ASSERT_EQ(zero_m2.size1(), dim3);
  ASSERT_EQ(zero_m2.size2(), dim3);

  for(size_t i = 0; i < dim1; i++){
    for(size_t j = 0; j < dim2; j++){
      EXPECT_EQ(zero_m1(i,j), 0);
    }
  }

  for(size_t i = 0; i < dim3; i++){
    for(size_t j = 0; j < dim3; j++){
      EXPECT_EQ(zero_m2(i,j), 0);
    }
  }

}

TEST(numerics, std_trace_exp_real){
  const int N = 3;
  std::vector<double> v(N);
  for(int i = 0; i < N ; i++)
    v[i] = i + 1;
  
  ublas::matrix<double> m(N, N);
  for(int i = 0; i < N*N; i++)
    m(i) = i + 1;
  
  double expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v[i]) * m(i, i);
  }

  EXPECT_DOUBLE_EQ(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, std_trace_exp_complex){
  const int N = 3;
  std::vector<std::complex<double>> v(N);
  for(int j = 0; j < N ; j++)
    v[j] = j + 1 + 3i * j;
  
  ublas::matrix<std::complex<double>> m(N, N);
  for(int j = 0; j < N*N; j++)
    m(j) = j + 3 + 2i * j;
  
  std::complex<double> expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v[i]) * m(i, i);
  }

  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, submatrix_ublas){
  const int N = 3; 
  ublas::matrix<double> a(N, N);
  for(int i = 0; i < N*N; i++) a(i) = i;
  
  ublas::matrix<double> dg(N - 1, N - 1);
  dg(0, 0) = 0;
  dg(0, 1) = 1;
  dg(1, 0) = 3;
  dg(1, 1) = 4;
  ublas::matrix<double> dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics, sum_of_exp){
  const int N = 3;
  double factor = 1.8;
  ublas::matrix<double> a(N, N);
  for(int i = 0; i < N*N; i++) a(i) = i;
  double result_ublas = sum_of_exp(a.data(), factor);
  double sum = 0;
  for(int i = 0; i < N*N; i++) sum += exp(-factor*a(i));
  compare(result_ublas, sum);
}

#endif