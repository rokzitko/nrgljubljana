#include <gtest/gtest.h>
#include "compare.hpp"

#ifdef USE_EIGEN

#include "numerics.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics, zero_matrix){
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  Eigen::MatrixXd zero_m1 = NRG::zero_matrix<double>(dim1,dim2);
  Eigen::MatrixXd zero_m2 = NRG::zero_matrix<double>(dim3);

  ASSERT_EQ(zero_m1.rows(), dim1);
  ASSERT_EQ(zero_m1.cols(), dim2);
  ASSERT_EQ(zero_m2.rows(), dim3);
  ASSERT_EQ(zero_m2.cols(), dim3);

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

TEST(numerics, Eigen_trace_exp_real){
  const int N = 3;
  Eigen::VectorXd v(N);
  for(int i = 0; i < N ; i++)
    v(i) = i + 1;
  
  Eigen::MatrixXd m(N, N);
  for(int i = 0; i < N*N; i++)
    m(i) = i + 1;
  
  double expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v[i]) * m(i, i);
  }

  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, Eigen_trace_exp_complex){
  const int N = 3;
  Eigen::VectorXcd v(N);
  for(int j = 0; j < N ; j++)
    v(j) = j + 1.0 + 3i * (double)j;
  
  Eigen::MatrixXcd m(N, N);
  for(int j = 0; j < N*N; j++)
    m(j) = j + 3.0 + 2i * (double)j;
  
  std::complex<double> expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v(i)) * m(i, i);
  }

  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, submatrix_eigen){
  const int N = 3; 
  Eigen::MatrixXd a(N,N);
  for(int i = 0; i < N*N; i++) a(i) = i;

  Eigen::MatrixXd dg(N - 1, N - 1);
  dg << 0, 3, 1, 4;

  Eigen::MatrixXd dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics, sum_of_exp){
  const int N = 3;
  double factor = 1.8;
  double sum = 0;
  Eigen::MatrixXd b(N, N);
  for(int i = 0; i < N*N; i++) b(i) = i;
  for(int i = 0; i < N*N; i++) sum += exp(-factor*b(i));

  double result_eigen = sum_of_exp(b, factor);
  compare(result_eigen, sum);

}

#endif