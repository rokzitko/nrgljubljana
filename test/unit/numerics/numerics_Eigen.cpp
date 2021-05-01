#include <gtest/gtest.h>

#define USE_EIGEN
#include <traits.hpp>
#include <numerics.hpp>

#include "compare.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(Eigen, default_constructor) {
  EigenMatrix<double> md;
  EigenMatrix<std::complex<double>> mc;
}

TEST(eigen, basic_vector) {
  Eigen::Vector3i a(5,7,9);
  std::vector b = {5,7,9};
  compare(a,b);
  Eigen::Vector3i c;
  c << 5,7,9;
  compare(c,b);
}

TEST(eigen, matrix_ops) {
  EigenMatrix<double> m(2,2);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  EXPECT_DOUBLE_EQ(m(0,0), 1);
  EXPECT_DOUBLE_EQ(m(0,1), 2);
  EXPECT_DOUBLE_EQ(m(1,0), 3);
  EXPECT_DOUBLE_EQ(m(1,1), 4);
}

TEST(numerics_Eigen, zero_matrix) {
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  EigenMatrix<double> zero_m1 = NRG::zero_matrix<double>(dim1,dim2);
  EigenMatrix<double> zero_m2 = NRG::zero_matrix<double>(dim3);
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
  using T = double;
  const size_t N = 3;
  Eigen::VectorX<T> v(N);
  for(size_t i = 0; i < N ; i++)
    v(i) = i + 1;
  EigenMatrix<T> m(N, N);
  for(size_t i = 0; i < N*N; i++)
    m(i) = i + 1;
  T expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(-2.5 * v[i]) * m(i, i);
  compare(expected, trace_exp(v, m, 2.5));
}

TEST(numerics_Eigen, trace_exp_complex) {
  using T = std::complex<double>;
  const size_t N = 3;
  Eigen::VectorX<T> v(N);
  for(size_t j = 0; j < N ; j++)
    v(j) = j + 1.0 + 3i * (double)j;
  EigenMatrix<T> m(N, N);
  for(size_t j = 0; j < N*N; j++)
    m(j) = j + 3.0 + 2i * (double)j;
  T expected = 0;
  for(size_t i = 0; i < N; i++)
    expected += exp(- 2.5 * v(i)) * m(i, i);
  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics_Eigen, submatrix1) {
  const size_t N = 3; 
  EigenMatrix<double> a(N,N);
  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      a(i, j) = i+10*j;
  EigenMatrix<double> dg(N-1, N-1);
  dg(0,0) = 0;
  dg(0,1) = 10;
  dg(1,0) = 1;
  dg(1,1) = 11;
  const EigenMatrix<double> dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics_Eigen, submatrix2) {
  EigenMatrix<double> m(2,2);
  m(0,0) = m(0,1) = m(1,0) = m(1,1) = 0;
  auto sub = submatrix(m, {1,2}, {1,2});
  EXPECT_EQ(size1(sub), 1);
  EXPECT_EQ(size2(sub), 1);
  sub(0,0) = 1;
  EXPECT_DOUBLE_EQ(m(1,1), 1);
}

TEST(numerics_Eigen, submatrix3) {
  using T = std::complex<double>;
  const size_t dim1 = 15;
  const size_t dim2 = 25;
  EigenMatrix<T> m(dim1, dim2);
  EXPECT_EQ(size1(m), dim1);
  EXPECT_EQ(size2(m), dim2);

  for (size_t i = 0; i < dim1; i++)
    for (size_t j = 0; j < dim2; j++)
      m(i, j) = T(i, j);

  for (size_t offset1 = 0; offset1 < dim1; offset1++) {
    for (size_t offset2 = 0; offset2 < dim2; offset2++) {
      for (size_t sz1 = 0; sz1 < std::min(dim1, dim1-offset1); sz1++) {
        for (size_t sz2 = 0; sz2 < std::min(dim2, dim2-offset2); sz2++) {
          const auto sm1 = m.block(offset1, offset2, sz1, sz2);
          //EXPECT_EQ(size1(sm1), sz1);
          //EXPECT_EQ(size2(sm1), sz2);
          for (size_t i = 0; i < sz1; i++)
            for (size_t j = 0; j < sz2; j++)
              EXPECT_EQ(sm1(i,j), T(offset1+i, offset2+j));

          const auto sm2 = submatrix_const(m, {offset1, offset1+sz1}, {offset2, offset2+sz2});
          EXPECT_EQ(size1(sm2), sz1);
          EXPECT_EQ(size2(sm2), sz2);
          for (size_t i = 0; i < sz1; i++)
            for (size_t j = 0; j < sz2; j++)
              EXPECT_EQ(sm2(i,j), T(offset1+i, offset2+j));
        }
      }
    }
  }
}

TEST(numerics_Eigen, data_r) {
  EigenMatrix<double> m(2, 3);
  m(0,0) = 1.;
  m(0,1) = 2.;
  m(0,2) = 3.;
  m(1,0) = 4.;
  m(1,1) = 5.;
  m(1,2) = 6.;
  double * d = data(m);
  if (is_row_ordered(m)) {
    EXPECT_EQ(*(d+0), 1.);
    EXPECT_EQ(*(d+1), 2.);
    EXPECT_EQ(*(d+2), 3.);
    EXPECT_EQ(*(d+3), 4.);
    EXPECT_EQ(*(d+4), 5.);
    EXPECT_EQ(*(d+5), 6.);
  } else {
    EXPECT_EQ(*(d+0), 1.);
    EXPECT_EQ(*(d+1), 4.);
    EXPECT_EQ(*(d+2), 2.);
    EXPECT_EQ(*(d+3), 5.);
    EXPECT_EQ(*(d+4), 3.);
    EXPECT_EQ(*(d+5), 6.);
  }
}

TEST(numerics_Eigen, data_c) {
  EigenMatrix<std::complex<double>> m(2, 3);
  m(0,0) = 1.*(1.+1.i);
  m(0,1) = 2.*(1.+1.i);
  m(0,2) = 3.*(1.+1.i);
  m(1,0) = 4.*(1.+1.i);
  m(1,1) = 5.*(1.+1.i);
  m(1,2) = 6.*(1.+1.i);
  std::complex<double> * d = data(m);
  if (is_row_ordered(m)) {
    EXPECT_EQ(*(d+0), 1.*(1.+1.i));
    EXPECT_EQ(*(d+1), 2.*(1.+1.i));
    EXPECT_EQ(*(d+2), 3.*(1.+1.i));
    EXPECT_EQ(*(d+3), 4.*(1.+1.i));
    EXPECT_EQ(*(d+4), 5.*(1.+1.i));
    EXPECT_EQ(*(d+5), 6.*(1.+1.i));
  } else {
    EXPECT_EQ(*(d+0), 1.*(1.+1.i));
    EXPECT_EQ(*(d+1), 4.*(1.+1.i));
    EXPECT_EQ(*(d+2), 2.*(1.+1.i));
    EXPECT_EQ(*(d+3), 5.*(1.+1.i));
    EXPECT_EQ(*(d+4), 3.*(1.+1.i));
    EXPECT_EQ(*(d+5), 6.*(1.+1.i));
  }
}
