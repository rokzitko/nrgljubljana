#include <gtest/gtest.h>

#define USE_EIGEN
#include <blas_gemm.hpp>
#include <numerics.hpp>
#include <traits.hpp>

using namespace NRG;
using namespace std::complex_literals;

template<scalar T>
T test_value(const size_t i, const size_t j) {
  if constexpr (std::is_same_v<T, double>) return 10.0 * static_cast<double>(i + 1) + static_cast<double>(j + 1);
  else return {10.0 * static_cast<double>(i + 1) + static_cast<double>(j + 1), static_cast<double>(i) - 2.0 * static_cast<double>(j)};
}

template<scalar T>
EigenMatrix<T> make_matrix(const size_t rows, const size_t cols) {
  EigenMatrix<T> matrix(rows, cols);
  for (const auto i : range0(rows))
    for (const auto j : range0(cols))
      matrix(i, j) = test_value<T>(i, j);
  return matrix;
}

template<scalar T>
void expect_near(const T actual, const T expected) {
  if constexpr (std::is_same_v<T, double>) {
    EXPECT_NEAR(actual, expected, 1e-10);
  } else {
    EXPECT_NEAR(actual.real(), expected.real(), 1e-10);
    EXPECT_NEAR(actual.imag(), expected.imag(), 1e-10);
  }
}

template<scalar T>
void expect_matrix_near(const EigenMatrix<T> &actual, const EigenMatrix<T> &expected) {
  ASSERT_EQ(size1(actual), size1(expected));
  ASSERT_EQ(size2(actual), size2(expected));
  for (const auto i : range0(size1(actual)))
    for (const auto j : range0(size2(actual)))
      expect_near(actual(i, j), expected(i, j));
}

template<scalar T>
T maybe_conj(const T value) {
  if constexpr (std::is_same_v<T, double>) return value;
  else return std::conj(value);
}

template<scalar T>
EigenMatrix<T> reference_nn(const T alpha, const EigenMatrix<T> &A, const EigenMatrix<T> &B, const T beta, const EigenMatrix<T> &C0) {
  EigenMatrix<T> result = C0;
  for (const auto i : range0(size1(result))) {
    for (const auto j : range0(size2(result))) {
      T sum{};
      for (const auto k : range0(size2(A))) sum += A(i, k) * B(k, j);
      result(i, j) = alpha * sum + beta * C0(i, j);
    }
  }
  return result;
}

template<scalar T>
EigenMatrix<T> reference_nh(const T alpha, const EigenMatrix<T> &A, const EigenMatrix<T> &B, const T beta, const EigenMatrix<T> &C0) {
  EigenMatrix<T> result = C0;
  for (const auto i : range0(size1(result))) {
    for (const auto j : range0(size2(result))) {
      T sum{};
      for (const auto k : range0(size2(A))) sum += A(i, k) * maybe_conj(B(j, k));
      result(i, j) = alpha * sum + beta * C0(i, j);
    }
  }
  return result;
}

template<scalar T>
EigenMatrix<T> reference_hn(const T alpha, const EigenMatrix<T> &A, const EigenMatrix<T> &B, const T beta, const EigenMatrix<T> &C0) {
  EigenMatrix<T> result = C0;
  for (const auto i : range0(size1(result))) {
    for (const auto j : range0(size2(result))) {
      T sum{};
      for (const auto k : range0(size1(A))) sum += maybe_conj(A(k, i)) * B(k, j);
      result(i, j) = alpha * sum + beta * C0(i, j);
    }
  }
  return result;
}

template<scalar T>
void test_gemm_nn() {
  const auto A = make_matrix<T>(2, 3);
  const auto B = make_matrix<T>(3, 4);
  auto C = make_matrix<T>(2, 4);
  const auto C0 = C;
  const T alpha = T{2.0};
  const T beta = T{-0.5};
  gemm_nn<T>(alpha, A, B, beta, C);
  expect_matrix_near(C, reference_nn(alpha, A, B, beta, C0));
}

template<scalar T>
void test_gemm_nh() {
  const auto A = make_matrix<T>(2, 3);
  const auto B = make_matrix<T>(4, 3);
  auto C = make_matrix<T>(2, 4);
  const auto C0 = C;
  const T alpha = T{1.5};
  const T beta = T{0.25};
  gemm_nh<T>(alpha, A, B, beta, C);
  expect_matrix_near(C, reference_nh(alpha, A, B, beta, C0));
}

template<scalar T>
void test_gemm_hn() {
  const auto A = make_matrix<T>(3, 2);
  const auto B = make_matrix<T>(3, 4);
  auto C = make_matrix<T>(2, 4);
  const auto C0 = C;
  const T alpha = T{-1.25};
  const T beta = T{0.75};
  gemm_hn<T>(alpha, A, B, beta, C);
  expect_matrix_near(C, reference_hn(alpha, A, B, beta, C0));
}

TEST(blas_gemm, real_nn) { test_gemm_nn<double>(); }
TEST(blas_gemm, real_nh) { test_gemm_nh<double>(); }
TEST(blas_gemm, real_hn) { test_gemm_hn<double>(); }

TEST(blas_gemm, complex_nn) { test_gemm_nn<std::complex<double>>(); }
TEST(blas_gemm, complex_nh) { test_gemm_nh<std::complex<double>>(); }
TEST(blas_gemm, complex_hn) { test_gemm_hn<std::complex<double>>(); }

TEST(blas_gemm, matrix_product_abstractions) {
  const auto A = make_matrix<double>(2, 3);
  const auto B = make_matrix<double>(3, 4);
  const EigenMatrix<double> C0 = EigenMatrix<double>::Zero(2, 4);
  expect_matrix_near(matrix_prod<double>(A, B), reference_nn(1.0, A, B, 0.0, C0));
}

TEST(blas_gemm, matrix_adjoint_product_abstraction) {
  const auto A = make_matrix<std::complex<double>>(3, 2);
  const auto B = make_matrix<std::complex<double>>(3, 4);
  const EigenMatrix<std::complex<double>> C0 = EigenMatrix<std::complex<double>>::Zero(2, 4);
  expect_matrix_near(matrix_adj_prod<std::complex<double>>(A, B), reference_hn(1.0 + 0.0i, A, B, 0.0 + 0.0i, C0));
}

TEST(blas_gemm, product_transform_rotate_abstractions) {
  const auto A = make_matrix<double>(2, 3);
  const auto B = make_matrix<double>(4, 3);
  auto M = make_matrix<double>(2, 4);
  const auto M0 = M;
  product_Eigen<double>(M, 2.0, A, B);
  expect_matrix_near(M, reference_nh(2.0, A, B, 1.0, M0));

  const auto O = make_matrix<double>(3, 5);
  const auto C = make_matrix<double>(4, 5);
  auto T = make_matrix<double>(2, 4);
  const auto T0 = T;
  const EigenMatrix<double> AO0 = EigenMatrix<double>::Zero(2, 5);
  const auto AO = reference_nn(1.0, A, O, 0.0, AO0);
  transform_Eigen<double>(T, -0.5, A, O, C);
  expect_matrix_near(T, reference_nh(-0.5, AO, C, 1.0, T0));

  const auto U = make_matrix<double>(3, 2);
  const auto Rop = make_matrix<double>(3, 3);
  auto R = make_matrix<double>(2, 2);
  const auto R0 = R;
  const EigenMatrix<double> UO0 = EigenMatrix<double>::Zero(2, 3);
  const auto UO = reference_hn(1.0, U, Rop, 0.0, UO0);
  rotate<double>(R, 1.25, U, Rop);
  expect_matrix_near(R, reference_nn(1.25, UO, U, 1.0, R0));
}

TEST(blas_gemm, rotate_diagonal_abstraction) {
  const auto U = make_matrix<std::complex<double>>(3, 2);
  EigenMatrix<std::complex<double>> O = EigenMatrix<std::complex<double>>::Zero(3, 3);
  O(0, 0) = 2.0 + 1.0i;
  O(1, 1) = -0.5 + 0.25i;
  O(2, 2) = 1.5 - 0.75i;
  auto M = make_matrix<std::complex<double>>(2, 2);
  const auto M0 = M;
  EigenMatrix<std::complex<double>> scaled(3, 2);
  for (const auto i : range0(size1(U)))
    for (const auto j : range0(size2(U)))
      scaled(i, j) = O(i, i) * U(i, j);
  rotate_diagonal<std::complex<double>>(M, 0.5 - 0.25i, U, O);
  expect_matrix_near(M, reference_hn(0.5 - 0.25i, U, scaled, 1.0 + 0.0i, M0));
}
