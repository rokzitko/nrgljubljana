// Explicit BLAS GEMM wrappers for row-major Eigen matrices.

#ifndef _blas_gemm_hpp_
#define _blas_gemm_hpp_

#include <cassert>
#include <complex>
#include <limits>
#include <stdexcept>

#include <fmt/format.h>

#include "linalg.hpp"
#include "traits.hpp"

namespace NRG {

enum class GemmTranspose { none, transpose, conjugate_transpose };

[[nodiscard]] inline lapack_int checked_blas_int(const size_t value, const char *what) {
  if (value > static_cast<size_t>(std::numeric_limits<lapack_int>::max()))
    throw std::runtime_error(fmt::format("{}={} exceeds BLAS integer range", what, value));
  return static_cast<lapack_int>(value);
}

[[nodiscard]] inline char gemm_transpose_char(const GemmTranspose trans) {
  switch (trans) {
  case GemmTranspose::none: return 'N';
  case GemmTranspose::transpose: return 'T';
  case GemmTranspose::conjugate_transpose: return 'C';
  }
  return 'N';
}

template<typename M>
[[nodiscard]] auto blas_outer_stride(const M &matrix, const char *name) {
  if (matrix.innerStride() != 1) throw std::runtime_error(fmt::format("{} has non-unit inner stride", name));
  return checked_blas_int(static_cast<size_t>(matrix.outerStride()), name);
}

template<typename T>
void call_gemm(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
               const T *alpha, const T *a, const lapack_int *lda, const T *b, const lapack_int *ldb, const T *beta, T *c,
               const lapack_int *ldc) {
  if constexpr (std::is_same_v<T, double>) {
    LAPACK_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  } else {
    LAPACK_zgemm(transa, transb, m, n, k, reinterpret_cast<const double *>(alpha), reinterpret_cast<const double *>(a), lda,
                 reinterpret_cast<const double *>(b), ldb, reinterpret_cast<const double *>(beta), reinterpret_cast<double *>(c), ldc);
  }
}

template<scalar T, Eigen_matrix AType, Eigen_matrix BType, Eigen_matrix CType>
void gemm_rowmajor(const T alpha, const AType &A, const GemmTranspose transA, const BType &B, const GemmTranspose transB, const T beta, CType &C) {
  const size_t a_rows = size1(A);
  const size_t a_cols = size2(A);
  const size_t b_rows = size1(B);
  const size_t b_cols = size2(B);
  const size_t m = transA == GemmTranspose::none ? a_rows : a_cols;
  const size_t k = transA == GemmTranspose::none ? a_cols : a_rows;
  const size_t b_k = transB == GemmTranspose::none ? b_rows : b_cols;
  const size_t n = transB == GemmTranspose::none ? b_cols : b_rows;

  if (k != b_k) throw std::runtime_error("GEMM inner dimensions do not match");
  assert(size1(C) == m && size2(C) == n);

  if (m == 0 || n == 0) return;
  if (k == 0) {
    if (beta == T{}) C.setZero();
    else C *= beta;
    return;
  }

  const lapack_int mm = checked_blas_int(m, "GEMM M");
  const lapack_int nn = checked_blas_int(n, "GEMM N");
  const lapack_int kk = checked_blas_int(k, "GEMM K");
  const lapack_int lda = blas_outer_stride(A, "GEMM A leading dimension");
  const lapack_int ldb = blas_outer_stride(B, "GEMM B leading dimension");
  const lapack_int ldc = blas_outer_stride(C, "GEMM C leading dimension");

  // For row-major storage, compute C^T = op(B)^T * op(A)^T using Fortran GEMM.
  const char transa = gemm_transpose_char(transB);
  const char transb = gemm_transpose_char(transA);
  call_gemm(&transa, &transb, &nn, &mm, &kk, &alpha, B.data(), &ldb, A.data(), &lda, &beta, C.data(), &ldc);
}

template<scalar T, Eigen_matrix AType, Eigen_matrix BType, Eigen_matrix CType>
void gemm_nn(const T alpha, const AType &A, const BType &B, const T beta, CType &C) {
  gemm_rowmajor(alpha, A, GemmTranspose::none, B, GemmTranspose::none, beta, C);
}

template<scalar T, Eigen_matrix AType, Eigen_matrix BType, Eigen_matrix CType>
void gemm_nh(const T alpha, const AType &A, const BType &B, const T beta, CType &C) {
  gemm_rowmajor(alpha, A, GemmTranspose::none, B, GemmTranspose::conjugate_transpose, beta, C);
}

template<scalar T, Eigen_matrix AType, Eigen_matrix BType, Eigen_matrix CType>
void gemm_hn(const T alpha, const AType &A, const BType &B, const T beta, CType &C) {
  gemm_rowmajor(alpha, A, GemmTranspose::conjugate_transpose, B, GemmTranspose::none, beta, C);
}

} // namespace NRG

#endif
