// diag.h - Diagonalisation code
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _diag_hpp_
#define _diag_hpp_

#include <type_traits> // is_same_v
#include <complex>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip> // std::setprecision
#include <stdexcept>
#include <limits>
#include <string>

#include "traits.hpp"
#include "params.hpp"
#include "eigen.hpp"
#include "time_mem.hpp"
#include "debug.hpp" // nrglogdp
#include "numerics.hpp" // is_matrix_upper

#include <fmt/format.h>

#include "lapack.h"

#ifndef NRG_ENABLE_CUDA
#define NRG_ENABLE_CUDA 0
#endif

#if NRG_ENABLE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

namespace NRG {

#if NRG_ENABLE_CUDA
inline void cuda_check(const cudaError_t status, const char *what) {
  if (status != cudaSuccess)
    throw std::runtime_error(fmt::format("{} failed: {}", what, cudaGetErrorString(status)));
}

inline void cusolver_check(const cusolverStatus_t status, const char *what) {
  if (status != CUSOLVER_STATUS_SUCCESS)
    throw std::runtime_error(fmt::format("{} failed. status={}", what, static_cast<int>(status)));
}

inline auto cuda_available_device_count() {
  int count = 0;
  const auto status = cudaGetDeviceCount(&count);
  if (status != cudaSuccess) {
    cudaGetLastError(); // clear the sticky error for callers that want to skip CUDA work
    return 0;
  }
  return count;
}

inline auto cuda_diag_requested(const std::string &diag) { return diag == "cuda_dsyevd" || diag == "cuda_zheevd"; }

inline void validate_cuda_diagonalisation_request(const std::string &diag) {
  if (!cuda_diag_requested(diag)) return;
  const auto count = cuda_available_device_count();
  if (count == 0) throw std::runtime_error(fmt::format("{} requested, but no CUDA accelerator was detected", diag));
}

class CudaSolverHandle {
 public:
  CudaSolverHandle() { cusolver_check(cusolverDnCreate(&handle_), "cusolverDnCreate"); }
  ~CudaSolverHandle() { if (handle_ != nullptr) cusolverDnDestroy(handle_); }
  CudaSolverHandle(const CudaSolverHandle &) = delete;
  CudaSolverHandle &operator=(const CudaSolverHandle &) = delete;
  operator cusolverDnHandle_t() const { return handle_; }

 private:
  cusolverDnHandle_t handle_{};
};

template<typename T> class CudaDeviceBuffer {
 public:
  explicit CudaDeviceBuffer(const size_t size) : size_(size) {
    if (size_ != 0) cuda_check(cudaMalloc(reinterpret_cast<void **>(&data_), size_ * sizeof(T)), "cudaMalloc");
  }
  ~CudaDeviceBuffer() { if (data_ != nullptr) cudaFree(data_); }
  CudaDeviceBuffer(const CudaDeviceBuffer &) = delete;
  CudaDeviceBuffer &operator=(const CudaDeviceBuffer &) = delete;
  T *data() { return data_; }
  const T *data() const { return data_; }
  size_t size() const { return size_; }

 private:
  T *data_{};
  size_t size_{};
};
#else
inline auto cuda_diag_requested(const std::string &diag) { return diag == "cuda_dsyevd" || diag == "cuda_zheevd"; }

inline void validate_cuda_diagonalisation_request(const std::string &diag) {
  if (cuda_diag_requested(diag)) throw std::runtime_error(fmt::format("{} requested, but CUDA support was not enabled at build time", diag));
}
#endif

[[nodiscard]] inline auto checked_lapack_int(const size_t value, const char *what) {
  if (value > static_cast<size_t>(std::numeric_limits<lapack_int>::max()))
    throw std::runtime_error(fmt::format("{}={} exceeds LAPACK integer range", what, value));
  return static_cast<lapack_int>(value);
}

template<vector SV, vector DV> requires std::is_convertible_v<typename SV::value_type, typename DV::value_type>
void copy_val(const SV &source, DV &dest, const size_t M) {
  using S = typename SV::value_type;
  my_assert(source.size() >= M);
  if (std::adjacent_find(source.begin(), source.begin() + M, std::greater<S>()) != source.begin() + M)
    std::cout << "WARNING: Values are not in ascending order. Bug in LAPACK dsyev* routines." << std::endl;
  dest.resize(M);
  std::copy_n(source.begin(), M, dest.begin());
}

template<typename U, matrix MM> // U may be _lapack_complex_double
void copy_vec(U* eigenvectors, MM & diagvectors, const size_t dim, const size_t M)
{
  diagvectors.resize(M, dim);
  for (const auto r : range0(M))
    for (const auto j : range0(dim)) diagvectors(r, j) = eigenvectors[dim * r + j];
  // this works correctly for both row-order and column-order storage types (diagvectors)
}

template<scalar S, vector V, typename U>
auto copy_results(const V &eigenvalues, U* eigenvectors, const char jobz, const size_t dim, const size_t M)
{
  RawEigen<S> d(M, dim);
  copy_val(eigenvalues, d.val, M);
  if (jobz == 'V') copy_vec(eigenvectors, d.vec, dim, M);
  my_assert(d.val.size() == nrvec(d.vec));
  return d;
}

// Perform diagonalisation: wrappers for LAPACK. jobz: 'N' for values only, 'V' for values and vectors
template<real_matrix RM>
auto diagonalise_dsyev(RM &m, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array
  double WORK0 = 0;         // on exit: optimal WORK size
  // Step 1: determine optimal LWORK
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK    = int(WORK0);
  std::vector<double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevd(RM &m, const char jobz = 'V')
{
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham       = data(m);
  std::vector<double> eigenvalues(dim);
  char UPLO  = 'L';
  int NN     = dim;
  int LDA    = dim;
  int INFO   = 0;
  int LWORK  = -1;
  int LIWORK = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  int IWORK0 = 0;   // on exit: optimal IWORK size
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK, &IWORK0, &LIWORK, &INFO);
  my_assert(INFO == 0);
  LWORK      = int(WORK0);
  LIWORK     = IWORK0;
  std::vector<double> WORK(LWORK);
  std::vector<int> IWORK(LIWORK);
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) {
    // dsyevd sometimes fails to converge (INFO>0). In such cases we do not trigger
    // an error but return 0, to permit error recovery.
    if (INFO > 0)
      return RawEigen<double>();
    else
      throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  }
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevr(RM &m, const double ratio = 1.0, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  // M is the number of the eigenvalues that we will attempt to
  // calculate using dsyevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = ceil(ratio * M); // round up
    M     = std::clamp<int>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';     // lower triangle of a is stored
  int NN        = dim;     // the order of the matrix
  int LDA       = dim;     // the leading dimension of the array a
  int INFO      = 0;       // 0 on successful exit
  double VL     = 0;       // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  int IL        = 1; // index range
  int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  int MM{}; // total number of eigenvalues found
  int LDZ = dim;
  std::vector<int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices
  //  indicating the nonzero elements in Z.  The i-th eigenvector is
  //  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<double> Z(LDZ * M); // eigenvectors
  int LWORK0  = -1;
  int LIWORK0 = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  int IWORK0 = 0;   // on exist: optimal IWORK size
  // Step 1: determine optimal LWORK and LIWORK
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK  = int(WORK0);
  int LIWORK = IWORK0;
  std::vector<double> WORK(LWORK);
  std::vector<int> IWORK(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(),
                &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  if (MM != int(M)) {
    std::cout << "dsyevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<double>(eigenvalues, Z.data(), jobz, dim, M);
}

template<complex_matrix CM>
auto diagonalise_zheev(CM &m, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  int RWORKdim = std::max(1, 3 * dim - 2);
  std::vector<double> RWORK(RWORKdim);
  // Step 1: determine optimal LWORK
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, RWORK.data(), &INFO);
  my_assert(INFO == 0);
  int LWORK  = int(WORK0.real());
  std::vector<lapack_complex_double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_zheevd(CM &m, const char jobz = 'V')
{
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  char UPLO  = 'L';
  int NN     = dim;
  int LDA    = dim;
  int INFO   = 0;
  int LWORK  = -1;
  int LRWORK = -1;
  int LIWORK = -1;
  lapack_complex_double WORK0{}; // on exit: optimal WORK size
  double RWORK0 = 0;             // on exit: optimal RWORK size
  int IWORK0 = 0;                // on exit: optimal IWORK size
  LAPACK_zheevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK, &RWORK0, &LRWORK, &IWORK0, &LIWORK, &INFO);
  my_assert(INFO == 0);
  LWORK  = int(WORK0.real());
  LRWORK = int(RWORK0);
  LIWORK = IWORK0;
  std::vector<lapack_complex_double> WORK(LWORK);
  std::vector<double> RWORK(LRWORK);
  std::vector<int> IWORK(LIWORK);
  LAPACK_zheevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, RWORK.data(), &LRWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) {
    if (INFO > 0)
      return RawEigen<std::complex<double>>();
    else
      throw std::runtime_error(fmt::format("zheevd failed. INFO={}", INFO));
  }
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_zheevr(CM &m, const double ratio = 1.0, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  // M is the number of the eigenvalues that we will attempt to
  // calculate using zheevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = ceil(ratio * M); // round up
    M     = std::clamp<int>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';      // lower triangle of a is stored
  int NN        = dim;      // the order of the matrix
  int LDA       = dim;      // the leading dimension of the array a
  int INFO      = 0;        // 0 on successful exit
  double VL     = 0;        // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  int IL        = 1; // index range
  int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  int MM = 0; // total number of eigenvalues found
  int LDZ = dim;
  std::vector<int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices indicating the nonzero elements in Z.  The i-th
  //  eigenvector is nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<lapack_complex_double> Z(LDZ * M); // eigenvectors
  int LWORK0 = -1;                 // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  int LRWORK0 = -1;  // query
  double RWORK0 = 0; // on exit: optimal RWORK size
  int LIWORK0 = -1;  // query
  int IWORK0 = 0;    // on exit: optimal IWORK size
  // Step 1: determine optimal LWORK, LRWORK, and LIWORK
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &RWORK0, &LRWORK0, &IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK   = int(WORK0.real());
  std::vector<lapack_complex_double> WORK(LWORK);
  int LRWORK  = int(RWORK0);
  std::vector<double> RWORK(LRWORK);
  int LIWORK  = IWORK0;
  std::vector<int> IWORK(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(), &LWORK,
                RWORK.data(), &LRWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheevr failed. INFO={}", INFO));
  if (MM != int(M)) {
    std::cout << "zheevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<std::complex<double>>(eigenvalues, Z.data(), jobz, dim, M);
}

#if NRG_ENABLE_CUDA
template<real_matrix RM>
auto diagonalise_cuda_dsyevd(RM &m, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  const auto elements = static_cast<size_t>(dim) * static_cast<size_t>(dim);
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  CudaSolverHandle solver;
  CudaDeviceBuffer<double> d_ham(elements);
  CudaDeviceBuffer<double> d_eigenvalues(dim);
  CudaDeviceBuffer<int> d_info(1);
  cuda_check(cudaMemcpy(d_ham.data(), ham, elements * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy H2D matrix");
  const auto mode = jobz == 'V' ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR;
  int lwork = 0;
  cusolver_check(cusolverDnDsyevd_bufferSize(solver, mode, CUBLAS_FILL_MODE_LOWER, dim, d_ham.data(), dim, d_eigenvalues.data(), &lwork),
                 "cusolverDnDsyevd_bufferSize");
  CudaDeviceBuffer<double> d_work(lwork);
  cusolver_check(cusolverDnDsyevd(solver, mode, CUBLAS_FILL_MODE_LOWER, dim, d_ham.data(), dim, d_eigenvalues.data(), d_work.data(), lwork, d_info.data()),
                 "cusolverDnDsyevd");
  int info = 0;
  cuda_check(cudaMemcpy(&info, d_info.data(), sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy D2H info");
  if (info != 0) {
    if (info > 0) return RawEigen<double>();
    throw std::runtime_error(fmt::format("cuda_dsyevd failed. INFO={}", info));
  }
  cuda_check(cudaMemcpy(eigenvalues.data(), d_eigenvalues.data(), eigenvalues.size() * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy D2H eigenvalues");
  if (jobz == 'V') cuda_check(cudaMemcpy(ham, d_ham.data(), elements * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy D2H eigenvectors");
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_cuda_zheevd(CM &m, const char jobz = 'V') {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  const auto elements = static_cast<size_t>(dim) * static_cast<size_t>(dim);
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  CudaSolverHandle solver;
  CudaDeviceBuffer<cuDoubleComplex> d_ham(elements);
  CudaDeviceBuffer<double> d_eigenvalues(dim);
  CudaDeviceBuffer<int> d_info(1);
  std::vector<cuDoubleComplex> h_ham(elements);
  for (const auto i : range0(elements)) h_ham[i] = make_cuDoubleComplex(ham[i].real(), ham[i].imag());
  cuda_check(cudaMemcpy(d_ham.data(), h_ham.data(), elements * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice), "cudaMemcpy H2D matrix");
  const auto mode = jobz == 'V' ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR;
  int lwork = 0;
  cusolver_check(cusolverDnZheevd_bufferSize(solver, mode, CUBLAS_FILL_MODE_LOWER, dim, d_ham.data(), dim, d_eigenvalues.data(), &lwork),
                 "cusolverDnZheevd_bufferSize");
  CudaDeviceBuffer<cuDoubleComplex> d_work(lwork);
  cusolver_check(cusolverDnZheevd(solver, mode, CUBLAS_FILL_MODE_LOWER, dim, d_ham.data(), dim, d_eigenvalues.data(), d_work.data(), lwork, d_info.data()),
                 "cusolverDnZheevd");
  int info = 0;
  cuda_check(cudaMemcpy(&info, d_info.data(), sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy D2H info");
  if (info != 0) {
    if (info > 0) return RawEigen<std::complex<double>>();
    throw std::runtime_error(fmt::format("cuda_zheevd failed. INFO={}", info));
  }
  cuda_check(cudaMemcpy(eigenvalues.data(), d_eigenvalues.data(), eigenvalues.size() * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy D2H eigenvalues");
  if (jobz == 'V') {
    cuda_check(cudaMemcpy(h_ham.data(), d_ham.data(), elements * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost), "cudaMemcpy D2H eigenvectors");
    for (const auto i : range0(elements)) ham[i] = std::complex<double>(cuCreal(h_ham[i]), cuCimag(h_ham[i]));
  }
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}
#endif

// Wrapper for the diagonalization of the Hamiltonian matrix. The number of eigenpairs returned does NOT need to be
// equal to the dimension of the matrix h. Matrix m is destroyed in the process, thus no const attribute!
template<matrix M> auto diagonalise(M &m, const DiagParams &DP, const int myrank) {
  using S = typename M::value_type;
  const std::string rank_string = myrank >= 0 ? " [rank=" + std::to_string(myrank) + "]" : "";
  mpilog("diagonalise " << size1(m) << "x" << size2(m) << " " << DP.diag << " " << DP.diagratio);
  nrglogdp('@', "diagonalise() - size(m)=" << size1(m) << rank_string);
  Timing timer;
  my_assert(is_matrix_upper(m));
  RawEigen<S> d;
  if constexpr (std::is_same_v<S, double>) {
    if (DP.diag == "dsyev"s) d = diagonalise_dsyev(m);
    if (DP.diag == "dsyevd"s  || DP.diag == "default"s) {
      d = diagonalise_dsyevd(m);
      if (d.getnrcomputed() == 0) {
        std::cout << "dsyevd failed, falling back to dsyev" << std::endl;
        d = diagonalise_dsyev(m);
      }
    }
    if (DP.diag == "dsyevr"s) d = diagonalise_dsyevr(m, DP.diagratio);
    if (DP.diag == "cuda_dsyevd"s) {
#if NRG_ENABLE_CUDA
      validate_cuda_diagonalisation_request(DP.diag);
      d = diagonalise_cuda_dsyevd(m);
      if (d.getnrcomputed() == 0) {
        std::cout << "cuda_dsyevd failed, falling back to dsyev" << std::endl;
        d = diagonalise_dsyev(m);
      }
#else
      throw std::runtime_error("cuda_dsyevd requested, but CUDA support was not enabled at build time");
#endif
    }
  }
  if constexpr (std::is_same_v<S, std::complex<double>>) {
    if (DP.diag == "zheev"s) d = diagonalise_zheev(m);
    if (DP.diag == "zheevd"s || DP.diag == "default"s) {
      d = diagonalise_zheevd(m);
      if (d.getnrcomputed() == 0) {
        std::cout << "zheevd failed, falling back to zheev" << std::endl;
        d = diagonalise_zheev(m);
      }
    }
    if (DP.diag == "zheevr"s) d = diagonalise_zheevr(m, DP.diagratio);
    if (DP.diag == "cuda_zheevd"s) {
#if NRG_ENABLE_CUDA
      validate_cuda_diagonalisation_request(DP.diag);
      d = diagonalise_cuda_zheevd(m);
      if (d.getnrcomputed() == 0) {
        std::cout << "cuda_zheevd failed, falling back to zheev" << std::endl;
        d = diagonalise_zheev(m);
      }
#else
      throw std::runtime_error("cuda_zheevd requested, but CUDA support was not enabled at build time");
#endif
    }
  }
  const auto nr_computed = d.getnrcomputed();
  my_assert(nr_computed > 0); // zero computed eigenvalues signals serious failure
  my_assert(has_lesseq_rows(d.vec, m)); // sanity check
  if (DP.logletter('e'))
    d.dump_eigenvalues();
  nrglogdp('A', "LAPACK, dim=" << dim(m) << " M=" << nr_computed << rank_string);
  nrglogdp('t', "Elapsed: " << std::setprecision(3) << timer.total_in_seconds() << rank_string);
  return d;
}

} // namespace

#endif
