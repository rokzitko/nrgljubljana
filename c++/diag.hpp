// diag.h - Diagonalisation code
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _diag_hpp_
#define _diag_hpp_

#include <type_traits> // is_same_v
#include <algorithm>
#include <complex>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip> // std::setprecision
#include <stdexcept>
#include <limits>
#include <string>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "traits.hpp"
#include "params.hpp"
#include "eigen.hpp"
#include "time_mem.hpp"
#include "debug.hpp" // nrglogdp
#include "numerics.hpp" // is_matrix_upper

#include <fmt/format.h>

#include "linalg.hpp"

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

inline auto cuda_diag_requested(const std::string &diag) { return diag == "cuda" || diag == "cuda_dsyevd" || diag == "cuda_zheevd"; }

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

class CudaSolverParams {
 public:
  CudaSolverParams() { cusolver_check(cusolverDnCreateParams(&params_), "cusolverDnCreateParams"); }
  ~CudaSolverParams() { if (params_ != nullptr) cusolverDnDestroyParams(params_); }
  CudaSolverParams(const CudaSolverParams &) = delete;
  CudaSolverParams &operator=(const CudaSolverParams &) = delete;
  operator cusolverDnParams_t() const { return params_; }

 private:
  cusolverDnParams_t params_{};
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
inline auto cuda_diag_requested(const std::string &diag) { return diag == "cuda" || diag == "cuda_dsyevd" || diag == "cuda_zheevd"; }

inline void validate_cuda_diagonalisation_request(const std::string &diag) {
  if (cuda_diag_requested(diag)) throw std::runtime_error(fmt::format("{} requested, but CUDA support was not enabled at build time", diag));
}
#endif

[[nodiscard]] inline auto checked_lapack_int(const size_t value, const char *what) {
  if (value > static_cast<size_t>(std::numeric_limits<lapack_int>::max()))
    throw std::runtime_error(fmt::format("{}={} exceeds LAPACK integer range", what, value));
  return static_cast<lapack_int>(value);
}

[[nodiscard]] inline auto checked_workspace_query(const double value, const char *routine, const char *name) {
  if (!std::isfinite(value) || value < 1.0) throw std::runtime_error(fmt::format("{} returned invalid {} workspace size {}", routine, name, value));
  if (value > static_cast<double>(std::numeric_limits<size_t>::max())) return std::numeric_limits<size_t>::max();
  return static_cast<size_t>(value);
}

[[nodiscard]] inline auto checked_workspace_formula(const std::uint64_t value, const char *routine, const char *name) {
  if (value < 1) throw std::runtime_error(fmt::format("{} computed invalid {} workspace size", routine, name));
  if (value > static_cast<std::uint64_t>(std::numeric_limits<size_t>::max()))
    throw std::runtime_error(fmt::format("{} computed {} workspace size larger than size_t", routine, name));
  return static_cast<size_t>(value);
}

[[nodiscard]] inline auto max_one(const std::uint64_t value) { return value > 1 ? value : std::uint64_t{1}; }

template<typename T> [[nodiscard]] auto workspace_cap() {
  return std::min(std::vector<T>().max_size(), static_cast<size_t>(std::numeric_limits<lapack_int>::max()));
}

template<typename T>
[[nodiscard]] auto select_workspace_size(const char *routine, const char *name, const size_t minimum, const size_t optimal, const bool saveram) {
  const auto cap = workspace_cap<T>();
  if (minimum > cap)
    throw std::runtime_error(fmt::format("{} minimum {} workspace size {} exceeds maximum usable size {}", routine, name, minimum, cap));
  const auto selected = saveram ? minimum : std::min(std::max(optimal, minimum), cap);
  return checked_lapack_int(selected, fmt::format("{} {} workspace size", routine, name).c_str());
}

template<typename T, typename I> [[nodiscard]] auto workspace_bytes(const I count) {
  return static_cast<size_t>(count) * sizeof(T);
}

[[nodiscard]] inline auto dsyev_min_lwork(const size_t n) {
  if (n == 0) return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(max_one(3 * nn - 1), "dsyev", "LWORK");
}

[[nodiscard]] inline auto dsyevd_min_lwork(const size_t n, const char jobz) {
  if (n <= 1) return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(jobz == 'N' ? 2 * nn + 1 : 1 + 6 * nn + 2 * nn * nn, "dsyevd", "LWORK");
}

[[nodiscard]] inline auto dsyevd_min_liwork(const size_t n, const char jobz) {
  if (n <= 1 || jobz == 'N') return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(3 + 5 * nn, "dsyevd", "LIWORK");
}

[[nodiscard]] inline auto dsyevr_min_lwork(const size_t n) { return checked_workspace_formula(max_one(26 * static_cast<std::uint64_t>(n)), "dsyevr", "LWORK"); }
[[nodiscard]] inline auto dsyevr_min_liwork(const size_t n) { return checked_workspace_formula(max_one(10 * static_cast<std::uint64_t>(n)), "dsyevr", "LIWORK"); }
[[nodiscard]] inline auto zheev_min_lwork(const size_t n) {
  if (n == 0) return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(max_one(2 * nn - 1), "zheev", "LWORK");
}
[[nodiscard]] inline auto zheev_min_rwork(const size_t n) {
  if (n == 0) return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(max_one(3 * nn - 2), "zheev", "RWORK");
}

[[nodiscard]] inline auto zheevd_min_lwork(const size_t n, const char jobz) {
  if (n <= 1) return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(jobz == 'N' ? nn + 1 : 2 * nn + nn * nn, "zheevd", "LWORK");
}

[[nodiscard]] inline auto zheevd_min_lrwork(const size_t n, const char jobz) {
  if (n <= 1) return size_t(1);
  if (jobz == 'N') return n;
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(1 + 5 * nn + 2 * nn * nn, "zheevd", "LRWORK");
}

[[nodiscard]] inline auto zheevd_min_liwork(const size_t n, const char jobz) {
  if (n <= 1 || jobz == 'N') return size_t(1);
  const auto nn = static_cast<std::uint64_t>(n);
  return checked_workspace_formula(3 + 5 * nn, "zheevd", "LIWORK");
}

[[nodiscard]] inline auto zheevr_min_lwork(const size_t n) { return checked_workspace_formula(max_one(2 * static_cast<std::uint64_t>(n)), "zheevr", "LWORK"); }
[[nodiscard]] inline auto zheevr_min_lrwork(const size_t n) { return checked_workspace_formula(max_one(24 * static_cast<std::uint64_t>(n)), "zheevr", "LRWORK"); }
[[nodiscard]] inline auto zheevr_min_liwork(const size_t n) { return checked_workspace_formula(max_one(10 * static_cast<std::uint64_t>(n)), "zheevr", "LIWORK"); }

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

template<vector V, typename U>
void sort_eigenpairs_by_value(V &eigenvalues, U *eigenvectors, const char jobz, const size_t dim, const size_t M) {
  my_assert(eigenvalues.size() >= M);
  if (M < 2) return;
  std::vector<size_t> p(M);
  std::iota(p.begin(), p.end(), 0);
  std::stable_sort(p.begin(), p.end(), [&eigenvalues](const auto a, const auto b) { return eigenvalues[a] < eigenvalues[b]; });
  const auto already_sorted = ranges::equal(p, range0(M));
  if (already_sorted) return;

  std::vector<typename V::value_type> eigenvalues_sorted(M);
  for (const auto r : range0(M)) eigenvalues_sorted[r] = std::move(eigenvalues[p[r]]);
  std::move(eigenvalues_sorted.begin(), eigenvalues_sorted.end(), eigenvalues.begin());

  if (jobz != 'V') return;
  std::vector<U> eigenvectors_sorted(M * dim);
  for (const auto r : range0(M))
    for (const auto j : range0(dim)) eigenvectors_sorted[dim * r + j] = std::move(eigenvectors[dim * p[r] + j]);
  std::move(eigenvectors_sorted.begin(), eigenvectors_sorted.end(), eigenvectors);
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
auto diagonalise_dsyev(RM &m, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  lapack_int NN     = dim;         // the order of the matrix
  lapack_int LDA    = dim;         // the leading dimension of the array a
  lapack_int INFO   = 0;           // 0 on successful exit
  lapack_int LWORK0 = -1;          // length of the WORK array
  double WORK0 = 0;         // on exit: optimal WORK size
  // Step 1: determine optimal LWORK
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev workspace query failed. INFO={}", INFO));
  const auto min_lwork = dsyev_min_lwork(dim);
  const auto opt_lwork = checked_workspace_query(WORK0, "dsyev", "LWORK");
  auto LWORK = select_workspace_size<double>("dsyev", "LWORK", min_lwork, opt_lwork, saveram);
  if (log_workspace) std::cout << "dsyev workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK)=" << min_lwork << " opt(LWORK)=" << opt_lwork << " use(LWORK)=" << LWORK << " workspace_bytes=" << format_memory_usage(workspace_bytes<double>(LWORK)) << std::endl;
  std::vector<double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevd(RM &m, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false)
{
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham       = data(m);
  std::vector<double> eigenvalues(dim);
  char UPLO  = 'L';
  lapack_int NN     = dim;
  lapack_int LDA    = dim;
  lapack_int INFO   = 0;
  lapack_int LWORK  = -1;
  lapack_int LIWORK = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  lapack_int IWORK0 = 0;   // on exit: optimal IWORK size
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK, &IWORK0, &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyevd workspace query failed. INFO={}", INFO));
  const auto min_lwork  = dsyevd_min_lwork(dim, jobz);
  const auto min_liwork = dsyevd_min_liwork(dim, jobz);
  const auto opt_lwork  = checked_workspace_query(WORK0, "dsyevd", "LWORK");
  const auto opt_liwork = checked_workspace_query(static_cast<double>(IWORK0), "dsyevd", "LIWORK");
  LWORK  = select_workspace_size<double>("dsyevd", "LWORK", min_lwork, opt_lwork, saveram);
  LIWORK = select_workspace_size<lapack_int>("dsyevd", "LIWORK", min_liwork, opt_liwork, saveram);
  if (log_workspace) std::cout << "dsyevd workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK,LIWORK)=" << min_lwork << "," << min_liwork << " opt(LWORK,LIWORK)=" << opt_lwork << "," << opt_liwork << " use(LWORK,LIWORK)=" << LWORK << "," << LIWORK << " workspace_bytes=" << format_memory_usage(workspace_bytes<double>(LWORK) + workspace_bytes<lapack_int>(LIWORK)) << std::endl;
  std::vector<double> WORK(LWORK);
  std::vector<lapack_int> IWORK(LIWORK);
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0)
    throw std::runtime_error(fmt::format("dsyevd failed. INFO={}. LAPACK may have overwritten the input matrix; not retrying fallback.", INFO));
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevr(RM &m, const double ratio = 1.0, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  // M is the number of the eigenvalues that we will attempt to
  // calculate using dsyevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = ceil(ratio * M); // round up
    M     = std::clamp<lapack_int>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';     // lower triangle of a is stored
  lapack_int NN        = dim;     // the order of the matrix
  lapack_int LDA       = dim;     // the leading dimension of the array a
  lapack_int INFO      = 0;       // 0 on successful exit
  double VL     = 0;       // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  lapack_int IL        = 1; // index range
  lapack_int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  lapack_int MM{}; // total number of eigenvalues found
  lapack_int LDZ = dim;
  std::vector<lapack_int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices
  //  indicating the nonzero elements in Z.  The i-th eigenvector is
  //  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<double> Z(LDZ * M); // eigenvectors
  lapack_int LWORK0  = -1;
  lapack_int LIWORK0 = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  lapack_int IWORK0 = 0;   // on exist: optimal IWORK size
  // Step 1: determine optimal LWORK and LIWORK
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &IWORK0, &LIWORK0, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyevr workspace query failed. INFO={}", INFO));
  const auto min_lwork  = dsyevr_min_lwork(dim);
  const auto min_liwork = dsyevr_min_liwork(dim);
  const auto opt_lwork  = checked_workspace_query(WORK0, "dsyevr", "LWORK");
  const auto opt_liwork = checked_workspace_query(static_cast<double>(IWORK0), "dsyevr", "LIWORK");
  auto LWORK  = select_workspace_size<double>("dsyevr", "LWORK", min_lwork, opt_lwork, saveram);
  auto LIWORK = select_workspace_size<lapack_int>("dsyevr", "LIWORK", min_liwork, opt_liwork, saveram);
  if (log_workspace) std::cout << "dsyevr workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK,LIWORK)=" << min_lwork << "," << min_liwork << " opt(LWORK,LIWORK)=" << opt_lwork << "," << opt_liwork << " use(LWORK,LIWORK)=" << LWORK << "," << LIWORK << " workspace_bytes=" << format_memory_usage(workspace_bytes<double>(LWORK) + workspace_bytes<lapack_int>(LIWORK)) << std::endl;
  std::vector<double> WORK(LWORK);
  std::vector<lapack_int> IWORK(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(),
                &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  if (MM != M) {
    std::cout << "dsyevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<double>(eigenvalues, Z.data(), jobz, dim, M);
}

template<complex_matrix CM>
auto diagonalise_zheev(CM &m, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  lapack_int NN     = dim;         // the order of the matrix
  lapack_int LDA    = dim;         // the leading dimension of the array a
  lapack_int INFO   = 0;           // 0 on successful exit
  lapack_int LWORK0 = -1;          // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  auto RWORKdim = checked_lapack_int(zheev_min_rwork(dim), "zheev RWORK workspace size");
  std::vector<double> RWORK(RWORKdim);
  // Step 1: determine optimal LWORK
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, RWORK.data(), &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheev workspace query failed. INFO={}", INFO));
  const auto min_lwork = zheev_min_lwork(dim);
  const auto opt_lwork = checked_workspace_query(WORK0.real(), "zheev", "LWORK");
  auto LWORK = select_workspace_size<lapack_complex_double>("zheev", "LWORK", min_lwork, opt_lwork, saveram);
  if (log_workspace) std::cout << "zheev workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK,RWORK)=" << min_lwork << "," << RWORKdim << " opt(LWORK)=" << opt_lwork << " use(LWORK,RWORK)=" << LWORK << "," << RWORKdim << " workspace_bytes=" << format_memory_usage(workspace_bytes<lapack_complex_double>(LWORK) + workspace_bytes<double>(RWORKdim)) << std::endl;
  std::vector<lapack_complex_double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_zheevd(CM &m, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false)
{
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  char UPLO  = 'L';
  lapack_int NN     = dim;
  lapack_int LDA    = dim;
  lapack_int INFO   = 0;
  lapack_int LWORK  = -1;
  lapack_int LRWORK = -1;
  lapack_int LIWORK = -1;
  lapack_complex_double WORK0{}; // on exit: optimal WORK size
  double RWORK0 = 0;             // on exit: optimal RWORK size
  lapack_int IWORK0 = 0;                // on exit: optimal IWORK size
  LAPACK_zheevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK, &RWORK0, &LRWORK, &IWORK0, &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheevd workspace query failed. INFO={}", INFO));
  const auto min_lwork  = zheevd_min_lwork(dim, jobz);
  const auto min_lrwork = zheevd_min_lrwork(dim, jobz);
  const auto min_liwork = zheevd_min_liwork(dim, jobz);
  const auto opt_lwork  = checked_workspace_query(WORK0.real(), "zheevd", "LWORK");
  const auto opt_lrwork = checked_workspace_query(RWORK0, "zheevd", "LRWORK");
  const auto opt_liwork = checked_workspace_query(static_cast<double>(IWORK0), "zheevd", "LIWORK");
  LWORK  = select_workspace_size<lapack_complex_double>("zheevd", "LWORK", min_lwork, opt_lwork, saveram);
  LRWORK = select_workspace_size<double>("zheevd", "LRWORK", min_lrwork, opt_lrwork, saveram);
  LIWORK = select_workspace_size<lapack_int>("zheevd", "LIWORK", min_liwork, opt_liwork, saveram);
  if (log_workspace) std::cout << "zheevd workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK,LRWORK,LIWORK)=" << min_lwork << "," << min_lrwork << "," << min_liwork << " opt(LWORK,LRWORK,LIWORK)=" << opt_lwork << "," << opt_lrwork << "," << opt_liwork << " use(LWORK,LRWORK,LIWORK)=" << LWORK << "," << LRWORK << "," << LIWORK << " workspace_bytes=" << format_memory_usage(workspace_bytes<lapack_complex_double>(LWORK) + workspace_bytes<double>(LRWORK) + workspace_bytes<lapack_int>(LIWORK)) << std::endl;
  std::vector<lapack_complex_double> WORK(LWORK);
  std::vector<double> RWORK(LRWORK);
  std::vector<lapack_int> IWORK(LIWORK);
  LAPACK_zheevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, RWORK.data(), &LRWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0)
    throw std::runtime_error(fmt::format("zheevd failed. INFO={}. LAPACK may have overwritten the input matrix; not retrying fallback.", INFO));
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_zheevr(CM &m, const double ratio = 1.0, const char jobz = 'V', const bool saveram = false, const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = checked_lapack_int(size1(m), "matrix dimension");
  // M is the number of the eigenvalues that we will attempt to
  // calculate using zheevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = ceil(ratio * M); // round up
    M     = std::clamp<lapack_int>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';      // lower triangle of a is stored
  lapack_int NN        = dim;      // the order of the matrix
  lapack_int LDA       = dim;      // the leading dimension of the array a
  lapack_int INFO      = 0;        // 0 on successful exit
  double VL     = 0;        // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  lapack_int IL        = 1; // index range
  lapack_int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  lapack_int MM = 0; // total number of eigenvalues found
  lapack_int LDZ = dim;
  std::vector<lapack_int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices indicating the nonzero elements in Z.  The i-th
  //  eigenvector is nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<lapack_complex_double> Z(LDZ * M); // eigenvectors
  lapack_int LWORK0 = -1;                 // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  lapack_int LRWORK0 = -1;  // query
  double RWORK0 = 0; // on exit: optimal RWORK size
  lapack_int LIWORK0 = -1;  // query
  lapack_int IWORK0 = 0;    // on exit: optimal IWORK size
  // Step 1: determine optimal LWORK, LRWORK, and LIWORK
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &RWORK0, &LRWORK0, &IWORK0, &LIWORK0, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheevr workspace query failed. INFO={}", INFO));
  const auto min_lwork  = zheevr_min_lwork(dim);
  const auto min_lrwork = zheevr_min_lrwork(dim);
  const auto min_liwork = zheevr_min_liwork(dim);
  const auto opt_lwork  = checked_workspace_query(WORK0.real(), "zheevr", "LWORK");
  const auto opt_lrwork = checked_workspace_query(RWORK0, "zheevr", "LRWORK");
  const auto opt_liwork = checked_workspace_query(static_cast<double>(IWORK0), "zheevr", "LIWORK");
  auto LWORK   = select_workspace_size<lapack_complex_double>("zheevr", "LWORK", min_lwork, opt_lwork, saveram);
  std::vector<lapack_complex_double> WORK(LWORK);
  auto LRWORK  = select_workspace_size<double>("zheevr", "LRWORK", min_lrwork, opt_lrwork, saveram);
  std::vector<double> RWORK(LRWORK);
  auto LIWORK  = select_workspace_size<lapack_int>("zheevr", "LIWORK", min_liwork, opt_liwork, saveram);
  std::vector<lapack_int> IWORK(LIWORK);
  if (log_workspace) std::cout << "zheevr workspace dim=" << dim << " jobz=" << jobz << " saveram=" << saveram << " min(LWORK,LRWORK,LIWORK)=" << min_lwork << "," << min_lrwork << "," << min_liwork << " opt(LWORK,LRWORK,LIWORK)=" << opt_lwork << "," << opt_lrwork << "," << opt_liwork << " use(LWORK,LRWORK,LIWORK)=" << LWORK << "," << LRWORK << "," << LIWORK << " workspace_bytes=" << format_memory_usage(workspace_bytes<lapack_complex_double>(LWORK) + workspace_bytes<double>(LRWORK) + workspace_bytes<lapack_int>(LIWORK)) << std::endl;
  // Step 2: perform the diagonalisation
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(), &LWORK,
                RWORK.data(), &LRWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheevr failed. INFO={}", INFO));
  if (MM != M) {
    std::cout << "zheevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<std::complex<double>>(eigenvalues, Z.data(), jobz, dim, M);
}

#if NRG_ENABLE_CUDA
template<real_matrix RM>
auto diagonalise_cuda_dsyevd(RM &m, const char jobz = 'V', const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = size1(m);
  const auto dim64 = static_cast<int64_t>(dim);
  const auto elements = static_cast<size_t>(dim) * static_cast<size_t>(dim);
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  CudaSolverHandle solver;
  CudaSolverParams params;
  CudaDeviceBuffer<double> d_ham(elements);
  CudaDeviceBuffer<double> d_eigenvalues(dim);
  CudaDeviceBuffer<int> d_info(1);
  cuda_check(cudaMemcpy(d_ham.data(), ham, elements * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy H2D matrix");
  const auto mode = jobz == 'V' ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR;
  size_t device_workspace_bytes = 0;
  size_t host_workspace_bytes = 0;
  cusolver_check(cusolverDnXsyevd_bufferSize(solver, params, mode, CUBLAS_FILL_MODE_LOWER, dim64, CUDA_R_64F, d_ham.data(), dim64,
                                             CUDA_R_64F, d_eigenvalues.data(), CUDA_R_64F, &device_workspace_bytes, &host_workspace_bytes),
                 "cusolverDnXsyevd_bufferSize");
  if (log_workspace) std::cout << "cuda_dsyevd workspace dim=" << dim << " jobz=" << jobz << " device_bytes=" << format_memory_usage(device_workspace_bytes) << " host_bytes=" << format_memory_usage(host_workspace_bytes) << std::endl;
  CudaDeviceBuffer<char> d_work(device_workspace_bytes);
  std::vector<char> h_work(host_workspace_bytes);
  cusolver_check(cusolverDnXsyevd(solver, params, mode, CUBLAS_FILL_MODE_LOWER, dim64, CUDA_R_64F, d_ham.data(), dim64,
                                  CUDA_R_64F, d_eigenvalues.data(), CUDA_R_64F, d_work.data(), device_workspace_bytes,
                                  h_work.data(), host_workspace_bytes, d_info.data()),
                 "cusolverDnXsyevd");
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
auto diagonalise_cuda_zheevd(CM &m, const char jobz = 'V', const bool log_workspace = false) {
  if (!is_row_ordered(m)) m = NRG::trans(m);
  const auto dim = size1(m);
  const auto dim64 = static_cast<int64_t>(dim);
  const auto elements = static_cast<size_t>(dim) * static_cast<size_t>(dim);
  auto ham = data(m);
  std::vector<double> eigenvalues(dim);
  CudaSolverHandle solver;
  CudaSolverParams params;
  CudaDeviceBuffer<cuDoubleComplex> d_ham(elements);
  CudaDeviceBuffer<double> d_eigenvalues(dim);
  CudaDeviceBuffer<int> d_info(1);
  std::vector<cuDoubleComplex> h_ham(elements);
  for (const auto i : range0(elements)) h_ham[i] = make_cuDoubleComplex(ham[i].real(), ham[i].imag());
  cuda_check(cudaMemcpy(d_ham.data(), h_ham.data(), elements * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice), "cudaMemcpy H2D matrix");
  const auto mode = jobz == 'V' ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR;
  size_t device_workspace_bytes = 0;
  size_t host_workspace_bytes = 0;
  cusolver_check(cusolverDnXsyevd_bufferSize(solver, params, mode, CUBLAS_FILL_MODE_LOWER, dim64, CUDA_C_64F, d_ham.data(), dim64,
                                             CUDA_R_64F, d_eigenvalues.data(), CUDA_C_64F, &device_workspace_bytes, &host_workspace_bytes),
                 "cusolverDnXsyevd_bufferSize");
  if (log_workspace) std::cout << "cuda_zheevd workspace dim=" << dim << " jobz=" << jobz << " device_bytes=" << format_memory_usage(device_workspace_bytes) << " host_bytes=" << format_memory_usage(host_workspace_bytes) << std::endl;
  CudaDeviceBuffer<char> d_work(device_workspace_bytes);
  std::vector<char> h_work(host_workspace_bytes);
  cusolver_check(cusolverDnXsyevd(solver, params, mode, CUBLAS_FILL_MODE_LOWER, dim64, CUDA_C_64F, d_ham.data(), dim64,
                                  CUDA_R_64F, d_eigenvalues.data(), CUDA_C_64F, d_work.data(), device_workspace_bytes,
                                  h_work.data(), host_workspace_bytes, d_info.data()),
                 "cusolverDnXsyevd");
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
  sort_eigenpairs_by_value(eigenvalues, ham, jobz, dim, dim);
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
  const auto log_workspace = DP.logletter('#');
  Timing timer;
  my_assert(is_matrix_upper(m));
  RawEigen<S> d;
  if constexpr (std::is_same_v<S, double>) {
    if (DP.diag == "dsyev"s) d = diagonalise_dsyev(m, 'V', DP.saveram, log_workspace);
    if (DP.diag == "dsyevd"s  || DP.diag == "default"s) {
      d = diagonalise_dsyevd(m, 'V', DP.saveram, log_workspace);
    }
    if (DP.diag == "dsyevr"s) d = diagonalise_dsyevr(m, DP.diagratio, 'V', DP.saveram, log_workspace);
    if (DP.diag == "cuda"s || DP.diag == "cuda_dsyevd"s) {
#if NRG_ENABLE_CUDA
      validate_cuda_diagonalisation_request(DP.diag);
      d = diagonalise_cuda_dsyevd(m, 'V', log_workspace);
      if (d.getnrcomputed() == 0) {
        std::cout << "cuda_dsyevd failed, falling back to dsyev" << std::endl;
        d = diagonalise_dsyev(m, 'V', DP.saveram, log_workspace);
      }
#else
      throw std::runtime_error(fmt::format("{} requested, but CUDA support was not enabled at build time", DP.diag));
#endif
    }
  }
  if constexpr (std::is_same_v<S, std::complex<double>>) {
    if (DP.diag == "zheev"s) d = diagonalise_zheev(m, 'V', DP.saveram, log_workspace);
    if (DP.diag == "zheevd"s || DP.diag == "default"s) {
      d = diagonalise_zheevd(m, 'V', DP.saveram, log_workspace);
    }
    if (DP.diag == "zheevr"s) d = diagonalise_zheevr(m, DP.diagratio, 'V', DP.saveram, log_workspace);
    if (DP.diag == "cuda"s || DP.diag == "cuda_zheevd"s) {
#if NRG_ENABLE_CUDA
      validate_cuda_diagonalisation_request(DP.diag);
      d = diagonalise_cuda_zheevd(m, 'V', log_workspace);
      if (d.getnrcomputed() == 0) {
        std::cout << "cuda_zheevd failed, falling back to zheev" << std::endl;
        d = diagonalise_zheev(m, 'V', DP.saveram, log_workspace);
      }
#else
      throw std::runtime_error(fmt::format("{} requested, but CUDA support was not enabled at build time", DP.diag));
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
