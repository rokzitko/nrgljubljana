// Don't include directly; include numerics.hpp.

#ifndef _NUMERICS_CUDA_HPP_
#define _NUMERICS_CUDA_HPP_

#ifndef NRG_ENABLE_CUDA
#define NRG_ENABLE_CUDA 0
#endif

#if NRG_ENABLE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#endif

inline auto cuda_mult_requested(const std::string &mult) { return mult == "cuda"; }

#if NRG_ENABLE_CUDA
inline void cuda_mult_check(const cudaError_t status, const char *what) {
  if (status != cudaSuccess)
    throw std::runtime_error(fmt::format("{} failed: {}", what, cudaGetErrorString(status)));
}

inline void cublas_check(const cublasStatus_t status, const char *what) {
  if (status != CUBLAS_STATUS_SUCCESS)
    throw std::runtime_error(fmt::format("{} failed. status={}", what, static_cast<int>(status)));
}

inline auto cuda_mult_available_device_count() {
  int count = 0;
  const auto status = cudaGetDeviceCount(&count);
  if (status != cudaSuccess) {
    cudaGetLastError();
    return 0;
  }
  return count;
}

inline void validate_cuda_multiplication_request(const std::string &mult) {
  if (!cuda_mult_requested(mult)) return;
  if (cuda_mult_available_device_count() == 0) throw std::runtime_error("mult=cuda requested, but no CUDA accelerator was detected");
}

class CublasHandle {
 public:
  CublasHandle() { cublas_check(cublasCreate(&handle_), "cublasCreate"); }
  ~CublasHandle() { if (handle_ != nullptr) cublasDestroy(handle_); }
  CublasHandle(const CublasHandle &) = delete;
  CublasHandle &operator=(const CublasHandle &) = delete;
  operator cublasHandle_t() const { return handle_; }

 private:
  cublasHandle_t handle_{};
};

template<typename T> class CudaMultBuffer {
 public:
  explicit CudaMultBuffer(const size_t size) : size_(size) {
    if (size_ != 0) cuda_mult_check(cudaMalloc(reinterpret_cast<void **>(&data_), size_ * sizeof(T)), "cudaMalloc");
  }
  ~CudaMultBuffer() { if (data_ != nullptr) cudaFree(data_); }
  CudaMultBuffer(const CudaMultBuffer &) = delete;
  CudaMultBuffer &operator=(const CudaMultBuffer &) = delete;
  T *data() { return data_; }
  const T *data() const { return data_; }
  size_t bytes() const { return size_ * sizeof(T); }

 private:
  T *data_{};
  size_t size_{};
};

template<scalar S> struct CudaMultTraits;

template<> struct CudaMultTraits<double> {
  using cuda_type = double;
  static auto convert(const double x) { return x; }
  static auto convert_back(const double x) { return x; }
  static void gemm(cublasHandle_t handle, const cublasOperation_t transa, const cublasOperation_t transb,
                   const int m, const int n, const int k, const double *alpha, const double *A, const int lda,
                   const double *B, const int ldb, const double *beta, double *C, const int ldc) {
    cublas_check(cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc), "cublasDgemm");
  }
  static constexpr auto adjoint_op = CUBLAS_OP_T;
};

template<> struct CudaMultTraits<std::complex<double>> {
  using cuda_type = cuDoubleComplex;
  static auto convert(const std::complex<double> x) { return make_cuDoubleComplex(x.real(), x.imag()); }
  static auto convert_back(const cuDoubleComplex x) { return std::complex<double>(cuCreal(x), cuCimag(x)); }
  static void gemm(cublasHandle_t handle, const cublasOperation_t transa, const cublasOperation_t transb,
                   const int m, const int n, const int k, const cuDoubleComplex *alpha, const cuDoubleComplex *A, const int lda,
                   const cuDoubleComplex *B, const int ldb, const cuDoubleComplex *beta, cuDoubleComplex *C, const int ldc) {
    cublas_check(cublasZgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc), "cublasZgemm");
  }
  static constexpr auto adjoint_op = CUBLAS_OP_C;
};

template<scalar S, Eigen_matrix EM>
auto cuda_mult_copy_matrix(const EM &M) {
  using CT = typename CudaMultTraits<S>::cuda_type;
  std::vector<CT> host(size1(M) * size2(M));
  for (const auto i : range0(size1(M)))
    for (const auto j : range0(size2(M))) host[i * size2(M) + j] = CudaMultTraits<S>::convert(M(i, j));
  return host;
}

template<scalar S, Eigen_matrix EM>
void cuda_mult_copy_back(EM &M, const std::vector<typename CudaMultTraits<S>::cuda_type> &host) {
  for (const auto i : range0(size1(M)))
    for (const auto j : range0(size2(M))) M(i, j) = CudaMultTraits<S>::convert_back(host[i * size2(M) + j]);
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product_CUDA(EM &M, const t_coef factor, const EM &A, const EM &B) {
  if (!finite_size(A) || !finite_size(B)) return;
  assert(size1(M) == size1(A) && size2(A) == size2(B) && size1(B) == size2(M));
  assert(my_isfinite(factor));
  using Traits = CudaMultTraits<S>;
  using CT = typename Traits::cuda_type;
  const auto m = static_cast<int>(size1(A));
  const auto n = static_cast<int>(size1(B));
  const auto k = static_cast<int>(size2(A));
  auto h_M = cuda_mult_copy_matrix<S>(M);
  auto h_A = cuda_mult_copy_matrix<S>(A);
  auto h_B = cuda_mult_copy_matrix<S>(B);
  CudaMultBuffer<CT> d_M(h_M.size()), d_A(h_A.size()), d_B(h_B.size());
  cuda_mult_check(cudaMemcpy(d_M.data(), h_M.data(), d_M.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D M");
  cuda_mult_check(cudaMemcpy(d_A.data(), h_A.data(), d_A.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D A");
  cuda_mult_check(cudaMemcpy(d_B.data(), h_B.data(), d_B.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D B");
  const auto alpha = Traits::convert(S(factor));
  const auto beta = Traits::convert(S(1.0));
  CublasHandle handle;
  Traits::gemm(handle, Traits::adjoint_op, CUBLAS_OP_N, n, m, k, &alpha, d_B.data(), k, d_A.data(), k, &beta, d_M.data(), n);
  cuda_mult_check(cudaMemcpy(h_M.data(), d_M.data(), d_M.bytes(), cudaMemcpyDeviceToHost), "cudaMemcpy D2H M");
  cuda_mult_copy_back<S>(M, h_M);
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform_CUDA(EM &M, const t_coef factor, const EM &A, const EM &O, const EM &B) {
  if (!finite_size(A) || !finite_size(B)) return;
  assert(size1(M) == size1(A) && size2(A) == size1(O) && size2(O) == size2(B) && size1(B) == size2(M));
  assert(my_isfinite(factor));
  using Traits = CudaMultTraits<S>;
  using CT = typename Traits::cuda_type;
  const auto m = static_cast<int>(size1(A));
  const auto a = static_cast<int>(size2(A));
  const auto b = static_cast<int>(size2(O));
  const auto n = static_cast<int>(size1(B));
  auto h_M = cuda_mult_copy_matrix<S>(M);
  auto h_A = cuda_mult_copy_matrix<S>(A);
  auto h_O = cuda_mult_copy_matrix<S>(O);
  auto h_B = cuda_mult_copy_matrix<S>(B);
  CudaMultBuffer<CT> d_M(h_M.size()), d_A(h_A.size()), d_O(h_O.size()), d_B(h_B.size()), d_tmp(static_cast<size_t>(b) * static_cast<size_t>(m));
  cuda_mult_check(cudaMemcpy(d_M.data(), h_M.data(), d_M.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D M");
  cuda_mult_check(cudaMemcpy(d_A.data(), h_A.data(), d_A.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D A");
  cuda_mult_check(cudaMemcpy(d_O.data(), h_O.data(), d_O.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D O");
  cuda_mult_check(cudaMemcpy(d_B.data(), h_B.data(), d_B.bytes(), cudaMemcpyHostToDevice), "cudaMemcpy H2D B");
  const auto one = Traits::convert(S(1.0));
  const auto zero = Traits::convert(S(0.0));
  const auto alpha = Traits::convert(S(factor));
  CublasHandle handle;
  Traits::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, b, m, a, &one, d_O.data(), b, d_A.data(), a, &zero, d_tmp.data(), b);
  Traits::gemm(handle, Traits::adjoint_op, CUBLAS_OP_N, n, m, b, &alpha, d_B.data(), b, d_tmp.data(), b, &one, d_M.data(), n);
  cuda_mult_check(cudaMemcpy(h_M.data(), d_M.data(), d_M.bytes(), cudaMemcpyDeviceToHost), "cudaMemcpy D2H M");
  cuda_mult_copy_back<S>(M, h_M);
}

#else
inline void validate_cuda_multiplication_request(const std::string &mult) {
  if (cuda_mult_requested(mult)) throw std::runtime_error("mult=cuda requested, but CUDA support was not enabled at build time");
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product_CUDA(EM &, const t_coef, const EM &, const EM &) {
  throw std::runtime_error("product_CUDA requested, but CUDA support was not enabled at build time");
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform_CUDA(EM &, const t_coef, const EM &, const EM &, const EM &) {
  throw std::runtime_error("transform_CUDA requested, but CUDA support was not enabled at build time");
}
#endif

#endif
