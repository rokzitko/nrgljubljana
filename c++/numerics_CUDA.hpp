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
#include <memory>
#include <unordered_map>
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
   CudaMultBuffer() = default;
   explicit CudaMultBuffer(const size_t size) : size_(size) {
     if (size_ != 0) cuda_mult_check(cudaMalloc(reinterpret_cast<void **>(&data_), size_ * sizeof(T)), "cudaMalloc");
   }
   ~CudaMultBuffer() { if (data_ != nullptr) cudaFree(data_); }
   CudaMultBuffer(const CudaMultBuffer &) = delete;
   CudaMultBuffer &operator=(const CudaMultBuffer &) = delete;
   CudaMultBuffer(CudaMultBuffer &&other) noexcept : data_(other.data_), size_(other.size_) {
     other.data_ = nullptr;
     other.size_ = 0;
   }
   CudaMultBuffer &operator=(CudaMultBuffer &&other) noexcept {
     if (this != &other) {
       if (data_ != nullptr) cudaFree(data_);
       data_ = other.data_;
       size_ = other.size_;
       other.data_ = nullptr;
       other.size_ = 0;
     }
     return *this;
   }
   T *data() { return data_; }
   const T *data() const { return data_; }
   void ensure(const size_t size) {
     if (size <= size_) return;
     if (data_ != nullptr) cuda_mult_check(cudaFree(data_), "cudaFree");
     size_ = size;
     cuda_mult_check(cudaMalloc(reinterpret_cast<void **>(&data_), size_ * sizeof(T)), "cudaMalloc");
   }
   size_t bytes() const { return bytes(size_); }
   static constexpr size_t bytes(const size_t size) { return size * sizeof(T); }

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

template<scalar S> class CudaRecalcScope {
 public:
  using Traits = CudaMultTraits<S>;
  using CT = typename Traits::cuda_type;

  CudaRecalcScope() : previous_(current_) { current_ = this; }
  ~CudaRecalcScope() { current_ = previous_; }
  CudaRecalcScope(const CudaRecalcScope &) = delete;
  CudaRecalcScope &operator=(const CudaRecalcScope &) = delete;

  template<Eigen_matrix EM>
  void register_matrix(const EM &M) {
    if (!finite_size(M)) return;
    const auto key = make_key(M);
    if (cache_.contains(key)) return;
    auto host = cuda_mult_copy_matrix<S>(M);
    auto cached = std::make_unique<CachedMatrix>(host.size());
    cached->rows = size1(M);
    cached->cols = size2(M);
    cuda_mult_check(cudaMemcpy(cached->buffer.data(), host.data(), CudaMultBuffer<CT>::bytes(host.size()), cudaMemcpyHostToDevice), "cudaMemcpy H2D cached matrix");
    cache_.emplace(key, std::move(cached));
  }

  template<Eigen_matrix EM>
  const CT *cached_data(const EM &M) const {
    if (!finite_size(M)) return nullptr;
    const auto it = cache_.find(make_key(M));
    return it != cache_.end() ? it->second->buffer.data() : nullptr;
  }

  CublasHandle &handle() { return handle_; }
  CudaMultBuffer<CT> &output_buffer(const size_t size) { output_.ensure(size); return output_; }
  CudaMultBuffer<CT> &temp_buffer(const size_t size) { temp_.ensure(size); return temp_; }
  CudaMultBuffer<CT> &operand_a_buffer(const size_t size) { operand_a_.ensure(size); return operand_a_; }
  CudaMultBuffer<CT> &operand_b_buffer(const size_t size) { operand_b_.ensure(size); return operand_b_; }
  CudaMultBuffer<CT> &operand_o_buffer(const size_t size) { operand_o_.ensure(size); return operand_o_; }

  template<Eigen_matrix EM>
  const CT *upload_uncached(const EM &M, CudaMultBuffer<CT> &buffer, const char *what) {
    auto host = cuda_mult_copy_matrix<S>(M);
    buffer.ensure(host.size());
    cuda_mult_check(cudaMemcpy(buffer.data(), host.data(), CudaMultBuffer<CT>::bytes(host.size()), cudaMemcpyHostToDevice), what);
    return buffer.data();
  }

  static CudaRecalcScope *current() { return current_; }

 private:
  struct MatrixKey {
    const void *data{};
    size_t rows{};
    size_t cols{};
    bool operator==(const MatrixKey &) const = default;
  };

  struct MatrixKeyHash {
    size_t operator()(const MatrixKey &key) const noexcept {
      const auto h1 = std::hash<const void *>{}(key.data);
      const auto h2 = std::hash<size_t>{}(key.rows);
      const auto h3 = std::hash<size_t>{}(key.cols);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };

  struct CachedMatrix {
    explicit CachedMatrix(const size_t size) : buffer(size) {}
    size_t rows{};
    size_t cols{};
    CudaMultBuffer<CT> buffer;
  };

  template<Eigen_matrix EM>
  static MatrixKey make_key(const EM &M) {
    return {M.data(), size1(M), size2(M)};
  }

  inline static thread_local CudaRecalcScope *current_ = nullptr;

  CudaRecalcScope *previous_{};
  CublasHandle handle_;
  std::unordered_map<MatrixKey, std::unique_ptr<CachedMatrix>, MatrixKeyHash> cache_;
  CudaMultBuffer<CT> output_;
  CudaMultBuffer<CT> temp_;
  CudaMultBuffer<CT> operand_a_;
  CudaMultBuffer<CT> operand_b_;
  CudaMultBuffer<CT> operand_o_;
};

template<scalar S, Eigen_matrix EM> class CudaRecalcAccumulator {
 public:
  using Traits = CudaMultTraits<S>;
  using CT = typename Traits::cuda_type;

  explicit CudaRecalcAccumulator(EM &M) : M_(M), ctx_(CudaRecalcScope<S>::current()) {
    if (ctx_ == nullptr) return;
    host_ = cuda_mult_copy_matrix<S>(M_);
    device_ = &ctx_->output_buffer(host_.size());
    cuda_mult_check(cudaMemcpy(device_->data(), host_.data(), CudaMultBuffer<CT>::bytes(host_.size()), cudaMemcpyHostToDevice), "cudaMemcpy H2D accumulator");
  }

  [[nodiscard]] bool active() const noexcept { return ctx_ != nullptr; }
  CudaRecalcScope<S> &context() { return *ctx_; }
  CT *data() { return device_->data(); }
  EM &matrix() { return M_; }

  void copy_back() {
    if (!active()) return;
    cuda_mult_check(cudaMemcpy(host_.data(), device_->data(), CudaMultBuffer<CT>::bytes(host_.size()), cudaMemcpyDeviceToHost), "cudaMemcpy D2H accumulator");
    cuda_mult_copy_back<S>(M_, host_);
  }

 private:
  EM &M_;
  CudaRecalcScope<S> *ctx_{};
  CudaMultBuffer<CT> *device_{};
  std::vector<CT> host_;
};

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product_CUDA_accumulate(CudaRecalcAccumulator<S, EM> &accumulator, const t_coef factor, const EM &A, const EM &B) {
  if (!finite_size(A) || !finite_size(B)) return;
  [[maybe_unused]] auto &M = accumulator.matrix();
  assert(size1(M) == size1(A) && size2(A) == size2(B) && size1(B) == size2(M));
  assert(my_isfinite(factor));
  using Traits = CudaMultTraits<S>;
  const auto m = static_cast<int>(size1(A));
  const auto n = static_cast<int>(size1(B));
  const auto k = static_cast<int>(size2(A));
  auto &ctx = accumulator.context();
  const auto *d_A = ctx.cached_data(A);
  if (d_A == nullptr) d_A = ctx.upload_uncached(A, ctx.operand_a_buffer(size1(A) * size2(A)), "cudaMemcpy H2D A");
  const auto *d_B = ctx.cached_data(B);
  if (d_B == nullptr) d_B = ctx.upload_uncached(B, ctx.operand_b_buffer(size1(B) * size2(B)), "cudaMemcpy H2D B");
  const auto alpha = Traits::convert(S(factor));
  const auto beta = Traits::convert(S(1.0));
  Traits::gemm(ctx.handle(), Traits::adjoint_op, CUBLAS_OP_N, n, m, k, &alpha, d_B, k, d_A, k, &beta, accumulator.data(), n);
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform_CUDA_accumulate(CudaRecalcAccumulator<S, EM> &accumulator, const t_coef factor, const EM &A, const EM &O, const EM &B) {
  if (!finite_size(A) || !finite_size(B)) return;
  [[maybe_unused]] auto &M = accumulator.matrix();
  assert(size1(M) == size1(A) && size2(A) == size1(O) && size2(O) == size2(B) && size1(B) == size2(M));
  assert(my_isfinite(factor));
  using Traits = CudaMultTraits<S>;
  const auto m = static_cast<int>(size1(A));
  const auto a = static_cast<int>(size2(A));
  const auto b = static_cast<int>(size2(O));
  const auto n = static_cast<int>(size1(B));
  auto &ctx = accumulator.context();
  const auto *d_A = ctx.cached_data(A);
  if (d_A == nullptr) d_A = ctx.upload_uncached(A, ctx.operand_a_buffer(size1(A) * size2(A)), "cudaMemcpy H2D A");
  const auto *d_B = ctx.cached_data(B);
  if (d_B == nullptr) d_B = ctx.upload_uncached(B, ctx.operand_b_buffer(size1(B) * size2(B)), "cudaMemcpy H2D B");
  const auto *d_O = ctx.upload_uncached(O, ctx.operand_o_buffer(size1(O) * size2(O)), "cudaMemcpy H2D O");
  auto &d_tmp = ctx.temp_buffer(static_cast<size_t>(b) * static_cast<size_t>(m));
  const auto one = Traits::convert(S(1.0));
  const auto zero = Traits::convert(S(0.0));
  const auto alpha = Traits::convert(S(factor));
  Traits::gemm(ctx.handle(), CUBLAS_OP_N, CUBLAS_OP_N, b, m, a, &one, d_O, b, d_A, a, &zero, d_tmp.data(), b);
  Traits::gemm(ctx.handle(), Traits::adjoint_op, CUBLAS_OP_N, n, m, b, &alpha, d_B, b, d_tmp.data(), b, &one, accumulator.data(), n);
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
  if (auto *ctx = CudaRecalcScope<S>::current(); ctx != nullptr) {
    auto &d_M = ctx->output_buffer(h_M.size());
    cuda_mult_check(cudaMemcpy(d_M.data(), h_M.data(), CudaMultBuffer<CT>::bytes(h_M.size()), cudaMemcpyHostToDevice), "cudaMemcpy H2D M");
    const auto *d_A = ctx->cached_data(A);
    if (d_A == nullptr) d_A = ctx->upload_uncached(A, ctx->operand_a_buffer(size1(A) * size2(A)), "cudaMemcpy H2D A");
    const auto *d_B = ctx->cached_data(B);
    if (d_B == nullptr) d_B = ctx->upload_uncached(B, ctx->operand_b_buffer(size1(B) * size2(B)), "cudaMemcpy H2D B");
    const auto alpha = Traits::convert(S(factor));
    const auto beta = Traits::convert(S(1.0));
    Traits::gemm(ctx->handle(), Traits::adjoint_op, CUBLAS_OP_N, n, m, k, &alpha, d_B, k, d_A, k, &beta, d_M.data(), n);
    cuda_mult_check(cudaMemcpy(h_M.data(), d_M.data(), CudaMultBuffer<CT>::bytes(h_M.size()), cudaMemcpyDeviceToHost), "cudaMemcpy D2H M");
    cuda_mult_copy_back<S>(M, h_M);
    return;
  }
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
  if (auto *ctx = CudaRecalcScope<S>::current(); ctx != nullptr) {
    auto &d_M = ctx->output_buffer(h_M.size());
    cuda_mult_check(cudaMemcpy(d_M.data(), h_M.data(), CudaMultBuffer<CT>::bytes(h_M.size()), cudaMemcpyHostToDevice), "cudaMemcpy H2D M");
    const auto *d_A = ctx->cached_data(A);
    if (d_A == nullptr) d_A = ctx->upload_uncached(A, ctx->operand_a_buffer(size1(A) * size2(A)), "cudaMemcpy H2D A");
    const auto *d_B = ctx->cached_data(B);
    if (d_B == nullptr) d_B = ctx->upload_uncached(B, ctx->operand_b_buffer(size1(B) * size2(B)), "cudaMemcpy H2D B");
    const auto *d_O = ctx->upload_uncached(O, ctx->operand_o_buffer(size1(O) * size2(O)), "cudaMemcpy H2D O");
    auto &d_tmp = ctx->temp_buffer(static_cast<size_t>(b) * static_cast<size_t>(m));
    const auto one = Traits::convert(S(1.0));
    const auto zero = Traits::convert(S(0.0));
    const auto alpha = Traits::convert(S(factor));
    Traits::gemm(ctx->handle(), CUBLAS_OP_N, CUBLAS_OP_N, b, m, a, &one, d_O, b, d_A, a, &zero, d_tmp.data(), b);
    Traits::gemm(ctx->handle(), Traits::adjoint_op, CUBLAS_OP_N, n, m, b, &alpha, d_B, b, d_tmp.data(), b, &one, d_M.data(), n);
    cuda_mult_check(cudaMemcpy(h_M.data(), d_M.data(), CudaMultBuffer<CT>::bytes(h_M.size()), cudaMemcpyDeviceToHost), "cudaMemcpy D2H M");
    cuda_mult_copy_back<S>(M, h_M);
    return;
  }
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

template<scalar S> class CudaRecalcScope {
 public:
  CudaRecalcScope() = default;
  template<Eigen_matrix EM> void register_matrix(const EM &) {}
};

template<scalar S, Eigen_matrix EM> class CudaRecalcAccumulator {
 public:
  explicit CudaRecalcAccumulator(EM &) {}
  [[nodiscard]] bool active() const noexcept { return false; }
  void copy_back() {}
};

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product_CUDA_accumulate(CudaRecalcAccumulator<S, EM> &, const t_coef, const EM &, const EM &) {
  throw std::runtime_error("product_CUDA_accumulate requested, but CUDA support was not enabled at build time");
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform_CUDA_accumulate(CudaRecalcAccumulator<S, EM> &, const t_coef, const EM &, const EM &, const EM &) {
  throw std::runtime_error("transform_CUDA_accumulate requested, but CUDA support was not enabled at build time");
}
#endif

#endif
