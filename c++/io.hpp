#ifndef _io_hpp_
#define _io_hpp_

#include <string>
#include <fstream>
#include <iostream>
#include <complex>
#include "params.hpp"
#include "numerics.hpp"// reim

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

namespace fmt {
template <typename S, typename... Args,
  FMT_ENABLE_IF(detail::is_string<S>::value)>
    auto color_print(bool enable_color, const text_style& ts, const S& format_str, const Args&... args) {
      return enable_color ? print(stdout, ts, format_str, args...)
                          : print(stdout, format_str, args...);
    }
} // namespace fmt

namespace NRG {

using namespace fmt::literals;

template <typename T>
inline std::string formatted_output(const T x, const Params &P) {
  return fmt::format("{x:>{width}}", "x"_a=x, "width"_a=P.width_custom);
}

inline std::string formatted_output(const double x, const Params &P) {
  return fmt::format("{x:>{width}.{prec}}", "x"_a=x, "prec"_a=P.prec_custom, "width"_a=P.width_custom);
}

inline std::string formatted_output(const std::complex<double> z, const Params &P) { // XXX
  return fmt::format("{x:>{width}.{prec}}", "x"_a=z.real(), "prec"_a=P.prec_custom, "width"_a=P.width_custom);
}

inline void outputxy(std::ostream &F, const double x, const std::complex<double> z, const bool imagpart, const double clip_tol_imag = 1e-10) {
  const auto [r, i] = reim(z);
  F << x << " " << r;
  if (imagpart) F << " " << (abs(i)>abs(r)*clip_tol_imag ? i : 0);
  F << std::endl;
}

template <matrix M> std::ostream &operator<<(std::ostream &os, const M &m) {
  for (auto r1 = 0; r1 < size1(m); r1++) {
    for (auto r2 = 0; r2 < size2(m); r2++)
      os << m(r1, r2) << ' ';
    os << std::endl;
  }
  return os;
}

// Read dim1 x dim2 matrix from stream. Use function next_value to extract consecutive values. 'gen' generates the matrix.
template<typename GEN, typename FNC>
auto read_matrix_data(GEN && generate_matrix, FNC && next_value, const size_t dim1, const size_t dim2, const bool check_is_finite = true) {
  auto M = generate_matrix(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto x = next_value();
      if (check_is_finite && !std::isfinite(x)) throw std::runtime_error("Non-finite number detected.");
      M(i, j) = x;
    }
  }
  return M;
}

// Read a matrix from stream (text)
template<typename GEN>
auto read_matrix_text(GEN && generate_matrix, const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, false);
  const auto [dim1, dim2] = get_dims(F);
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  F.clear();
  F.seekg (0, std::ios::beg);
  return read_matrix_data(generate_matrix, [&F]() { return read_one<double>(F); }, dim1, dim2);
}

// Read a matrix from stream (binary). Format: two unit32_t for matrix size, followed by
// dim1 x dim2 double values.
template<typename GEN>
auto read_matrix_bin(GEN && generate_matrix, const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, true);
  uint32_t dim1, dim2;
  F.read((char *)&dim1, sizeof(uint32_t));
  F.read((char *)&dim2, sizeof(uint32_t));
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  return read_matrix_data(generate_matrix, [&F]() { double x; F.read((char *)&x, sizeof(double)); return x; }, dim1, dim2);
}

#ifdef INCL_UBLAS
inline auto read_matrix_ublas(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? read_matrix_bin(generate_ublas<double>, filename, verbose) 
               : read_matrix_text(generate_ublas<double>, filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}
#endif

#ifdef INCL_EIGEN
inline auto read_matrix_Eigen(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? read_matrix_bin(generate_Eigen<double>, filename, verbose) 
               : read_matrix_text(generate_Eigen<double>, filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}
#endif

inline auto read_matrix(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? read_matrix_bin(generate_matrix<double>, filename, verbose) 
               : read_matrix_text(generate_matrix<double>, filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}

template<matrix M>
inline void save_matrix(const std::string &filename, const M &m, const bool verbose = false,
                        const double chop_tol = 1e-14, const int output_prec = 18)
{
  if (verbose) std::cout << "Saving result to " << filename << std::endl;
  auto F = safe_open(filename);
  F << std::setprecision(output_prec);
  for (auto i = 0; i < size1(m); i++) {
    for (auto j = 0; j < size2(m); j++) {
      const auto val = m(i, j);
      F << (std::abs(val) > chop_tol ? val : 0.0) << (j != size2(m) - 1 ? " " : "");
    }
    F << std::endl;
  }
  F.close();
}

} // namespace

#endif
