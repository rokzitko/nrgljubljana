#ifndef _basicio_hpp_
#define _basicio_hpp_

#include <algorithm> // copy
#include <string>
#include <complex>
#include <set>
#include <fstream>
#include <iostream>
#include <iomanip> // set_precision
#include <iterator> // ostream_iterator
#include "portabil.hpp"

#define HIGHPREC(val) std::setprecision(std::numeric_limits<double>::max_digits10) << (val)

#include <boost/lexical_cast.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include <boost/lexical_cast.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <Eigen/Dense>

namespace NRG {

using namespace fmt::literals;
using namespace boost::numeric;

template <class T>
inline T from_string(const std::string &str) {
  T result;
  try {
    result = boost::lexical_cast<T>(str);
  } catch (boost::bad_lexical_cast &) { throw std::runtime_error(fmt::format("Lexical cast [{}] failed.", str)); }
  return result;
}

template <>
inline bool from_string(const std::string &str) { return (strcasecmp(str.c_str(), "true") == 0 ? true : false); }

// for T=int, std::to_string is used
template <class T>
inline std::string to_string(const T &val) { return boost::lexical_cast<std::string>(val); }

inline std::string to_string(const std::complex<double> &z) {
  std::ostringstream s;
  s << z;
  return s.str();
}

template <typename T>
std::ostream & operator<<(std::ostream &os, const std::set<T> &x) {
  std::copy(x.cbegin(), x.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) {
  return os << p.first << ' ' << p.second;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
  for (const auto &x : vec) os << x << " ";
  return os;
}
  
template <typename T> std::ostream &operator<<(std::ostream &os, const ublas::vector<T> &vec){
  for (const auto &x : vec) os << x << " ";
  return os;
}

template <typename T> std::ostream &operator<<(std::ostream &os, const ublas::matrix<T> &m) {
  for (auto r1 = 0; r1 < m.size1(); r1++) {
    for (auto r2 = 0; r2 < m.size2(); r2++)
      os << m(r1, r2) << ' ';
    os << std::endl;
  }
  return os;
}

// Returns a string with a floating value in fixed (non-exponential) format with N digits of precision after the
// decimal point.
inline std::string prec(const double x, const int N)
{
  std::ostringstream s;
  s << std::fixed << std::setprecision(N) << x;
  return s.str();
}
inline std::string prec3(const double x) { return prec(x, 3); }

template <typename T>
inline bool negligible_imag_part(const std::complex<T> &z, const double output_imag_eps = 1e-13) {
  return abs(z.imag()) < abs(z.real()) * output_imag_eps;
}

inline std::ofstream safe_open(const std::string &filename, const bool binary = false) {
  my_assert(filename != "");
  std::ios::openmode flag = std::ios::out;
  if (binary) flag |= std::ios::binary;
  std::ofstream F(filename, flag);
  if (!F) throw std::runtime_error(fmt::format("Can't open {} for writing", filename));
  return F;
}

inline std::ifstream safe_open_for_reading(const std::string &filename, const bool binary = false) {
  my_assert(filename != "");
  std::ios::openmode flag = std::ios::in;
  if (binary) flag |= std::ios::binary;
  std::ifstream F(filename, flag);
  if (!F) throw std::runtime_error(fmt::format("Can't open {} for reading", filename));
  return F;
}

inline bool file_exists(const std::string &fn)
{
   std::ofstream F(fn, std::ios::binary | std::ios::out);
   return bool(F);
}

inline auto count_words_in_string(const std::string &s) {
  std::stringstream stream(s);
  return std::distance(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>());
}

// Read one object of type S from an object F which has an extractor operator. Use as read_one<S>(F).
template<typename S, typename T>
inline auto read_one(T &F) {
  S value;
  F >> value;
  return value;
}
  
// Determine the matrix dimensions from a stream of rows of white-space-separated tabulated values
inline auto get_dims(std::istream &F) {
  auto dim1 = 0; // number of rows
  auto dim2 = 0; // number of columns
  while (F.good()) {
    std::string s;
    std::getline(F, s);
    if (!F.fail()) {
      auto n = count_words_in_string(s);
      if (dim2 > 0 && dim2 != n) throw std::runtime_error("All matrix rows must be equally long");
      dim2 = n;
      dim1++;
    }
  }
  return std::make_pair(dim1, dim2);
}

// Read dim1 x dim2 matrix from stream. Use function next_value to extract consecutive values.
template<typename FNC>
auto read_matrix_data(FNC next_value, const size_t dim1, const size_t dim2, const bool check_is_finite = true) {
  ublas::matrix<double> M(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto x = next_value();
      if (check_is_finite && !std::isfinite(x)) throw std::runtime_error("Non-finite number detected.");
      M(i, j) = x;
    }
  }
  // if (F.fail()) throw std::runtime_error("read_matrix_text() failed. Input corrupted?");
  return M;
}

// Read dim1 x dim2 matrix from stream. Use function next_value to extract consecutive values.
template<typename FNC>
auto _eigen_read_matrix_data(FNC next_value, const size_t dim1, const size_t dim2, const bool check_is_finite = true) {
  Eigen::MatrixXd M(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto x = next_value();
      if (check_is_finite && !finite(x)) throw std::runtime_error("Non-finite number detected.");      
      M(i, j) = x;
    }
  }
  // if (F.fail()) throw std::runtime_error("read_matrix_text() failed. Input corrupted?");
  return M;
}


// Read a matrix from stream (text)
inline auto read_matrix_text(const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, false);
  const auto [dim1, dim2] = get_dims(F);
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  F.clear();
  F.seekg (0, std::ios::beg);
  return read_matrix_data([&F]() { return read_one<double>(F); }, dim1, dim2);
}

// Read a matrix from stream (text)
inline auto _eigen_read_matrix_text(const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, false);
  const auto [dim1, dim2] = get_dims(F);
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  F.clear();
  F.seekg (0, std::ios::beg);
  return _eigen_read_matrix_data([&F]() { return read_one<double>(F); }, dim1, dim2);
}

// Read a matrix from stream (binary). Format: two unit32_t for matrix size, followed by
// dim1 x dim2 double values;
inline auto read_matrix_bin(const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, true);
  uint32_t dim1, dim2;
  F.read((char *)&dim1, sizeof(uint32_t));
  F.read((char *)&dim2, sizeof(uint32_t));
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  return read_matrix_data([&F]() { double x; F.read((char *)&x, sizeof(double)); return x; }, dim1, dim2);
}

// Read a matrix from stream (binary). Format: two unit32_t for matrix size, followed by
// dim1 x dim2 double values;
inline auto _eigen_read_matrix_bin(const std::string &filename, const bool verbose = false) {
  auto F = safe_open_for_reading(filename, true);
  uint32_t dim1, dim2;
  F.read((char *)&dim1, sizeof(uint32_t));
  F.read((char *)&dim2, sizeof(uint32_t));
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  return _eigen_read_matrix_data([&F]() { double x; F.read((char *)&x, sizeof(double)); return x; }, dim1, dim2);
}

inline auto read_matrix(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? read_matrix_bin(filename, verbose) : read_matrix_text(filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}

inline auto _eigen_read_matrix(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? _eigen_read_matrix_bin(filename, verbose) : _eigen_read_matrix_text(filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}

inline void save_matrix(const std::string &filename, const ublas::matrix<double> &M, const bool verbose = false,
                        const double chop_tol = 1e-14, const int output_prec = 18)
{
  if (verbose) std::cout << "Saving result to " << filename << std::endl;
  auto F = safe_open(filename);
  F << std::setprecision(output_prec);
  const auto dim1 = M.size1();
  const auto dim2 = M.size2();
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto val = M(i, j);
      F << (std::abs(val) > chop_tol ? val : 0.0) << (j != dim2 - 1 ? " " : "");
    }
    F << std::endl;
  }
  F.close();
}

inline void _eigen_save_matrix(const std::string &filename, const Eigen::MatrixXd &M, const bool verbose = false,
                        const double chop_tol = 1e-14, const int output_prec = 18)
{
  if (verbose) std::cout << "Saving result to " << filename << std::endl;
  auto F = safe_open(filename);
  F << std::setprecision(output_prec) << M;
  F.close();
}

} // namespace

#endif
