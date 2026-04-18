#ifndef _basicio_hpp_
#define _basicio_hpp_

#include <algorithm> // copy
#include <string>
#include <complex>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip> // set_precision
#include <iterator> // ostream_iterator
#include <stdexcept>
#include <cctype>
#include <strings.h> // strcasecmp

#include <fmt/format.h>

#include <boost/lexical_cast.hpp>

#include "portabil.hpp"
#include "traits.hpp"

#define HIGHPREC(val) std::setprecision(std::numeric_limits<double>::max_digits10) << (val)

namespace NRG {

inline std::string strip_leading_whitespace(const std::string &in) {
  const auto first = std::find_if_not(in.begin(), in.end(), [](const unsigned char ch) { return std::isspace(ch); });
  return std::string(first, in.end());
}

inline std::string strip_trailing_whitespace(const std::string &in) {
  auto s(in);
  auto it = s.rbegin();
  while (it != s.rend() && std::isspace(static_cast<unsigned char>(*it))) {
    s.erase(--it.base());
    it = s.rbegin();
  }
  return s;
}

inline std::string strip_whitespace(const std::string &in) {
  return strip_trailing_whitespace(strip_leading_whitespace(in));
}

template <class T>
inline T from_string(const std::string &str) {
  T result;
  try {
    result = boost::lexical_cast<T>(str);
  } catch (boost::bad_lexical_cast &) { throw std::runtime_error(fmt::format("Lexical cast [{}] failed.", str)); }
  return result;
}

template <>
inline bool from_string(const std::string &str) {
  const auto trimmed = strip_whitespace(str);
  if (trimmed == "1") return true;
  if (trimmed == "0") return false;
  if (strcasecmp(trimmed.c_str(), "true") == 0 || strcasecmp(trimmed.c_str(), "yes") == 0) return true;
  if (strcasecmp(trimmed.c_str(), "false") == 0 || strcasecmp(trimmed.c_str(), "no") == 0) return false;
  throw std::runtime_error(fmt::format("Lexical cast [{}] failed.", str));
}

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
   std::ifstream F(fn, std::ios::binary | std::ios::in);
   return bool(F);
}

inline size_t count_words_in_string(const std::string &s) { // size_t because it should be unsigned
  std::stringstream stream(s);
  return std::distance(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>());
}

inline bool is_blank_or_comment_line(const std::string &s) {
  const auto first = std::find_if_not(s.begin(), s.end(), [](const unsigned char ch) { return std::isspace(ch); });
  return first == s.end() || *first == '#';
}

// Read one object of type S from an object F which has an extractor operator. Use as read_one<S>(F).
template<typename S, typename T>
inline auto read_one(T &F) {
  S value{};
  F >> value;
  if constexpr (requires { F.fail(); })
    if (F.fail()) throw std::runtime_error("read_one() error. Input is corrupted.");
  return value;
}
  
// Determine the matrix dimensions from a stream of rows of white-space-separated tabulated values
inline auto get_dims(std::istream &F) {
  size_t dim1 = 0; // number of rows
  size_t dim2 = 0; // number of columns
  std::string s;
  while (std::getline(F, s)) {
    if (is_blank_or_comment_line(s)) continue;
    const auto n = count_words_in_string(s);
    if (dim2 > 0 && dim2 != n) throw std::runtime_error("All matrix rows must be equally long");
    dim2 = n;
    dim1++;
  }
  return std::make_pair(dim1, dim2);
}

// Read 'size' values of type T into a std::vector<T>.
template <scalar T> auto read_std_vector(std::istream &F, const size_t size) {
  std::vector<T> vec(size);
  for (auto j = 0; j < size; j++)
    vec[j] = read_one<T>(F);
  if (F.fail()) throw std::runtime_error("read_std_vector() error. Input file is corrupted.");
  return vec;
}

// Read values of type T into a std::vector<T>. First value to be read, 'nr', is either vector dimension or
// the value of maximum index.
template <scalar T> auto read_std_vector(std::istream &F, const bool nr_is_max_index = false) {
  const auto nr = read_one<size_t>(F);
  const auto len = nr_is_max_index ? nr+1 : nr;
  return read_std_vector<T>(F, len);
}

} // namespace

#endif
