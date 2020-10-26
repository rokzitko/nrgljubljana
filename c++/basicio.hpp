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
using namespace fmt::literals;

#include <boost/lexical_cast.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric;

namespace NRG {

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
inline std::string to_string(const T val) { return boost::lexical_cast<std::string>(val); }

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


} // namespace

#endif
