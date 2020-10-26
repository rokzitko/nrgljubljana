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

namespace NRG {

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
