#ifndef _io_hpp_
#define _io_hpp_

#include <string>
#include <fstream>
#include <iostream>
#include <complex>
#include "params.hpp"

#ifdef USE_UBLAS
#include "numerics_ublas.hpp"  // reim
#else
#include "numerics_Eigen.hpp"  // reim
#endif

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

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

} // namespace

#endif
