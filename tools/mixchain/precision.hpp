// Channel-mixing discretization for NRG
// ** Multiprecision arithmetic for the block Lanczos step

#ifndef _mixchain_precision_hpp_
#define _mixchain_precision_hpp_

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "types.hpp"

namespace NRG::MixChain {

// Multiple precision for Lanczos tridiagonalization

// The scalar types come from Boost.Multiprecision, which tools/common/piecewise_polynomial.hpp already uses, rather
// than from GMP: the block Lanczos needs complex arithmetic and a Hermitian square root, and Boost's complex wrapper
// around the GMP backend does not compile (it wants eval_signbit, which that backend does not provide), while the
// alternative complex backend is MPC, which is not a dependency of this project. The cpp_bin_float family is header
// only, works with Eigen's decompositions, and needs nothing new.
//
// Its precision is a compile-time constant, whereas preccpp is a runtime parameter, so the chain is instantiated for
// a few precisions and the smallest one that covers the request is used. Two requests that resolve to the same rung
// give bit-identical results; resolve_precision() reports the rung so that this is visible rather than silent.
//
// Eigen and Boost both have expression templates, and they interact badly, so the numbers are et_off. cpp_complex is
// et_off already.
template<unsigned Digits>
using WideReal = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<Digits>,
                                               boost::multiprecision::et_off>;
template<unsigned Digits> using WideComplex = boost::multiprecision::cpp_complex<Digits>;

// Decimal digits. Each rung instantiates the whole recursion twice, once real and once complex, so they are few and
// far apart. The last one covers the default preccpp of 2000 bits, which is 603 digits.
inline constexpr unsigned precision_ladder[] = {50, 200, 800};

// The decimal digits needed to carry 'bits' binary digits: rounded up, since a partial digit does not cover them.
inline constexpr auto digits_for_bits(const unsigned int bits) {
  return static_cast<unsigned>(std::ceil(static_cast<double>(bits) * 0.30102999566398119521));
}

// The binary digits a rung of 'digits' decimal digits guarantees: rounded down, for the same reason in reverse. The
// two are then inverse in the sense that resolve_precision(bits_for_digits(rung)) is that rung.
inline constexpr auto bits_for_digits(const unsigned digits) {
  return static_cast<unsigned int>(std::floor(static_cast<double>(digits) / 0.30102999566398119521));
}

// The rung to use for a request of 'bits' binary digits, in decimal digits.
inline auto resolve_precision(const unsigned int bits) {
  if (bits <= 10) throw std::invalid_argument("preccpp must be greater than 10.");
  const auto requested = digits_for_bits(bits);
  for (const auto rung : precision_ladder)
    if (rung >= requested) return rung;
  throw std::invalid_argument("preccpp=" + std::to_string(bits) + " asks for " + std::to_string(requested)
                              + " decimal digits, more than the " + std::to_string(precision_ladder[std::size(precision_ladder) - 1])
                              + " digits this tool is compiled for.");
}

// Call f.template operator()<Scalar>() with the multiprecision type of the rung that covers 'bits'. The caller
// writes a templated lambda, and the chain code is instantiated once per rung from there.
template<typename F> auto with_precision(const unsigned int bits, const bool complex_data, F &&f) {
  switch (resolve_precision(bits)) {
    case precision_ladder[0]:
      return complex_data ? f.template operator()<WideComplex<precision_ladder[0]>>()
                          : f.template operator()<WideReal<precision_ladder[0]>>();
    case precision_ladder[1]:
      return complex_data ? f.template operator()<WideComplex<precision_ladder[1]>>()
                          : f.template operator()<WideReal<precision_ladder[1]>>();
    default:
      return complex_data ? f.template operator()<WideComplex<precision_ladder[2]>>()
                          : f.template operator()<WideReal<precision_ladder[2]>>();
  }
}

// The same, for a caller that knows at compile time whether its data are complex: the scalar is real or complex as
// Like is, and only that kind is instantiated. This is what building a chain from a star needs, since a complex star
// cannot be widened into a real type, and with_precision() above would instantiate that combination too.
template<typename Like, typename F> auto with_precision_like(const unsigned int bits, F &&f) {
  const auto call = [&f]<unsigned Digits>() {
    if constexpr (is_complex_v<Like>)
      return f.template operator()<WideComplex<Digits>>();
    else
      return f.template operator()<WideReal<Digits>>();
  };
  switch (resolve_precision(bits)) {
    case precision_ladder[0]: return call.template operator()<precision_ladder[0]>();
    case precision_ladder[1]: return call.template operator()<precision_ladder[1]>();
    default: return call.template operator()<precision_ladder[2]>();
  }
}

} // namespace NRG::MixChain

#endif
