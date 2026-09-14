// Channel-mixing discretization for NRG
// ** Matrix and scalar types

#ifndef _mixchain_types_hpp_
#define _mixchain_types_hpp_

#include <complex>
#include <string>

#include <Eigen/Dense>

namespace NRG::MixChain {

// The same shapes as NRG::EigenMatrix and NRG::EigenVector in c++/traits.hpp, but without the 'scalar' concept: the
// block Lanczos stage instantiates these with a multiprecision scalar, which is neither a floating point type nor a
// std::complex of one. A column vector must stay column-major.
template<typename S> using Matrix = Eigen::Matrix<S, -1, -1, Eigen::RowMajor>;
template<typename S> using Vector = Eigen::Matrix<S, -1, 1>;

template<typename T> struct scalar_traits {
  using real_type                   = T;
  static constexpr bool is_complex = false;
};

template<typename T> struct scalar_traits<std::complex<T>> {
  using real_type                   = T;
  static constexpr bool is_complex = true;
};

template<typename S> using real_type = typename scalar_traits<S>::real_type;
template<typename S> inline constexpr bool is_complex_v = scalar_traits<S>::is_complex;

// Assemble a scalar from its real and imaginary parts. For a real S the imaginary part is dropped, and the caller is
// responsible for having checked that it vanishes.
template<typename S> inline S make_scalar(const real_type<S> re, [[maybe_unused]] const real_type<S> im) {
  if constexpr (is_complex_v<S>)
    return S(re, im);
  else
    return re;
}

// Positive and negative frequencies are discretized separately, as in adapt and nrgchain.
enum class Sign { POS, NEG };

inline auto sign_name(const Sign sign) { return sign == Sign::POS ? std::string("POS") : std::string("NEG"); }
inline auto sign_value(const Sign sign) { return sign == Sign::POS ? 1.0 : -1.0; }

} // namespace NRG::MixChain

#endif
