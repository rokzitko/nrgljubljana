// Channel-mixing discretization for NRG
// ** Writing of the Wilson chain

#ifndef _mixchain_chain_io_hpp_
#define _mixchain_chain_io_hpp_

#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include "../common/io.hpp"
#include "chain.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// THE CHAIN FILE
//
// A text file describing the Wilson chain
//
//   H = sum_ij ( V_ij d_i^dag f_{0j} + h.c. )
//       + sum_n sum_ij (E_n)_ij f_{ni}^dag f_{nj}
//       + sum_n sum_ij ( (T_n)_ij f_{n+1,i}^dag f_{nj} + h.c. ),
//
// with one row per matrix element. It is written by the 'l' mode and by the default mode.
//
// Lines beginning with '#' are comments. The second line is the header, as whitespace-separated key=value pairs:
//
//   channels     the dimension N of every block
//   Nmax         the last site; the chain has the sites 0..Nmax
//   z, Lambda    the discretization the star was built with
//   bandrescale  the band rescaling applied to Gamma; the coefficients are in the rescaled band, whose edge is 1
//   complex      1 if the coefficients are complex, 0 if real. It fixes the number of columns.
//   digits       the decimal digits of the arithmetic the recursion ran in
//
// The third line holds the diagnostics of the recursion. Columns of a data row:
//
//   block  V, E or T
//   n      the site: 0 for V, 0..Nmax for E, 0..Nmax-1 for T
//   i, j   the matrix indices, 1..channels, as in Gamma_ij
//   value  the element; for complex=1 a pair Re Im
//
// V_ij multiplies d_i^dag f_{0j}, (E_n)_ij multiplies f_{ni}^dag f_{nj}, and (T_n)_ij multiplies f_{n+1,i}^dag f_{nj}.
// In the polar gauge V and every T_n are Hermitian positive semidefinite and every E_n is Hermitian.
//
// The recursion runs in multiprecision because the late coefficients fall off as Lambda^(-n/2), but the result is
// written as double with 18 significant digits, as nrgchain writes xi.dat: the extra digits are needed to get the
// recursion right, not to use its result.
//
// V is in the normalization of the input, V^2 = Theta = pi sum_k V_k V_k^dag for the physical couplings V_k. A
// writer for a particular symmetry type of nrg applies the factor sqrt(1/pi) where it needs the physical amplitude.

inline constexpr auto chain_default_filename = "chain.dat";

// What the file records that the chain itself does not carry.
struct ChainFileHeader {
  double z{};
  double Lambda{};
  double bandrescale{1.0};
  unsigned digits{};
};

namespace detail {

// One element as double: its value, or its real and imaginary parts. Works for double, std::complex<double> and the
// wide types alike, so a chain can be written in whatever arithmetic it was computed in.
template<typename S> void write_element(std::ostream &out, const S &x) {
  out << " " << static_cast<double>(Eigen::numext::real(x));
  if constexpr (is_complex_v<S>) out << " " << static_cast<double>(Eigen::numext::imag(x));
}

template<typename S>
void write_block(std::ostream &out, const char *name, const unsigned int site, const Matrix<S> &block) {
  for (Eigen::Index i = 0; i < block.rows(); i++) {
    for (Eigen::Index j = 0; j < block.cols(); j++) {
      out << name << " " << site << " " << i + 1 << " " << j + 1;
      write_element(out, block(i, j));
      out << "\n";
    }
  }
}

} // namespace detail

template<typename S> void save_chain(const Chain<S> &chain, const ChainFileHeader &header, std::ostream &out) {
  const auto &d = chain.diagnostics;
  out << std::setprecision(18);
  out << "# mixchain Wilson chain" << std::endl;
  out << "# channels=" << chain.channels << " Nmax=" << chain.Nmax << " z=" << header.z << " Lambda=" << header.Lambda
      << " bandrescale=" << header.bandrescale << " complex=" << (is_complex_v<S> ? 1 : 0)
      << " digits=" << header.digits << std::endl;
  out << "# theta_condition=" << d.theta_condition << " max_antihermitian=" << d.max_antihermitian
      << " max_reorthogonalization=" << d.max_reorthogonalization
      << " min_residual_condition=" << d.min_residual_condition << std::endl;
  out << (is_complex_v<S> ? "# block n i j Re Im" : "# block n i j value") << std::endl;

  detail::write_block(out, "V", 0, chain.V);
  for (unsigned int n = 0; n < chain.E.size(); n++) detail::write_block(out, "E", n, chain.E[n]);
  for (unsigned int n = 0; n < chain.T.size(); n++) detail::write_block(out, "T", n, chain.T[n]);
  out.flush();
}

template<typename S>
void save_chain(const Chain<S> &chain, const ChainFileHeader &header, const std::string &filename) {
  std::ofstream F;
  NRG::Tools::open_output(F, filename, 18);
  save_chain(chain, header, F);
  F.close();
  if (!F) throw std::runtime_error("Error writing " + filename + ".");
}

} // namespace NRG::MixChain

#endif
