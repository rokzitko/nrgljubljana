// Channel-mixing discretization for NRG
// ** Writing of the Wilson chain

#ifndef _mixchain_chain_io_hpp_
#define _mixchain_chain_io_hpp_

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "../common/io.hpp"
#include "blocks.hpp"
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
//   Nmax         the last site; the chain has the sites 0..Nmax, and one hopping per site, T_0..T_Nmax, the last one
//                leading out of the chain, as the coefficient tables of nrg are indexed 0..Nmax
//   z, Lambda    the discretization the star was built with
//   bandrescale  the band rescaling that was applied to Gamma when it was read
//   complex      1 if the coefficients are complex, 0 if real. It fixes the number of columns.
//   digits       the decimal digits of the arithmetic the recursion ran in
//   gauge        'polar', where V and every T_n are Hermitian positive semidefinite, or 'nambu', where the hole
//                component of every second site is flipped so that each block has the Nambu structure
//
// When the star has more than one block, the header is followed by the blocks line of the star file, as in
// "# blocks= {1,3} {2}": each block was mapped onto a chain of its own, and every element between channels of
// different blocks is exactly zero. The next line holds the diagnostics of the recursion over the whole chain.
// Columns of a data row:
//
//   block  V, E or T
//   n      the site: 0 for V, 0..Nmax for E and for T
//   i, j   the matrix indices, 1..channels, as in Gamma_ij
//   value  the element; for complex=1 a pair Re Im
//
// V_ij multiplies d_i^dag f_{0j}, (E_n)_ij multiplies f_{ni}^dag f_{nj}, and (T_n)_ij multiplies f_{n+1,i}^dag f_{nj}.
// In the polar gauge V and every T_n are Hermitian positive semidefinite and every E_n is Hermitian.
//
// Units: the discretization runs in the band rescaled to the edge 1, but E_n and T_n are written multiplied by
// bandrescale, in the units of the input. That is what nrgchain writes into xi.dat and zeta.dat, which apply the
// same factor, and what nrg expects, since its SCALE(N) carries bandrescale as well. V needs no factor: rescaling
// omega and Gamma leaves the integral of Gamma, hence Theta and V, unchanged.
//
// The recursion runs in multiprecision because the late coefficients fall off as Lambda^(-n/2), but the result is
// written as double with 18 significant digits, as nrgchain writes xi.dat: the extra digits are needed to get the
// recursion right, not to use its result.
//
// V is in the normalization of the input: V^2 = Theta = int Gamma domega, whatever Gamma was given. With the
// convention of the dos file of adapt, where Gamma is pi times the spectral function, that is pi sum_k V_k V_k^dag
// for the physical couplings V_k, and a writer for a symmetry type of nrg applies the factor sqrt(1/pi) where it
// needs the physical amplitude.

inline constexpr auto chain_default_filename = "chain.dat";

// What the file records that the chain itself does not carry.
struct ChainFileHeader {
  double z{};
  double Lambda{};
  double bandrescale{1.0};
  unsigned digits{};
  double innermost_input{}; // as recorded by the star, in the rescaled band; 0 if unknown
};

namespace detail {

// One element as double: its value, or its real and imaginary parts. Works for double, std::complex<double> and the
// wide types alike, so a chain can be written in whatever arithmetic it was computed in.
template<typename S> void write_element(std::ostream &out, const S &x) {
  out << " " << static_cast<double>(Eigen::numext::real(x));
  if constexpr (is_complex_v<S>) out << " " << static_cast<double>(Eigen::numext::imag(x));
}

template<typename S>
void write_block(std::ostream &out, const char *name, const unsigned int site, const Matrix<S> &block,
                 const double factor = 1.0) {
  for (Eigen::Index i = 0; i < block.rows(); i++) {
    for (Eigen::Index j = 0; j < block.cols(); j++) {
      out << name << " " << site << " " << i + 1 << " " << j + 1;
      write_element(out, make_scalar<S>(static_cast<real_type<S>>(factor), 0) * block(i, j));
      out << "\n";
    }
  }
}

} // namespace detail

template<typename S> void save_chain(const Chain<S> &chain, const ChainFileHeader &header, std::ostream &out) {
  const auto &d        = chain.diagnostics;
  const auto continued = first_continued_site(chain, header.innermost_input);
  out << std::setprecision(18);
  out << "# mixchain Wilson chain" << std::endl;
  out << "# channels=" << chain.channels << " Nmax=" << chain.Nmax << " z=" << header.z << " Lambda=" << header.Lambda
      << " bandrescale=" << header.bandrescale << " complex=" << (is_complex_v<S> ? 1 : 0)
      << " digits=" << header.digits << " gauge=" << chain_gauge_name(chain.gauge) << std::endl;
  if (chain.blocks.size() > 1) out << "# blocks= " << blocks_name(chain.blocks) << std::endl;
  out << "# levels=" << d.levels << " coupled_levels=" << d.coupled_levels << " theta_rank=" << d.theta_rank
      << " min_rank=" << d.min_rank << " rank_drop_site="
      << (d.rank_drop_site ? std::to_string(*d.rank_drop_site) : std::string("none"))
      << " continued_from_site="
      << (continued ? std::to_string(*continued) : std::string("none"))
      << " theta_condition=" << d.theta_condition << " max_antihermitian=" << d.max_antihermitian
      << " max_reorthogonalization=" << d.max_reorthogonalization
      << " min_residual_condition=" << d.min_residual_condition << std::endl;
  out << (is_complex_v<S> ? "# block n i j Re Im" : "# block n i j value") << std::endl;

  // The recursion runs in the rescaled band; E_n and T_n are written in the units of the input. V is invariant under
  // the rescaling and is written as it is.
  detail::write_block(out, "V", 0, chain.V);
  for (unsigned int n = 0; n < chain.E.size(); n++) detail::write_block(out, "E", n, chain.E[n], header.bandrescale);
  for (unsigned int n = 0; n < chain.T.size(); n++) detail::write_block(out, "T", n, chain.T[n], header.bandrescale);
  out.flush();
}

// THE MATRIX FILES
//
// The same chain, one file per matrix element: V11.dat, V12.dat, ..., E11.dat, ..., T11.dat, ..., in the directory of
// chain.dat. V holds one row, E and T the sites 0..Nmax, one row each. The rows are plain
// numbers with no header, a single column for a real chain and the pair "Re Im" for a complex one, in the units of
// chain.dat: E and T carry bandrescale, V does not. Every element is written, including those that are exactly zero
// between blocks, so the set of N*N files is always complete.
//
// This is the form the coefficient readers of nrg take. Which element belongs to which coefficient set of a given
// symmetry type is up to the writer that stages them.
template<typename S>
void save_chain_matrix_files(const Chain<S> &chain, const ChainFileHeader &header,
                             const std::filesystem::path &directory) {
  const auto element = [](std::ostream &out, const S &x, const double factor) {
    const S scaled = make_scalar<S>(static_cast<real_type<S>>(factor), 0) * x;
    out << static_cast<double>(Eigen::numext::real(scaled));
    if constexpr (is_complex_v<S>) out << " " << static_cast<double>(Eigen::numext::imag(scaled));
    out << std::endl;
  };
  for (Eigen::Index i = 0; i < chain.V.rows(); i++)
    for (Eigen::Index j = 0; j < chain.V.cols(); j++) {
      const auto indices = std::to_string(i + 1) + std::to_string(j + 1);
      for (const auto &[name, blocks, factor] :
           {std::tuple{"V", std::vector<Matrix<S>>{chain.V}, 1.0}, std::tuple{"E", chain.E, header.bandrescale},
            std::tuple{"T", chain.T, header.bandrescale}}) {
        const auto filename = (directory / (name + indices + ".dat")).string();
        std::ofstream F;
        NRG::Tools::open_output(F, filename, 18);
        for (const auto &block : blocks) element(F, block(i, j), factor);
        F.close();
        if (!F) throw std::runtime_error("Error writing " + filename + ".");
      }
    }
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
