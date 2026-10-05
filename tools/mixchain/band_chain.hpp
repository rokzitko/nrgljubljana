// Channel-mixing discretization for NRG
// ** Reduction of the star to a band matrix by plane rotations

#ifndef _mixchain_band_chain_hpp_
#define _mixchain_band_chain_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "types.hpp"

namespace NRG::MixChain {

// THE BAND REDUCTION
//
// G. S. Ammar and W. B. Gragg, SIAM J. Matrix Anal. Appl. 12, 426 (1991), Algorithm 1: the generalization to p
// channels of the rotations of Rutishauser that scalar_star_to_chain() of the nrg library uses for one. The star, M
// levels with energies E_k and an M x p start matrix A whose row k holds the couplings of level k, is the bordered
// matrix
//
//   [ 0   A^dag ]
//   [ A   D     ],   D = diag(E_k),
//
// and a unitary transformation of the levels, which leaves the first p coordinates alone, brings it to a band matrix
// of half-bandwidth p: the chain, with p x p blocks,
//
//   A = Q_0 R,   E_n = Q_n^dag D Q_n,   T_n = Q_{n+1}^dag D Q_n,
//
// in the gauge where R and every T_n are upper triangular with a real nonnegative diagonal. The transformation is
// built from plane rotations, so no orthogonality is lost and double precision is enough.
//
// The levels are added one at a time. A new level is a row with p border entries and its energy on the diagonal; each
// border entry in turn is annihilated by a rotation against the outermost band element of an earlier row, the pivot,
// and the fill-in this creates is the next entry to annihilate, one column further on. The row has at most 2p+1
// nonzero entries at any time.
//
// Only the leading 'blocks' blocks are kept. Once the band is full, a further level is chased through it in the same
// way and then dropped: the rotations that would bring it into band form act on rows beyond the kept ones, so the
// kept blocks are those of the full reduction.
//
// The result does not depend on the order of the levels mathematically, but its rounding error does; see
// chain.hpp for the order to hand them over in.
//
// Nothing is decided about ranks here. A level without coupling passes through untouched and takes a row of the band,
// so the caller removes those; a pivot that is zero when the entry to annihilate is not is an exact exchange of the
// two rows. If the star has no more levels than the band has rows, the rows that came last have not served as pivots,
// and their outermost elements are not real nonnegative.

template<typename S> struct BandChain {
  Matrix<S> R;              // p x p: A = Q_0 R
  std::vector<Matrix<S>> E; // E_0..E_{blocks-1}, Hermitian
  std::vector<Matrix<S>> T; // T_0..T_{blocks-2}
  std::size_t rows{};       // the rows of the band that were filled, at most p*blocks; fewer if the star is short
};

template<typename S>
BandChain<S> band_star_to_chain(const std::vector<double> &energies, const Matrix<S> &start, const std::size_t blocks) {
  using std::abs;
  using Eigen::numext::conj;
  using Eigen::numext::real;

  const auto levels = start.rows();
  const auto p      = start.cols();
  if (p < 1) throw std::invalid_argument("band: the start matrix has no columns.");
  if (blocks < 1) throw std::invalid_argument("band: the number of blocks must be positive.");
  if (static_cast<Eigen::Index>(energies.size()) != levels)
    throw std::invalid_argument("band: one energy is needed for every row of the start matrix.");
  double largest = 0.0;
  for (Eigen::Index k = 0; k < levels; k++) {
    if (!std::isfinite(energies[static_cast<std::size_t>(k)])) throw std::invalid_argument("band: energies must be finite.");
    for (Eigen::Index i = 0; i < p; i++) {
      const auto size = abs(start(k, i));
      if (!std::isfinite(size)) throw std::invalid_argument("band: couplings must be finite.");
      largest = std::max(largest, size);
    }
  }
  // A common power of two in the couplings only scales R. Taking it out keeps the norms below away from the ends of
  // the range of double.
  int exponent = 0;
  std::frexp(largest, &exponent);
  const auto scale = std::scalbn(1.0, 1 - exponent);

  const auto kept = p * static_cast<Eigen::Index>(blocks); // rows of the band
  const auto n    = p + kept;                              // with the p coordinates of the border
  Matrix<S> band  = Matrix<S>::Zero(n, n);                 // Hermitian, both triangles kept
  Eigen::Index filled = 0;
  Vector<S> row(n); // the level being added: its elements in the columns of the band
  for (Eigen::Index level = 0; level < levels; level++) {
    row.setZero();
    for (Eigen::Index i = 0; i < p; i++) row(i) = scale * start(level, i);
    auto diagonal = energies[static_cast<std::size_t>(level)];

    const auto last = p + filled - 1; // the last row of the band so far
    for (Eigen::Index j = p; j <= last; j++) {
      const auto l = j - p;
      const S g    = row(l);
      if (g == S(0)) continue;
      const S f = band(j, l);
      // The unitary [[a, b], [-conj(b), conj(a)]] takes (f, g) to (rho, 0).
      const auto rho = std::hypot(abs(f), abs(g));
      S a = conj(f) / rho, b = conj(g) / rho;
      if (rho < std::numeric_limits<double>::min()) {
        // A subnormal norm has few bits; take the direction from operands scaled to order one.
        const auto size = std::max(abs(f), abs(g));
        const S fs = f / size, gs = g / size;
        const auto norm = std::hypot(abs(fs), abs(gs));
        a = conj(fs) / norm;
        b = conj(gs) / norm;
      }
      for (Eigen::Index c = l + 1; c <= std::min(j + p, last); c++) {
        if (c == j) continue;
        const S upper = band(j, c), lower = row(c);
        band(j, c) = a * upper + b * lower;
        band(c, j) = conj(band(j, c));
        row(c)     = -conj(b) * upper + conj(a) * lower;
      }
      const auto u = real(band(j, j)), v = diagonal;
      const S w        = row(j);
      const auto cross = 2.0 * real(a * conj(b) * conj(w));
      const auto aa = abs(a) * abs(a), bb = abs(b) * abs(b);
      band(j, j) = make_scalar<S>(aa * u + bb * v + cross, 0);
      diagonal   = bb * u + aa * v - cross;
      row(j)     = conj(a) * conj(b) * (v - u) + conj(a) * conj(a) * w - conj(b) * conj(b) * conj(w);
      row(l)     = S(0);
      band(j, l) = make_scalar<S>(rho, 0);
      band(l, j) = band(j, l);
      if (!std::isfinite(rho) || !std::isfinite(real(band(j, j))) || !std::isfinite(diagonal))
        throw std::runtime_error("band: numerical breakdown or nonrepresentable intermediate coefficient.");
    }

    if (filled == kept) continue; // the band is full: the level has been chased through it and is dropped
    const auto k = p + filled++;
    for (Eigen::Index c = 0; c < k; c++) {
      band(k, c) = row(c);
      band(c, k) = conj(row(c));
    }
    band(k, k) = make_scalar<S>(diagonal, 0);
  }

  BandChain<S> result;
  result.rows = static_cast<std::size_t>(filled);
  // By the real factor: dividing by a complex number squares its modulus on the way, which over- or underflows here.
  result.R = band.block(p, 0, p, p) * (1.0 / scale);
  for (std::size_t site = 0; site < blocks; site++) {
    const auto offset = p + p * static_cast<Eigen::Index>(site);
    result.E.push_back(band.block(offset, offset, p, p));
    if (site + 1 < blocks) result.T.push_back(band.block(offset + p, offset, p, p));
  }
  return result;
}

} // namespace NRG::MixChain

#endif
