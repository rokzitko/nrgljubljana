// Channel-mixing discretization for NRG
// ** Plane rotations: from the star to the Wilson chain, in double precision

#ifndef _mixchain_chain_rkpw_hpp_
#define _mixchain_chain_rkpw_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include <star-to-chain.hpp>

#include "band_chain.hpp"
#include "blocks.hpp"
#include "chain.hpp"
#include "star.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// The chain of chain.hpp, built by adding the levels of the star one at a time and restoring the form of the chain
// with plane rotations, instead of by the Lanczos recursion. The rotations are orthogonal transformations of the
// bath, so nothing is lost to cancellation and double precision is enough.
//
// A block of one channel is the scalar problem, and goes through scalar_star_to_chain() of the nrg library, the
// Rutishauser-Kahan-Pal-Walker rotations that nrgchain uses. A block of several channels is reduced to a band matrix
// by band_star_to_chain(), and the chain is then read off that band in a second stage, which is where the ranks are
// decided and the polar gauge is fixed.

namespace detail {

// The levels of a block in the order the rotations take them: interval by interval from the band edge inwards, the
// two frequency branches alternating. Levels without an interval index keep the order they came in.
template<typename S0> void sort_for_insertion(std::vector<const StarLevel<S0> *> &levels) {
  std::stable_sort(levels.begin(), levels.end(), [](const StarLevel<S0> *a, const StarLevel<S0> *b) {
    const auto key = [](const StarLevel<S0> *level) {
      return std::tuple(level->m, level->sign != Sign::POS, level->branch);
    };
    return key(a) < key(b);
  });
}

// The chain of one channel: xi[n] couples the sites n and n+1.
struct ScalarBlockChain {
  double V{};
  std::vector<double> zeta, xi; // Nmax+1 of each, zero beyond the support of the star
  ChainDiagnostics diagnostics;
};

// 'levels' are those of the block, 'channel' the one channel it has.
//
// The levels are handed over interval by interval from the band edge inwards, with the two frequency branches
// alternating, whatever their order in the star: the result does not depend on the order mathematically, but its
// rounding error does, and this is the order that keeps it small at the end of the chain. Levels that carry no
// interval index, as in a star that was not produced by the star stage, keep the order they came in.
template<typename S0>
ScalarBlockChain scalar_block_chain(std::vector<const StarLevel<S0> *> levels, const Eigen::Index channel,
                                    const ChainOptions &options) {
  using std::abs;
  const auto sites  = static_cast<std::size_t>(options.Nmax) + 1;
  const auto needed = sites + 1;
  if (levels.size() < needed)
    throw std::invalid_argument("The star has " + std::to_string(levels.size()) + " levels, but a chain of "
                                + std::to_string(sites) + " sites with one channel needs at least "
                                + std::to_string(needed)
                                + ", one block more than the sites, for the hopping out of the last site"
                                + ". Increase mMAX or decrease Nmax.");

  ScalarBlockChain result;
  result.zeta.assign(sites, 0.0);
  result.xi.assign(sites, 0.0);
  auto &diagnostics  = result.diagnostics;
  diagnostics.levels = static_cast<int>(levels.size());

  std::erase_if(levels, [channel](const StarLevel<S0> *level) { return abs(level->coupling(channel)) == 0.0; });
  diagnostics.coupled_levels = static_cast<int>(levels.size());
  sort_for_insertion(levels);

  std::vector<NRG::StarPoint> points;
  points.reserve(levels.size());
  std::set<double> energies; // levels of the same energy are one pole of the hybridization
  double largest = 0.0;
  for (const auto *level : levels) {
    points.push_back({level->energy, abs(level->coupling(channel))});
    energies.insert(level->energy);
    largest = std::max(largest, points.back().amplitude);
  }
  const auto support = energies.size();

  diagnostics.hopping_ranks.assign(sites, 0);
  if (support == 0) return result; // nothing couples: Theta is zero, and so is the chain
  diagnostics.theta_rank      = 1;
  diagnostics.theta_condition = 1.0;

  // V^2 = Theta = sum_k |v_k|^2, relative to the largest term so that the squares neither overflow nor underflow.
  double sum = 0.0;
  for (const auto &point : points) sum += (point.amplitude / largest) * (point.amplitude / largest);
  result.V = largest * std::sqrt(sum);

  // A star with fewer poles than the chain has sites ends early: the hopping out of its last site is exactly zero,
  // and so is everything beyond.
  const auto count = std::min(sites, support);
  const auto chain = NRG::scalar_star_to_chain(points, count);
  std::copy(chain.zeta.begin(), chain.zeta.end(), result.zeta.begin());
  std::copy(chain.xi.begin(), chain.xi.end(), result.xi.begin());

  diagnostics.min_rank = 1;
  for (std::size_t n = 0; n < sites; n++) {
    const auto rank              = result.xi[n] > 0.0 ? 1 : 0;
    diagnostics.hopping_ranks[n] = rank;
    diagnostics.min_rank         = std::min(diagnostics.min_rank, rank);
    if (rank == 0 && !diagnostics.rank_drop_site) diagnostics.rank_drop_site = static_cast<unsigned int>(n);
  }
  return result;
}

// The chain of one block with its own channels 0..p-1.
template<typename S0> struct BlockChain {
  Matrix<S0> V;
  std::vector<Matrix<S0>> E, T; // Nmax+1 of each
  ChainDiagnostics diagnostics;
};

// THE SECOND STAGE
//
// band_star_to_chain() decides no rank. Where Theta is singular, or the Krylov space of the star runs out in some
// direction, a pivot of the band is rounding and the rows after it come in no particular order. The band is still a
// unitary transformation of the star, exact to rounding, and site s of the chain lies within its first s+1 blocks; so
// the chain can be read off the kept band, site by site.
//
// M is the kept band as a dense Hermitian matrix. At every step there are the rows not yet given to a site, their
// coupling C to the previous site (at the start the block R, the coupling to the impurity), and an isometry Phi that
// says where the orbitals of the previous site sit among the channels (at the start the identity):
//
//   C = Q [R_1; 0]      Householder QR over the rows that couple to the site, and M -> Q^dag M Q on those rows,
//   R_1 = U S X^dag     the rank r counts the singular values that are not zero by rank_tolerance,
//   M -> U^dag M U      on the first rows, of which the first r are the new site,
//   Phi' = Phi X_r      hopping Phi' S Phi'^dag, on-site block Phi' M[site, site] Phi'^dag.
//
// The hopping is Hermitian positive semidefinite, which is the polar gauge of chain.hpp, and for a hopping of lower
// rank this is the pseudo-inverse convention of the Lanczos recursion: the chain is zero along the directions that
// are lost, from there on. The rows that are left over stay among those not yet given to a site, where a later site
// may still reach them through M. With full rank throughout, Q is trivial and a step is the SVD of one block.
//
// A singular value counts as zero when its square is below rank_tolerance times the square of the largest, as an
// eigenvalue of the Gram matrix does in the Lanczos recursion; and all of them do when the largest is rounding on the
// scale of the bath Hamiltonian applied to the site.
template<typename S0>
BlockChain<S0> matrix_block_chain(std::vector<const StarLevel<S0> *> levels, const Block &block,
                                  const ChainOptions &options) {
  using Dense = Eigen::Matrix<S0, Eigen::Dynamic, Eigen::Dynamic>;
  const auto p      = static_cast<Eigen::Index>(block.size());
  const auto sites  = static_cast<std::size_t>(options.Nmax) + 1;
  const auto needed = static_cast<std::size_t>(p) * (sites + 1);
  if (levels.size() < needed)
    throw std::invalid_argument("The star has " + std::to_string(levels.size()) + " levels, but a chain of "
                                + std::to_string(sites) + " sites with " + std::to_string(p)
                                + " channels needs at least " + std::to_string(needed)
                                + ", one block more than the sites, for the hopping out of the last site"
                                + ". Increase mMAX or decrease Nmax.");

  BlockChain<S0> result;
  result.V = Matrix<S0>::Zero(p, p);
  result.E.assign(sites, Matrix<S0>::Zero(p, p));
  result.T.assign(sites, Matrix<S0>::Zero(p, p));
  auto &diagnostics  = result.diagnostics;
  diagnostics.levels = static_cast<int>(levels.size());
  diagnostics.hopping_ranks.assign(sites, 0);

  // A level without coupling would pass through the rotations untouched and take a row of the band.
  std::erase_if(levels, [&block](const StarLevel<S0> *level) {
    for (const auto channel : block)
      if (level->coupling(channel) != S0(0)) return false;
    return true;
  });
  diagnostics.coupled_levels = static_cast<int>(levels.size());
  sort_for_insertion(levels);

  std::vector<double> energies;
  energies.reserve(levels.size());
  Matrix<S0> start(static_cast<Eigen::Index>(levels.size()), p);
  for (std::size_t k = 0; k < levels.size(); k++) {
    energies.push_back(levels[k]->energy);
    for (Eigen::Index i = 0; i < p; i++)
      start(static_cast<Eigen::Index>(k), i) = Eigen::numext::conj(levels[k]->coupling(block[static_cast<std::size_t>(i)]));
  }
  const auto band = band_star_to_chain(energies, start, sites + 1);
  const auto rows = static_cast<Eigen::Index>(band.rows);

  Dense M = Dense::Zero(rows, rows);
  Dense C = Dense::Zero(rows, p);
  for (Eigen::Index i = 0; i < rows; i++) {
    for (Eigen::Index j = 0; j < rows; j++) {
      const auto si = i / p, sj = j / p;
      if (si == sj) M(i, j) = band.E[static_cast<std::size_t>(si)](i % p, j % p);
      if (si == sj + 1) M(i, j) = band.T[static_cast<std::size_t>(sj)](i % p, j % p);
      if (sj == si + 1) M(i, j) = Eigen::numext::conj(band.T[static_cast<std::size_t>(si)](j % p, i % p));
    }
    if (i < p) C.row(i) = band.R.row(i);
  }

  const auto epsilon = std::numeric_limits<double>::epsilon();
  Dense Phi          = Dense::Identity(p, p);
  Eigen::Index first = 0; // the first row not yet given to a site
  double scale       = 0.0; // the norm of the bath Hamiltonian applied to the previous site; nothing for the impurity
  for (std::size_t step = 0; step <= sites; step++) { // step 0 gives V and site 0, step n+1 gives T_n and site n+1
    const auto rest     = rows - first;
    const auto previous = C.cols();
    Eigen::Index window = 0; // the rows that couple to the previous site are the first 'window' of the rest
    for (Eigen::Index i = 0; i < C.rows(); i++)
      if ((C.row(i).array() != S0(0)).any()) window = i + 1;

    Eigen::Index rank = 0;
    Dense Phi_new     = Dense::Zero(p, 0);
    Dense hopping     = Dense::Zero(p, p);
    if (window > 0 && previous > 0) {
      Eigen::HouseholderQR<Dense> qr(Dense(C.topRows(window)));
      const Dense Q = qr.householderQ();
      M.block(first, first, window, rest).applyOnTheLeft(Q.adjoint());
      M.block(first, first, rest, window).applyOnTheRight(Q);
      const auto k   = std::min(window, previous);
      const Dense R1 = qr.matrixQR().topRows(k).template triangularView<Eigen::Upper>();
      Eigen::JacobiSVD<Dense> svd(R1, Eigen::ComputeFullU | Eigen::ComputeFullV);
      const auto &sigma = svd.singularValues();
      if (sigma(0) > 0.0 && sigma(0) > 1000.0 * epsilon * scale)
        for (Eigen::Index i = 0; i < k; i++)
          if ((sigma(i) / sigma(0)) * (sigma(i) / sigma(0)) >= options.rank_tolerance) rank++;
      const Dense U = svd.matrixU();
      M.block(first, first, k, rest).applyOnTheLeft(U.adjoint());
      M.block(first, first, rest, k).applyOnTheRight(U);

      Phi_new = Phi * svd.matrixV().leftCols(rank);
      hopping = Phi_new * sigma.head(rank).template cast<S0>().asDiagonal() * Phi_new.adjoint();
      hopping = (make_scalar<S0>(0.5, 0) * (hopping + hopping.adjoint())).eval();
      if (rank > 0) {
        const auto condition = (sigma(rank - 1) / sigma(0)) * (sigma(rank - 1) / sigma(0));
        if (step == 0)
          diagnostics.theta_condition = condition;
        else
          diagnostics.min_residual_condition = std::min(diagnostics.min_residual_condition, condition);
      }
    }

    if (step == 0) {
      result.V               = hopping;
      diagnostics.theta_rank = static_cast<int>(rank);
      diagnostics.min_rank   = static_cast<int>(rank);
    } else {
      const auto n = step - 1;
      result.T[n]  = hopping;
      diagnostics.hopping_ranks[n] = static_cast<int>(rank);
      diagnostics.min_rank         = std::min(diagnostics.min_rank, static_cast<int>(rank));
      if (static_cast<int>(rank) < diagnostics.theta_rank && !diagnostics.rank_drop_site)
        diagnostics.rank_drop_site = static_cast<unsigned int>(n);
    }

    if (step < sites && rank > 0) {
      const Dense onsite = Phi_new * M.block(first, first, rank, rank) * Phi_new.adjoint();
      result.E[step]     = make_scalar<S0>(0.5, 0) * (onsite + onsite.adjoint());
      scale              = M.block(0, first, rows, rank).norm();
      C                  = M.block(first + rank, first, rest - rank, rank);
    } else {
      C = Dense::Zero(rest - rank, 0);
    }
    first += rank;
    Phi = Phi_new;
  }
  return result;
}

} // namespace detail

// The chain of a star in the arithmetic of the star: each block mapped onto a chain of its own, with exact zeros
// between channels of different blocks, as build_chain() of chain_lanczos.hpp does.
template<typename S0> Chain<S0> build_chain_rkpw(const Star<S0> &star, const ChainOptions &options) {
  const auto channels = static_cast<Eigen::Index>(star.channels);
  if (channels < 1) throw std::invalid_argument("The star has no channels.");
  if (options.Nmax < 1) throw std::invalid_argument("Nmax must be greater than 0.");

  Chain<S0> chain;
  chain.channels = star.channels;
  chain.Nmax     = options.Nmax;
  chain.blocks   = star.blocks;
  if (chain.blocks.empty()) { // a single block of all channels
    chain.blocks.emplace_back(static_cast<std::size_t>(channels));
    std::iota(chain.blocks.front().begin(), chain.blocks.front().end(), 0);
  }
  chain.V = Matrix<S0>::Zero(channels, channels);
  chain.E.assign(options.Nmax + 1, Matrix<S0>::Zero(channels, channels));
  chain.T.assign(options.Nmax + 1, Matrix<S0>::Zero(channels, channels));

  int offset = 0; // the first branch label of the current block
  for (const auto &block : chain.blocks) {
    const auto size = static_cast<int>(block.size());

    // As in build_chain(): the levels of a block are those of its branches, with or without coupling.
    std::vector<const StarLevel<S0> *> levels;
    for (const auto &level : star.levels)
      if (chain.blocks.size() == 1 || (level.branch >= offset && level.branch < offset + size))
        levels.push_back(&level);

    if (size == 1) {
      const auto channel = static_cast<Eigen::Index>(block.front());
      const auto piece   = detail::scalar_block_chain(std::move(levels), channel, options);
      chain.V(channel, channel) = make_scalar<S0>(piece.V, 0);
      for (unsigned int n = 0; n <= options.Nmax; n++) {
        chain.E[n](channel, channel) = make_scalar<S0>(piece.zeta[n], 0);
        chain.T[n](channel, channel) = make_scalar<S0>(piece.xi[n], 0);
      }
      chain.block_diagnostics.push_back(piece.diagnostics);
    } else {
      const auto piece = detail::matrix_block_chain(std::move(levels), block, options);
      const auto place = [&block, size](Matrix<S0> &target, const Matrix<S0> &source) {
        for (int i = 0; i < size; i++)
          for (int j = 0; j < size; j++)
            target(block[static_cast<std::size_t>(i)], block[static_cast<std::size_t>(j)]) = source(i, j);
      };
      place(chain.V, piece.V);
      for (unsigned int n = 0; n <= options.Nmax; n++) place(chain.E[n], piece.E[n]);
      for (unsigned int n = 0; n <= options.Nmax; n++) place(chain.T[n], piece.T[n]);
      chain.block_diagnostics.push_back(piece.diagnostics);
    }
    offset += size;
  }
  chain.diagnostics = detail::merge_diagnostics(chain.block_diagnostics, options.Nmax + 1);
  if (options.gauge == ChainGauge::nambu) detail::apply_nambu_gauge(chain, options.nambu_tolerance);
  return chain;
}

// THE SENSITIVITY OF THE CHAIN TO THE STAR
//
// The star is stored in double precision, so each of its numbers is known to one unit in the last place at best. For
// most stars that moves the chain by rounding. Where the mesh accumulates at a finite energy, as at a gap edge, the
// levels close to it differ in digits that double precision does not hold, and the late sites of the chain, which are
// built from those differences, move by many orders of magnitude more. No method can determine them better from such
// a star, in whatever arithmetic it runs: multiprecision then gives the exact chain of numbers that are not exact.
//
// This measures it: the chain is built again from the star with every energy and coupling moved by one unit in the
// last place, up or down at random, for a few fixed choices, and compared site by site, the hoppings relative to
// their largest element and the on-site blocks on the scale of their site. Always with the rotations, which cost
// next to nothing, and in the polar gauge.
struct StarSensitivity {
  double largest{};
  unsigned int largest_site{};
  std::optional<unsigned int> from_site; // the first site above the tolerance
};

template<typename S0> StarSensitivity star_sensitivity(const Star<S0> &star, ChainOptions options) {
  options.gauge = ChainGauge::polar;
  const auto size      = [](const Matrix<S0> &m) { return m.size() ? m.cwiseAbs().maxCoeff() : 0.0; };
  const auto reference = build_chain_rkpw(star, options);
  std::vector<double> moved(options.Nmax + 1, 0.0);
  for (const unsigned int seed : {1U, 2U, 3U, 4U}) {
    std::mt19937 generator(seed);
    const auto nudge = [&generator](const double x) { return std::nextafter(x, generator() % 2 ? 2.0 * x : 0.0); };
    auto other       = star;
    for (auto &level : other.levels) {
      level.energy = nudge(level.energy);
      for (Eigen::Index i = 0; i < level.coupling.size(); i++) {
        if constexpr (is_complex_v<S0>)
          level.coupling(i) = S0(nudge(level.coupling(i).real()), nudge(level.coupling(i).imag()));
        else
          level.coupling(i) = nudge(level.coupling(i));
      }
    }
    const auto chain = build_chain_rkpw(other, options);
    for (unsigned int n = 0; n <= options.Nmax; n++) {
      const auto hopping = size(reference.T[n]);
      const auto local   = std::max({size(reference.E[n]), hopping, n > 0 ? size(reference.T[n - 1]) : 0.0});
      if (hopping > 0.0) moved[n] = std::max(moved[n], size(chain.T[n] - reference.T[n]) / hopping);
      if (local > 0.0) moved[n] = std::max(moved[n], size(chain.E[n] - reference.E[n]) / local);
    }
  }
  StarSensitivity result;
  for (unsigned int n = 0; n <= options.Nmax; n++) {
    if (moved[n] > result.largest) {
      result.largest      = moved[n];
      result.largest_site = n;
    }
    if (moved[n] > options.sensitivity_tolerance && !result.from_site) result.from_site = n;
  }
  return result;
}

} // namespace NRG::MixChain

#endif
