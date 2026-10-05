// Channel-mixing discretization for NRG
// ** The Wilson chain: from the star to the chain by plane rotations

#ifndef _mixchain_chain_hpp_
#define _mixchain_chain_hpp_

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
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <star-to-chain.hpp>

#include "band_chain.hpp"
#include "blocks.hpp"
#include "star.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// THE CHAIN
//
// The star
//
//   H = sum_k E_k c_k^dag c_k + sum_k sum_i ( v_{k,i} d_i^dag c_k + h.c. )
//
// is mapped to the Wilson chain
//
//   H = sum_ij ( V_ij d_i^dag f_{0j} + h.c. )
//       + sum_n sum_ij (E_n)_ij f_{ni}^dag f_{nj}
//       + sum_n sum_ij ( (T_n)_ij f_{n+1,i}^dag f_{nj} + h.c. ),
//
// with N x N blocks: E_n Hermitian, and in the polar gauge V and every T_n Hermitian positive semidefinite.
//
// In the single-particle space of the bath the levels |k> are the eigenstates, H|k> = E_k|k>. Impurity orbital i
// couples to b_i = sum_k v_{k,i} c_k, so the state it reaches is b_i^dag|0> = sum_k conj(v_{k,i}) |k>: the starting
// block is the M x N matrix A with A[k,i] = conj(v_{k,i}), whose Gram matrix A^dag A = sum_k v_k v_k^dag is Theta.
// Without the conjugation it would be the transpose of Theta, which for a complex Gamma is a different model.

// The gauge the chain is written in. The star fixes each site only up to a unitary rotation of its N orbitals, and the
// blocks of V, E_n and T_n transform together, so every gauge describes the same bath.
//
//   polar: V and every T_n Hermitian positive semidefinite, the matrix analogue of choosing xi_n > 0. Nothing is
//          assumed about what the channels mean, and every element is written, so a consumer that reads the whole
//          matrix (a pol2x2 template of nrg, or one coefficient set per channel) can use this as it is.
//   nambu: for blocks of two channels read as (particle, hole). The consumer of a superconducting chain stores only
//          xi = T(1,1), sckappa = T(1,2), zeta = E(1,1), scdelta = E(1,2) and reconstructs the rest from the Nambu
//          structure, so the chain must be in the gauge where that structure holds: E(2,2) = -E(1,1) and
//          T(2,2) = -conj(T(1,1)). The polar gauge is not: it absorbs the sign of the hole component into the orbitals
//          of the site, which turns a constant gap into one alternating along the chain. Flipping the hole component of
//          every second site, U_n = diag(1, (-1)^n), puts it back.
enum class ChainGauge { polar, nambu };

inline auto chain_gauge_from_string(const std::string &value) {
  if (value == "polar") return ChainGauge::polar;
  if (value == "nambu") return ChainGauge::nambu;
  throw std::invalid_argument("chain_gauge must be either 'polar' or 'nambu'.");
}

inline auto chain_gauge_name(const ChainGauge gauge) {
  return gauge == ChainGauge::polar ? std::string("polar") : std::string("nambu");
}

struct ChainOptions {
  // The chain has the sites 0..Nmax, and one hopping per site, T_0..T_Nmax: the last leads out of the chain and is
  // there because the coefficient tables of nrg are indexed 0..Nmax, as nrgchain writes xi.dat and zeta.dat.
  unsigned int Nmax{0};
  // A direction of the coupling of a site to the next, or of the impurity to the first, counts as lost when the square
  // of its singular value is below this fraction of the square of the largest: when the eigenvalue of the Gram matrix
  // of the couplings is below this fraction of the largest one.
  double rank_tolerance{1e-20};
  ChainGauge gauge{ChainGauge::polar};
  // How far a block may depart from the Nambu structure before the nambu gauge refuses the chain, relative to the
  // largest element of that block.
  double nambu_tolerance{1e-8};
  // The chain counts as no longer determined by the star from the first site that moves by more than this when every
  // number of the star is changed by one unit in the last place; see star_sensitivity().
  double sensitivity_tolerance{1e-10};
};

struct ChainDiagnostics {
  int theta_rank{};                 // the number of combinations of the impurity orbitals that couple to the bath
  double theta_condition{};         // the smallest nonzero eigenvalue of Theta over its largest; 0 if Theta is zero
  // The smallest ratio of the squares of the smallest and the largest nonzero singular value of a hopping along the
  // chain. It says how close a direction came to being counted as zero.
  double min_residual_condition{1.0};
  int min_rank{};                            // the smallest rank of a hopping T_n, and never above theta_rank
  std::optional<unsigned int> rank_drop_site; // the first n at which the rank of T_n is below theta_rank
  std::vector<int> hopping_ranks;             // the rank of every T_n
  // The levels of the star, and those with a nonzero coupling. Only the latter enter the Krylov space, so a block of
  // size s spans at most coupled_levels/s full sites.
  int levels{};
  int coupled_levels{};
  // In the nambu gauge, the largest departure from the Nambu structure of a block, relative to its largest element.
  double max_nambu_deviation{};
  // How far the chain moves when every number of the star is changed by one unit in the last place: the largest
  // relative change of a block, the site where it occurs, and the first site where it exceeds sensitivity_tolerance.
  // It is a property of the star, not of the method: a chain cannot be known better than this from a star in double
  // precision. Filled in by the caller from star_sensitivity(); of the whole chain only.
  double max_star_sensitivity{};
  unsigned int max_star_sensitivity_site{};
  std::optional<unsigned int> sensitive_from_site;
};

template<typename S> struct Chain {
  int channels{};
  unsigned int Nmax{};
  Matrix<S> V;              // the impurity coupling, Theta^(1/2) in the polar gauge
  std::vector<Matrix<S>> E; // the on-site blocks E_0..E_Nmax
  std::vector<Matrix<S>> T; // the hoppings T_0..T_Nmax, the last one out of the chain
  Blocks blocks;            // the blocks of the star, each mapped onto a chain of its own
  ChainGauge gauge{ChainGauge::polar};
  ChainDiagnostics diagnostics;                    // of the whole chain
  std::vector<ChainDiagnostics> block_diagnostics; // one per block, in the order of 'blocks'
};

namespace detail {

// The diagnostics of the whole chain from those of its blocks. Ranks add up site by site; the ratios of eigenvalues
// are taken within each block, since comparing eigenvalues across independent blocks means nothing.
inline ChainDiagnostics merge_diagnostics(const std::vector<ChainDiagnostics> &parts, const unsigned int hoppings) {
  ChainDiagnostics merged;
  merged.hopping_ranks.assign(hoppings, 0);
  bool any_rank = false;
  for (const auto &part : parts) {
    merged.theta_rank += part.theta_rank;
    merged.levels += part.levels;
    merged.coupled_levels += part.coupled_levels;
    for (unsigned int n = 0; n < hoppings; n++) merged.hopping_ranks[n] += part.hopping_ranks[n];
    if (part.theta_rank > 0) {
      merged.theta_condition = any_rank ? std::min(merged.theta_condition, part.theta_condition) : part.theta_condition;
      any_rank               = true;
    }
    merged.min_residual_condition  = std::min(merged.min_residual_condition, part.min_residual_condition);
  }
  merged.min_rank = merged.theta_rank;
  for (unsigned int n = 0; n < hoppings; n++) {
    merged.min_rank = std::min(merged.min_rank, merged.hopping_ranks[n]);
    if (merged.hopping_ranks[n] < merged.theta_rank && !merged.rank_drop_site) merged.rank_drop_site = n;
  }
  return merged;
}

// Move the chain into the nambu gauge: flip the hole component of every second site, U_n = diag(1, (-1)^(n+1)), so
// that V -> V U_0 with U_0 = diag(1, -1), E_n -> U_n E_n U_n and T_n -> U_{n+1} T_n U_n. The blocks must be pairs of
// channels read as (particle, hole).
//
// U_0 is not the identity, and it must not be: only the chain orbitals are free, while the impurity index of V is
// physical, and in Nambu space a normal hybridization v enters as V = diag(v, -conj(v)), since the hole row is
// written with the creation operator. The polar gauge gives V = Theta^(1/2), positive in both slots, which satisfies
// Theta but has the hole coupling of the wrong sign; flipping the hole at the even sites fixes V and leaves the
// relative signs along the chain, which is where the Nambu structure of E_n and T_n lives, untouched.
//
// What comes out is checked against that structure, since the consumer of such a chain stores only the (1,1) and
// (1,2) elements of each block and reconstructs the rest from it.
template<typename S> void apply_nambu_gauge(Chain<S> &chain, const double tolerance) {
  using std::abs;
  for (const auto &block : chain.blocks) {
    if (block.size() != 2)
      throw std::invalid_argument("The nambu gauge needs blocks of two channels, read as particle and hole, but "
                                  + blocks_name({block}) + " has " + std::to_string(block.size()) + ".");
    const auto particle = block[0];
    const auto hole     = block[1];
    const auto flip     = [&](Matrix<S> &m, const int row_sign, const int column_sign) {
      if (column_sign < 0) m(particle, hole) = -m(particle, hole);
      if (row_sign < 0) m(hole, particle) = -m(hole, particle);
      if (row_sign * column_sign < 0) m(hole, hole) = -m(hole, hole);
    };
    const auto column_sign = [](const unsigned int n) { return n % 2 == 0 ? -1 : 1; }; // U_n
    flip(chain.V, 1, column_sign(0));                                                  // V U_0, the impurity index stays
    for (unsigned int n = 0; n <= chain.Nmax; n++) {
      flip(chain.E[n], column_sign(n), column_sign(n));      // U_n E_n U_n
      flip(chain.T[n], column_sign(n + 1), column_sign(n));  // U_{n+1} T_n U_n
    }

    // V(2,2) = -conj(V(1,1)), E(2,2) = -E(1,1) and T(2,2) = -conj(T(1,1)) are what the stored numbers rely on.
    auto &worst  = chain.diagnostics.max_nambu_deviation;
    const auto check = [&worst](const Matrix<S> &m, const S &deviation) {
      const auto scale = static_cast<double>(m.cwiseAbs().maxCoeff());
      if (scale > 0) worst = std::max(worst, static_cast<double>(abs(deviation)) / scale);
    };
    check(chain.V, chain.V(hole, hole) + Eigen::numext::conj(chain.V(particle, particle)));
    check(chain.V, chain.V(hole, particle) + Eigen::numext::conj(chain.V(particle, hole)));
    for (unsigned int n = 0; n <= chain.Nmax; n++) {
      check(chain.E[n], chain.E[n](hole, hole) + chain.E[n](particle, particle));
      check(chain.T[n], chain.T[n](hole, hole) + Eigen::numext::conj(chain.T[n](particle, particle)));
    }
    if (chain.diagnostics.max_nambu_deviation > tolerance)
      throw std::runtime_error("The chain of block " + blocks_name({block})
                               + " does not have the Nambu structure in the nambu gauge: one of V(2,2) + conj(V(1,1)), "
                                 "E(2,2) + E(1,1) and T(2,2) + conj(T(1,1)) reaches "
                               + std::to_string(chain.diagnostics.max_nambu_deviation)
                               + " of the largest element of its block. Is this a superconducting chain?");
  }
  chain.gauge = ChainGauge::nambu;
}

} // namespace detail

// The first site from which the chain samples the untabulated region of the input, the star's
// untabulated_from < |omega| < untabulated_to above the accumulation point of the mesh. The chain resolves ever
// smaller distances from the accumulation point as it goes; the scale of a site is taken as the norm of its hopping,
// and the site samples the region once that scale is below the region's width. Both in the rescaled band. Empty when
// the chain stays above it, when there is no such region, or when it is not known (width 0).
template<typename S> std::optional<unsigned int> first_continued_site(const Chain<S> &chain, const double width) {
  if (!(width > 0.0)) return std::nullopt;
  for (unsigned int n = 0; n < chain.T.size(); n++)
    if (static_cast<double>(chain.T[n].norm()) < width) return n;
  return std::nullopt;
}

// FROM THE STAR TO THE CHAIN
//
// The chain is built by adding the levels of the star one at a time and restoring the form of the chain with plane
// rotations. The Lanczos recursion gives the same chain in exact arithmetic, but loses the orthogonality of its
// vectors to rounding and needs multiprecision arithmetic to be usable. The rotations are unitary transformations of
// the bath, so nothing is lost to cancellation and double precision is enough.
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
// rank it is the pseudo-inverse convention: the chain is zero along the directions that
// are lost, from there on. The rows that are left over stay among those not yet given to a site, where a later site
// may still reach them through M. With full rank throughout, Q is trivial and a step is the SVD of one block.
//
// A singular value counts as zero when its square is below rank_tolerance times the square of the largest, as an
// eigenvalue of the Gram matrix of the couplings would; and all of them do when the largest is rounding on the
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

// The chain of a star, in the arithmetic of the star: each block mapped onto a chain of its own, with exact zeros
// between channels of different blocks.
//
// The levels are given to the blocks by their branch labels, not by where their couplings are nonzero, so that every
// block receives its levels with vanishing coupling as well and has 2*size*(mMAX+1) levels, in proportion to its size
// exactly as the whole star.
template<typename S0> Chain<S0> build_chain(const Star<S0> &star, const ChainOptions &options) {
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
// a star, in whatever arithmetic it runs: multiprecision would give the exact chain of numbers that are not exact.
//
// This measures it: the chain is built again from the star with every energy and coupling moved by one unit in the
// last place, up or down at random, for a few fixed choices, and compared site by site, the hoppings relative to
// their largest element and the on-site blocks on the scale of their site. Always in the polar gauge.
struct StarSensitivity {
  double largest{};
  unsigned int largest_site{};
  std::optional<unsigned int> from_site; // the first site above the tolerance
};

template<typename S0> StarSensitivity star_sensitivity(const Star<S0> &star, ChainOptions options) {
  options.gauge = ChainGauge::polar;
  const auto size      = [](const Matrix<S0> &m) { return m.size() ? m.cwiseAbs().maxCoeff() : 0.0; };
  const auto reference = build_chain(star, options);
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
    const auto chain = build_chain(other, options);
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
