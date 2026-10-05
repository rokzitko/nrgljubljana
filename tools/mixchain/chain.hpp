// Channel-mixing discretization for NRG
// ** The Wilson chain: its blocks, gauges and diagnostics

#ifndef _mixchain_chain_hpp_
#define _mixchain_chain_hpp_

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

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

// The gauge the chain is written in. Lanczos fixes each site only up to a unitary rotation of its N orbitals, and the
// blocks of V, E_n and T_n transform together, so every gauge describes the same bath.
//
//   polar: V and every T_n Hermitian positive semidefinite, the matrix analogue of choosing xi_n > 0. Nothing is
//          assumed about what the channels mean, and every element is written, so a consumer that reads the whole
//          matrix (a pol2x2 template of nrg, or one coefficient set per channel) can use this as it is.
//   nambu: for blocks of two channels read as (particle, hole). The consumer of a superconducting chain stores only
//          xi = T(1,1), sckappa = T(1,2), zeta = E(1,1), scdelta = E(1,2) and reconstructs the rest from the Nambu
//          structure, so the chain must be in the gauge where that structure holds: E(2,2) = -E(1,1) and
//          T(2,2) = -conj(T(1,1)). The polar gauge is not: it absorbs the sign of the hole component into the Lanczos
//          block, which turns a constant gap into one alternating along the chain. Flipping the hole component of
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
  // An eigenvalue of a Gram matrix, of Theta or of R^dag R for a residual block, counts as zero when it is below this
  // fraction of the largest one. Relative, so that it means the same at every precision.
  double rank_tolerance{1e-20};
  ChainGauge gauge{ChainGauge::polar};
  // How far a block may depart from the Nambu structure before the nambu gauge refuses the chain, relative to the
  // largest element of that block.
  double nambu_tolerance{1e-8};
  // The chain counts as no longer determined by the star from the first site that moves by more than this when every
  // number of the star is changed by one unit in the last place; see star_sensitivity() in chain_rkpw.hpp.
  double sensitivity_tolerance{1e-10};
};

struct ChainDiagnostics {
  int theta_rank{};                 // the number of combinations of the impurity orbitals that couple to the bath
  double theta_condition{};         // the smallest nonzero eigenvalue of Theta over its largest; 0 if Theta is zero
  double max_antihermitian{};       // the largest anti-Hermitian part removed from an on-site block, relative to it
  double max_reorthogonalization{}; // the largest component along the earlier blocks removed from a residual, relative
  // The smallest lambda_min/lambda_max over the nonzero eigenvalues of the Gram matrices R^dag R along the chain. It
  // says how close a direction came to being counted as zero.
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
    merged.max_antihermitian       = std::max(merged.max_antihermitian, part.max_antihermitian);
    merged.max_reorthogonalization = std::max(merged.max_reorthogonalization, part.max_reorthogonalization);
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

} // namespace NRG::MixChain

#endif
