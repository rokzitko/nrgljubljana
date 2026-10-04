// Channel-mixing discretization for NRG
// ** Plane rotations: from the star to the Wilson chain, in double precision

#ifndef _mixchain_chain_rkpw_hpp_
#define _mixchain_chain_rkpw_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include <star-to-chain.hpp>

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
// Rutishauser-Kahan-Pal-Walker rotations that nrgchain uses. Blocks of several channels are not handled yet.

namespace detail {

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
  std::stable_sort(levels.begin(), levels.end(), [](const StarLevel<S0> *a, const StarLevel<S0> *b) {
    const auto key = [](const StarLevel<S0> *level) {
      return std::tuple(level->m, level->sign != Sign::POS, level->branch);
    };
    return key(a) < key(b);
  });

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
    if (size != 1)
      throw std::invalid_argument("tridiag_method=rkpw handles only blocks of one channel so far, but "
                                  + blocks_name({block}) + " has " + std::to_string(size)
                                  + "; use tridiag_method=lanczos.");

    // As in build_chain(): the levels of a block are those of its branches, with or without coupling.
    std::vector<const StarLevel<S0> *> levels;
    for (const auto &level : star.levels)
      if (chain.blocks.size() == 1 || (level.branch >= offset && level.branch < offset + size))
        levels.push_back(&level);

    const auto channel = static_cast<Eigen::Index>(block.front());
    const auto piece   = detail::scalar_block_chain(std::move(levels), channel, options);
    chain.V(channel, channel) = make_scalar<S0>(piece.V, 0);
    for (unsigned int n = 0; n <= options.Nmax; n++) {
      chain.E[n](channel, channel) = make_scalar<S0>(piece.zeta[n], 0);
      chain.T[n](channel, channel) = make_scalar<S0>(piece.xi[n], 0);
    }
    chain.block_diagnostics.push_back(piece.diagnostics);
    offset += size;
  }
  chain.diagnostics = detail::merge_diagnostics(chain.block_diagnostics, options.Nmax + 1);
  if (options.gauge == ChainGauge::nambu) detail::apply_nambu_gauge(chain, options.nambu_tolerance);
  return chain;
}

} // namespace NRG::MixChain

#endif
