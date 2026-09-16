// Channel-mixing discretization for NRG
// ** Block Lanczos: from the star to the Wilson chain

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

struct ChainOptions {
  // The chain has the sites 0..Nmax, and one hopping per site, T_0..T_Nmax: the last leads out of the chain and is
  // there because the coefficient tables of nrg are indexed 0..Nmax, as nrgchain writes xi.dat and zeta.dat.
  unsigned int Nmax{0};
  // An eigenvalue of a Gram matrix, of Theta or of R^dag R for a residual block, counts as zero when it is below this
  // fraction of the largest one. Relative, so that it means the same at every precision.
  double rank_tolerance{1e-20};
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
};

template<typename S> struct Chain {
  int channels{};
  unsigned int Nmax{};
  Matrix<S> V;              // the impurity coupling, Theta^(1/2) in the polar gauge
  std::vector<Matrix<S>> E; // the on-site blocks E_0..E_Nmax
  std::vector<Matrix<S>> T; // the hoppings T_0..T_Nmax, the last one out of the chain
  Blocks blocks;            // the blocks of the star, each mapped onto a chain of its own
  ChainDiagnostics diagnostics;                    // of the whole chain
  std::vector<ChainDiagnostics> block_diagnostics; // one per block, in the order of 'blocks'
};

// A scalar of type To from one of type From, through its real and imaginary parts: between double and the wide
// types of precision.hpp, in either direction. Converting a complex scalar to a real one would silently drop data,
// so it is refused at compile time.
template<typename To, typename From> To convert_scalar(const From &x) {
  static_assert(is_complex_v<To> || !is_complex_v<From>, "a complex scalar cannot be converted to a real one");
  const auto re = static_cast<real_type<To>>(Eigen::numext::real(x));
  const auto im = static_cast<real_type<To>>(Eigen::numext::imag(x));
  return make_scalar<To>(re, im);
}

// Eigen's cast converts between double and the wide types in both directions; the assertion keeps the guard of
// convert_scalar().
template<typename To, typename From> Matrix<To> convert_matrix(const Matrix<From> &m) {
  static_assert(is_complex_v<To> || !is_complex_v<From>, "a complex matrix cannot be converted to a real one");
  return m.template cast<To>();
}

// The star in the arithmetic of the recursion: the level energies, which are real, and the starting block A, with the
// blocks and the branch label of every level, which says which block the level belongs to.
template<typename S> struct WideStar {
  std::vector<real_type<S>> energies; // E_k
  Matrix<S> start;                    // A[k,i] = conj(v_{k,i})
  Blocks blocks;                      // empty for a single block of all channels
  std::vector<int> branches;          // the branch label of every level
};

template<typename S, typename StarScalar> auto to_wide(const Star<StarScalar> &star) {
  static_assert(is_complex_v<S> || !is_complex_v<StarScalar>, "a complex star needs a complex chain");
  const auto levels   = star.levels.size();
  const auto channels = static_cast<Eigen::Index>(star.channels);

  WideStar<S> wide;
  wide.blocks = star.blocks;
  wide.energies.reserve(levels);
  wide.branches.reserve(levels);
  wide.start = Matrix<S>::Zero(static_cast<Eigen::Index>(levels), channels);
  for (std::size_t k = 0; k < levels; k++) {
    const auto &level = star.levels[k];
    wide.energies.push_back(static_cast<real_type<S>>(level.energy));
    wide.branches.push_back(level.branch);
    for (Eigen::Index i = 0; i < channels; i++)
      wide.start(static_cast<Eigen::Index>(k), i) = convert_scalar<S>(Eigen::numext::conj(level.coupling(i)));
  }
  return wide;
}

namespace detail {

// The square root of a Hermitian positive semidefinite matrix and its pseudo-inverse, with the number of eigenvalues
// kept as nonzero and the ratio of the smallest kept one to the largest.
template<typename S> struct HermitianRoot {
  Matrix<S> root;
  Matrix<S> inverse;
  int rank{};
  double condition{};
};

// The rank tolerance actually applied: the requested one, but never below what rounding alone produces in this
// arithmetic. In double precision an exactly singular matrix still shows a smallest eigenvalue of about 1e-16 of the
// largest, which a fixed tolerance of 1e-20 would take for a nonzero one; at 800 digits the requested tolerance
// governs.
template<typename S> double effective_tolerance(const double tolerance) {
  const auto epsilon = static_cast<double>(std::numeric_limits<real_type<S>>::epsilon());
  return std::max(tolerance, 1000.0 * epsilon);
}

// G^(1/2) and the pseudo-inverse G^(+1/2) from a single eigendecomposition, for the Gram matrix G of a block: Theta at
// the start, where the root is the impurity coupling V and the inverse turns A into the first Lanczos block, and
// R^dag R at every step, where the root is the hopping T_n and the inverse normalizes the residual R into the next
// block.
//
// Eigenvalues below the tolerance times the largest are set to zero in both, so the directions they belong to drop
// out of the next block and their part of the chain is zero from there on. A relative test cannot tell when every
// direction is rounding, as when the Krylov space of a single channel runs out; so the largest eigenvalue is also
// compared with 'scale', the squared norm of the block before its projection, H Q_n for a residual. Rounding makes
// the residual of order epsilon times that norm, and a hopping that is really there is far above it. For Theta the
// scale is zero: Theta is formed directly from the star, and is zero only when every coupling is.
template<typename S> auto hermitian_root(const Matrix<S> &gram, const real_type<S> &scale, const double tolerance) {
  using std::sqrt; // for double; the wide types are found by argument-dependent lookup
  const auto n = gram.rows();
  // G is Hermitian mathematically, but its (i,j) and (j,i) elements are different sums.
  const ColumnMajor<S> symmetric = (make_scalar<S>(0.5, 0) * (gram + gram.adjoint())).eval();
  Eigen::SelfAdjointEigenSolver<ColumnMajor<S>> solver(symmetric);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Diagonalization of a Gram matrix failed in the block Lanczos recursion.");

  const auto &lambda      = solver.eigenvalues(); // real and ascending
  const auto largest      = lambda(n - 1);
  const auto epsilon      = std::numeric_limits<real_type<S>>::epsilon();
  const auto rounding     = real_type<S>(1000) * epsilon;
  const bool all_rounding = !(largest > 0) || largest <= rounding * rounding * scale;
  const auto threshold    = real_type<S>(effective_tolerance<S>(tolerance)) * largest;

  Vector<real_type<S>> roots          = Vector<real_type<S>>::Zero(n);
  Vector<real_type<S>> inverse_roots  = Vector<real_type<S>>::Zero(n);
  HermitianRoot<S> result;
  if (!all_rounding) {
    for (Eigen::Index i = 0; i < n; i++) {
      if (lambda(i) < threshold) continue;
      if (result.rank == 0) result.condition = static_cast<double>(lambda(i) / largest); // the smallest one kept
      roots(i)         = sqrt(lambda(i));
      inverse_roots(i) = 1 / roots(i);
      result.rank++;
    }
  }
  const auto &U = solver.eigenvectors();
  // The root becomes V or T_n, which are Hermitian; the product U diag U^dag is so only up to rounding.
  const Matrix<S> root = U * roots.template cast<S>().asDiagonal() * U.adjoint();
  result.root          = make_scalar<S>(0.5, 0) * (root + root.adjoint());
  result.inverse       = U * inverse_roots.template cast<S>().asDiagonal() * U.adjoint();
  return result;
}

// The component of R along the earlier Lanczos blocks, sum_m Q_m (Q_m^dag R).
template<typename S> Matrix<S> component_along(const Matrix<S> &residual, const std::vector<Matrix<S>> &blocks) {
  Matrix<S> component = Matrix<S>::Zero(residual.rows(), residual.cols());
  for (const auto &block : blocks) component += block * (block.adjoint() * residual);
  return component;
}

} // namespace detail

namespace detail {

// Block Lanczos from the star to the chain, in the polar gauge:
//
//   start:   Theta = A^dag A,   V = Theta^(1/2),   Q_0 = A Theta^(-1/2),
//   site n:  E_n = Q_n^dag H Q_n,
//            R = H Q_n - Q_n E_n - Q_{n-1} T_{n-1}^dag,   then reorthogonalized against every earlier block,
//            T_n = (R^dag R)^(1/2),   Q_{n+1} = R (R^dag R)^(-1/2).
//
// H is the bath Hamiltonian, diagonal in the star levels. The recursion runs to n = Nmax, so it produces the sites
// 0..Nmax and the hoppings T_0..T_Nmax: one hopping per site, the last leading out of the chain. That is the length
// nrg reads, and what nrgchain writes into xi.dat and zeta.dat.
//
// Every block is kept, because the reorthogonalization needs them all. 'lanczos_blocks', if given, receives
// Q_0..Q_{Nmax+1}, which only a test has a use for: stacked side by side they are the unitary that maps the star
// onto the chain.
//
// The star is taken as a single block; build_chain() splits it first.
template<typename S>
Chain<S> block_lanczos(const WideStar<S> &star, const ChainOptions &options, std::vector<Matrix<S>> *lanczos_blocks) {
  using std::sqrt; // for double; the wide types are found by argument-dependent lookup
  const auto levels   = star.start.rows();
  const auto channels = star.start.cols();
  if (channels < 1) throw std::invalid_argument("The star has no channels.");
  if (options.Nmax < 1) throw std::invalid_argument("Nmax must be greater than 0.");
  // The sites 0..Nmax and the hopping out of the last one span a Krylov space of dimension channels*(Nmax+2), which
  // the star must be able to hold. Levels with vanishing coupling can make the space that is actually reached smaller
  // still; that shows up as a drop in the rank of a hopping, rank_drop_site in the diagnostics.
  const auto needed = channels * static_cast<Eigen::Index>(options.Nmax + 2);
  if (levels < needed)
    throw std::invalid_argument("The star has " + std::to_string(levels) + " levels, but a chain of "
                                + std::to_string(options.Nmax + 1) + " sites with " + std::to_string(channels)
                                + " channels needs at least " + std::to_string(needed)
                                + ", one block more than the sites, for the hopping out of the last site"
                                + ". Increase mMAX or decrease Nmax.");

  Chain<S> chain;
  chain.channels = static_cast<int>(channels);
  chain.Nmax     = options.Nmax;
  chain.E.reserve(options.Nmax + 1);
  chain.T.reserve(options.Nmax + 1);
  auto &diagnostics  = chain.diagnostics;
  diagnostics.levels = static_cast<int>(levels);
  for (Eigen::Index k = 0; k < levels; k++)
    if ((star.start.row(k).array() != S(0)).any()) diagnostics.coupled_levels++;

  // The energies are real, but a same-type product with the blocks needs them in S.
  Vector<S> energies(levels);
  for (Eigen::Index k = 0; k < levels; k++) energies(k) = make_scalar<S>(star.energies[static_cast<std::size_t>(k)], 0);

  const Matrix<S> theta = star.start.adjoint() * star.start;
  const auto start      = detail::hermitian_root<S>(theta, real_type<S>(0), options.rank_tolerance);
  chain.V                     = start.root;
  diagnostics.theta_rank      = start.rank;
  diagnostics.theta_condition = start.condition;
  diagnostics.min_rank        = start.rank;

  std::vector<Matrix<S>> blocks;
  blocks.reserve(options.Nmax + 2);
  blocks.push_back(star.start * start.inverse);

  const auto half = make_scalar<S>(0.5, 0);
  for (unsigned int n = 0; n <= options.Nmax; n++) {
    const Matrix<S> hq = energies.asDiagonal() * blocks[n];

    // E_n is Hermitian mathematically, but its (i,j) and (j,i) elements are different sums.
    const Matrix<S> onsite    = blocks[n].adjoint() * hq;
    const Matrix<S> hermitian = half * (onsite + onsite.adjoint());
    const auto onsite_norm    = hermitian.norm();
    if (onsite_norm > 0)
      diagnostics.max_antihermitian =
        std::max(diagnostics.max_antihermitian, static_cast<double>((onsite - hermitian).norm() / onsite_norm));
    chain.E.push_back(hermitian);

    Matrix<S> residual = hq - blocks[n] * chain.E[n];
    if (n > 0) residual -= blocks[n - 1] * chain.T[n - 1].adjoint();

    // Full reorthogonalization. The three-term recurrence is orthogonal to the earlier blocks only up to rounding, and
    // the loss accumulates along the chain. A second pass is made when the first removed a large part of the
    // residual, which is when a single classical Gram-Schmidt pass is known to be insufficient ("twice is enough",
    // Kahan and Parlett).
    const auto before    = residual.squaredNorm();
    const Matrix<S> once = detail::component_along(residual, blocks);
    residual -= once;
    if (before > 0)
      diagnostics.max_reorthogonalization =
        std::max(diagnostics.max_reorthogonalization, static_cast<double>(sqrt(once.squaredNorm() / before)));
    if (residual.squaredNorm() < before / 2) residual -= detail::component_along(residual, blocks);

    const Matrix<S> gram = residual.adjoint() * residual;
    const auto step      = detail::hermitian_root<S>(gram, hq.squaredNorm(), options.rank_tolerance);
    if (step.rank > 0)
      diagnostics.min_residual_condition = std::min(diagnostics.min_residual_condition, step.condition);
    diagnostics.min_rank = std::min(diagnostics.min_rank, step.rank);
    diagnostics.hopping_ranks.push_back(step.rank);
    if (step.rank < diagnostics.theta_rank && !diagnostics.rank_drop_site) diagnostics.rank_drop_site = n;
    chain.T.push_back(step.root);
    blocks.push_back(residual * step.inverse);
  }
  if (lanczos_blocks) *lanczos_blocks = std::move(blocks);
  return chain;
}

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

} // namespace detail

// The chain of a star: block Lanczos for each of its blocks, assembled into N x N blocks with exact zeros between
// channels of different blocks. With a single block this is block_lanczos() itself.
//
// The levels are given to the blocks by their branch labels, not by where their couplings are nonzero, so that every
// block receives its levels with vanishing coupling as well and has 2*size*(mMAX+1) levels, in proportion to its size
// exactly as the whole star. Within a block the levels keep their order in the star.
//
// 'lanczos_blocks', if given, receives Q_0..Q_{Nmax+1} of the whole star: the blocks of the parts, placed at the rows
// of their levels and the columns of their channels.
template<typename S>
auto build_chain(const WideStar<S> &star, const ChainOptions &options,
                 std::vector<Matrix<S>> *lanczos_blocks = nullptr) {
  const auto levels   = star.start.rows();
  const auto channels = static_cast<int>(star.start.cols());
  if (star.blocks.size() <= 1) {
    auto chain = detail::block_lanczos(star, options, lanczos_blocks);
    if (star.blocks.empty()) {
      chain.blocks.emplace_back(static_cast<std::size_t>(channels));
      std::iota(chain.blocks.front().begin(), chain.blocks.front().end(), 0);
    } else {
      chain.blocks = star.blocks;
    }
    chain.block_diagnostics = {chain.diagnostics};
    return chain;
  }
  if (star.branches.size() != static_cast<std::size_t>(levels))
    throw std::invalid_argument("A star with blocks needs the branch label of every level.");

  Chain<S> chain;
  chain.channels = channels;
  chain.Nmax     = options.Nmax;
  chain.blocks   = star.blocks;
  chain.V        = Matrix<S>::Zero(channels, channels);
  chain.E.assign(options.Nmax + 1, Matrix<S>::Zero(channels, channels));
  chain.T.assign(options.Nmax + 1, Matrix<S>::Zero(channels, channels));
  if (lanczos_blocks) lanczos_blocks->assign(options.Nmax + 2, Matrix<S>::Zero(levels, channels));

  int offset = 0; // the first branch label of the current block
  for (const auto &block : star.blocks) {
    const auto size = static_cast<Eigen::Index>(block.size());
    std::vector<Eigen::Index> rows; // the levels of this block, in their order in the star
    for (Eigen::Index k = 0; k < levels; k++) {
      const auto branch = star.branches[static_cast<std::size_t>(k)];
      if (branch >= offset && branch < offset + size) rows.push_back(k);
    }

    WideStar<S> part;
    part.start = Matrix<S>(static_cast<Eigen::Index>(rows.size()), size);
    part.energies.reserve(rows.size());
    for (std::size_t r = 0; r < rows.size(); r++) {
      part.energies.push_back(star.energies[static_cast<std::size_t>(rows[r])]);
      for (Eigen::Index i = 0; i < size; i++)
        part.start(static_cast<Eigen::Index>(r), i) = star.start(rows[r], block[static_cast<std::size_t>(i)]);
    }

    std::vector<Matrix<S>> part_blocks;
    const auto piece = detail::block_lanczos(part, options, lanczos_blocks ? &part_blocks : nullptr);

    const auto place = [&block, size](Matrix<S> &target, const Matrix<S> &source) {
      for (Eigen::Index i = 0; i < size; i++)
        for (Eigen::Index j = 0; j < size; j++)
          target(block[static_cast<std::size_t>(i)], block[static_cast<std::size_t>(j)]) = source(i, j);
    };
    place(chain.V, piece.V);
    for (unsigned int n = 0; n <= options.Nmax; n++) place(chain.E[n], piece.E[n]);
    for (unsigned int n = 0; n <= options.Nmax; n++) place(chain.T[n], piece.T[n]);
    if (lanczos_blocks)
      for (unsigned int n = 0; n <= options.Nmax + 1; n++)
        for (std::size_t r = 0; r < rows.size(); r++)
          for (Eigen::Index i = 0; i < size; i++)
            (*lanczos_blocks)[n](rows[r], block[static_cast<std::size_t>(i)]) =
              part_blocks[n](static_cast<Eigen::Index>(r), i);

    chain.block_diagnostics.push_back(piece.diagnostics);
    offset += static_cast<int>(size);
  }
  chain.diagnostics = detail::merge_diagnostics(chain.block_diagnostics, options.Nmax + 1);
  return chain;
}

// The same from a star in double precision, widened to S first.
template<typename S, typename StarScalar> auto build_chain(const Star<StarScalar> &star, const ChainOptions &options) {
  return build_chain<S>(to_wide<S>(star), options);
}

// The first site from which the chain is built on the continued density: the scale of a site is taken as the norm of
// its hopping, and 'innermost' is the smallest frequency the input tabulates, both in the rescaled band. Empty when
// the chain stays above it, or when the frequency is not known.
template<typename S>
std::optional<unsigned int> first_continued_site(const Chain<S> &chain, const double innermost) {
  if (!(innermost > 0.0)) return std::nullopt;
  for (unsigned int n = 0; n < chain.T.size(); n++)
    if (static_cast<double>(chain.T[n].norm()) < innermost) return n;
  return std::nullopt;
}

// The chain in another arithmetic: narrowing the result of the recursion to double for writing it out, or widening
// it in the tests.
template<typename To, typename From> auto convert_chain(const Chain<From> &chain) {
  Chain<To> result;
  result.channels    = chain.channels;
  result.Nmax        = chain.Nmax;
  result.V           = convert_matrix<To>(chain.V);
  result.blocks      = chain.blocks;
  result.diagnostics = chain.diagnostics;
  result.block_diagnostics = chain.block_diagnostics;
  result.E.reserve(chain.E.size());
  result.T.reserve(chain.T.size());
  for (const auto &block : chain.E) result.E.push_back(convert_matrix<To>(block));
  for (const auto &block : chain.T) result.T.push_back(convert_matrix<To>(block));
  return result;
}

} // namespace NRG::MixChain

#endif
