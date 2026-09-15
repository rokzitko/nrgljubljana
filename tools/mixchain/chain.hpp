// Channel-mixing discretization for NRG
// ** Block Lanczos: from the star to the Wilson chain

#ifndef _mixchain_chain_hpp_
#define _mixchain_chain_hpp_

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

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
  unsigned int Nmax{0}; // the chain has the sites 0..Nmax
  // A matrix that is inverted on the way, Theta or the Gram matrix of a residual block, counts as singular when its
  // smallest eigenvalue is below this fraction of its largest. Relative, so that it means the same at every
  // precision.
  double breakdown_tolerance{1e-20};
};

struct ChainDiagnostics {
  double theta_condition{};         // the smallest eigenvalue of Theta over its largest
  double max_antihermitian{};       // the largest anti-Hermitian part removed from an on-site block, relative to it
  double max_reorthogonalization{}; // the largest component along the earlier blocks removed from a residual, relative
  // The smallest lambda_min/lambda_max of the Gram matrix R^dag R of a residual over the whole chain. A breakdown is
  // raised only below the tolerance; this says how close a run that did not break down came to it.
  double min_residual_condition{1.0};
};

template<typename S> struct Chain {
  int channels{};
  unsigned int Nmax{};
  Matrix<S> V;              // the impurity coupling, Theta^(1/2) in the polar gauge
  std::vector<Matrix<S>> E; // the on-site blocks E_0..E_Nmax
  std::vector<Matrix<S>> T; // the hoppings T_0..T_{Nmax-1}
  ChainDiagnostics diagnostics;
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

// The star in the arithmetic of the recursion: the level energies, which are real, and the starting block A.
template<typename S> struct WideStar {
  std::vector<real_type<S>> energies; // E_k
  Matrix<S> start;                    // A[k,i] = conj(v_{k,i})
};

template<typename S, typename StarScalar> auto to_wide(const Star<StarScalar> &star) {
  static_assert(is_complex_v<S> || !is_complex_v<StarScalar>, "a complex star needs a complex chain");
  const auto levels   = star.levels.size();
  const auto channels = static_cast<Eigen::Index>(star.channels);

  WideStar<S> wide;
  wide.energies.reserve(levels);
  wide.start = Matrix<S>::Zero(static_cast<Eigen::Index>(levels), channels);
  for (std::size_t k = 0; k < levels; k++) {
    const auto &level = star.levels[k];
    wide.energies.push_back(static_cast<real_type<S>>(level.energy));
    for (Eigen::Index i = 0; i < channels; i++)
      wide.start(static_cast<Eigen::Index>(k), i) = convert_scalar<S>(Eigen::numext::conj(level.coupling(i)));
  }
  return wide;
}

// Raised when a matrix that the recursion has to invert is singular: Theta itself, when Gamma is rank deficient over
// the whole band, or the Gram matrix of a residual block part-way down the chain, when the Krylov space has run out.
class ChainBreakdown : public std::runtime_error {
 private:
  std::optional<unsigned int> site_; // empty for Theta, before the first site
  int rank_{};
  int channels_{};
  double condition_{};

 public:
  ChainBreakdown(const std::string &message, const std::optional<unsigned int> site, const int rank,
                 const int channels, const double condition)
    : std::runtime_error(message), site_(site), rank_(rank), channels_(channels), condition_(condition) {}

  [[nodiscard]] auto site() const { return site_; }
  [[nodiscard]] auto rank() const { return rank_; }
  [[nodiscard]] auto channels() const { return channels_; }
  [[nodiscard]] auto condition() const { return condition_; }
};

namespace detail {

// The square root of a Hermitian positive semidefinite matrix, the inverse of that root, and the ratio of its
// smallest to its largest eigenvalue.
template<typename S> struct HermitianRoot {
  Matrix<S> root;
  Matrix<S> inverse;
  double condition{};
};

// The breakdown tolerance actually applied: the requested one, but never below what rounding alone produces in this
// arithmetic. In double precision an exactly singular matrix still shows a smallest eigenvalue of about 1e-16 of the
// largest, which a fixed tolerance of 1e-20 would take for a regular one; at 800 digits the requested tolerance
// governs.
template<typename S> double effective_tolerance(const double tolerance) {
  const auto epsilon = static_cast<double>(std::numeric_limits<real_type<S>>::epsilon());
  return std::max(tolerance, 1000.0 * epsilon);
}

// G^(1/2) and G^(-1/2) from a single eigendecomposition, for the Gram matrix G of a block: Theta at the start, where
// the root is the impurity coupling V and the inverse turns A into the first Lanczos block, and R^dag R at every
// step, where the root is the hopping T_n and the inverse normalizes the residual R into the next block. 'site'
// names where this happens for the breakdown message, and is empty for Theta.
template<typename S>
auto hermitian_root(const Matrix<S> &gram, const double tolerance, const std::optional<unsigned int> site) {
  const auto n = gram.rows();
  // G is Hermitian mathematically, but its (i,j) and (j,i) elements are different sums.
  const ColumnMajor<S> symmetric = (make_scalar<S>(0.5, 0) * (gram + gram.adjoint())).eval();
  Eigen::SelfAdjointEigenSolver<ColumnMajor<S>> solver(symmetric);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Diagonalization of a Gram matrix failed in the block Lanczos recursion.");

  const auto &lambda   = solver.eigenvalues(); // real and ascending
  const auto largest   = lambda(n - 1);
  const auto condition = largest > 0 ? static_cast<double>(lambda(0) / largest) : 0.0;
  const auto threshold = effective_tolerance<S>(tolerance);

  if (!(largest > 0) || condition < threshold) {
    int rank = 0;
    for (Eigen::Index i = 0; i < n; i++)
      if (largest > 0 && static_cast<double>(lambda(i) / largest) >= threshold) rank++;
    const auto channels = static_cast<int>(n);
    std::ostringstream ratio; // std::to_string would print a ratio of 1e-17 as 0.000000
    ratio << std::setprecision(3) << condition;
    const auto figures = "rank " + std::to_string(rank) + " of " + std::to_string(channels) + ", smallest eigenvalue "
                         + ratio.str() + " times the largest";
    if (!site)
      throw ChainBreakdown("Theta is singular (" + figures + "): Gamma is rank deficient over the whole band, so the "
                           "bath couples to fewer combinations of the impurity orbitals than there are channels.",
                           site, rank, channels, condition);
    throw ChainBreakdown("The block Lanczos recursion broke down at site " + std::to_string(*site) + " (" + figures
                           + "): the Krylov space of the star is exhausted. The star needs at least channels*(Nmax+1) "
                             "levels with nonzero coupling.",
                         site, rank, channels, condition);
  }

  const auto roots = lambda.cwiseSqrt().eval();
  HermitianRoot<S> result;
  // The root becomes V or T_n, which are Hermitian; the product U diag U^dag is so only up to rounding.
  const Matrix<S> root = solver.eigenvectors() * roots.template cast<S>().asDiagonal() * solver.eigenvectors().adjoint();
  result.root          = make_scalar<S>(0.5, 0) * (root + root.adjoint());
  result.inverse   = solver.eigenvectors() * roots.cwiseInverse().template cast<S>().asDiagonal()
                   * solver.eigenvectors().adjoint();
  result.condition = condition;
  return result;
}

// The component of R along the earlier Lanczos blocks, sum_m Q_m (Q_m^dag R).
template<typename S> Matrix<S> component_along(const Matrix<S> &residual, const std::vector<Matrix<S>> &blocks) {
  Matrix<S> component = Matrix<S>::Zero(residual.rows(), residual.cols());
  for (const auto &block : blocks) component += block * (block.adjoint() * residual);
  return component;
}

} // namespace detail

// Block Lanczos from the star to the chain, in the polar gauge:
//
//   start:   Theta = A^dag A,   V = Theta^(1/2),   Q_0 = A Theta^(-1/2),
//   site n:  E_n = Q_n^dag H Q_n,
//            R = H Q_n - Q_n E_n - Q_{n-1} T_{n-1}^dag,   then reorthogonalized against every earlier block,
//            T_n = (R^dag R)^(1/2),   Q_{n+1} = R (R^dag R)^(-1/2).
//
// H is the bath Hamiltonian, diagonal in the star levels. Every block is kept, because the reorthogonalization needs
// them all. 'lanczos_blocks', if given, receives Q_0..Q_Nmax, which only a test has a use for: stacked side by side
// they are the unitary that maps the star onto the chain.
template<typename S>
auto build_chain(const WideStar<S> &star, const ChainOptions &options,
                 std::vector<Matrix<S>> *lanczos_blocks = nullptr) {
  using std::sqrt; // for double; the wide types are found by argument-dependent lookup
  const auto levels   = star.start.rows();
  const auto channels = star.start.cols();
  if (channels < 1) throw std::invalid_argument("The star has no channels.");
  if (options.Nmax < 1) throw std::invalid_argument("Nmax must be greater than 0.");
  // A chain of Nmax+1 sites spans a Krylov space of dimension channels*(Nmax+1), which the star must be able to hold.
  // Levels with vanishing coupling can make the space that is actually reached smaller still; that shows up as a
  // breakdown at a definite site.
  const auto needed = channels * static_cast<Eigen::Index>(options.Nmax + 1);
  if (levels < needed)
    throw std::invalid_argument("The star has " + std::to_string(levels) + " levels, but a chain of "
                                + std::to_string(options.Nmax + 1) + " sites with " + std::to_string(channels)
                                + " channels needs at least " + std::to_string(needed)
                                + ". Increase mMAX or decrease Nmax.");

  Chain<S> chain;
  chain.channels = static_cast<int>(channels);
  chain.Nmax     = options.Nmax;
  chain.E.reserve(options.Nmax + 1);
  chain.T.reserve(options.Nmax);
  auto &diagnostics = chain.diagnostics;

  // The energies are real, but a same-type product with the blocks needs them in S.
  Vector<S> energies(levels);
  for (Eigen::Index k = 0; k < levels; k++) energies(k) = make_scalar<S>(star.energies[static_cast<std::size_t>(k)], 0);

  const Matrix<S> theta = star.start.adjoint() * star.start;
  const auto start      = detail::hermitian_root<S>(theta, options.breakdown_tolerance, std::nullopt);
  chain.V                     = start.root;
  diagnostics.theta_condition = start.condition;

  std::vector<Matrix<S>> blocks;
  blocks.reserve(options.Nmax + 1);
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

    if (n == options.Nmax) break; // the last site has no outgoing hopping

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
    const auto step      = detail::hermitian_root<S>(gram, options.breakdown_tolerance, n);
    diagnostics.min_residual_condition = std::min(diagnostics.min_residual_condition, step.condition);
    chain.T.push_back(step.root);
    blocks.push_back(residual * step.inverse);
  }
  if (lanczos_blocks) *lanczos_blocks = std::move(blocks);
  return chain;
}

// The same from a star in double precision, widened to S first.
template<typename S, typename StarScalar> auto build_chain(const Star<StarScalar> &star, const ChainOptions &options) {
  return build_chain<S>(to_wide<S>(star), options);
}

// The chain in another arithmetic: narrowing the result of the recursion to double for writing it out, or widening
// it in the tests.
template<typename To, typename From> auto convert_chain(const Chain<From> &chain) {
  Chain<To> result;
  result.channels    = chain.channels;
  result.Nmax        = chain.Nmax;
  result.V           = convert_matrix<To>(chain.V);
  result.diagnostics = chain.diagnostics;
  result.E.reserve(chain.E.size());
  result.T.reserve(chain.T.size());
  for (const auto &block : chain.E) result.E.push_back(convert_matrix<To>(block));
  for (const auto &block : chain.T) result.T.push_back(convert_matrix<To>(block));
  return result;
}

} // namespace NRG::MixChain

#endif
