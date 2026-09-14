// Channel-mixing discretization for NRG
// ** Eigenvalue branches of Gamma(omega)

#ifndef _mixchain_branches_hpp_
#define _mixchain_branches_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "load.hpp"
#include "mesh.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// How the eigenvalues of Gamma at neighbouring nodes are matched up into branches.
enum class BranchOrdering {
  tracked, // by eigenvector overlap, so that a branch survives a crossing
  sorted   // by eigenvalue, descending
};

inline auto branch_ordering_from_string(const std::string &value) {
  if (value == "tracked") return BranchOrdering::tracked;
  if (value == "sorted") return BranchOrdering::sorted;
  throw std::invalid_argument("branch_ordering must be either 'tracked' or 'sorted'.");
}

inline auto branch_ordering_name(const BranchOrdering ordering) {
  return ordering == BranchOrdering::tracked ? std::string("tracked") : std::string("sorted");
}

struct BranchOptions {
  BranchOrdering ordering{BranchOrdering::tracked};
  double degeneracy_tolerance{1e-8}; // relative to the largest eigenvalue at that node
  double psd_tolerance{1e-8};        // a more negative eigenvalue is an error rather than rounding
};

// Gamma(omega) = sum_a rho_a(omega) u_a(omega) u_a(omega)^dagger on one frequency branch, with the labels a assigned
// consistently across the nodes.
template<typename S> struct BranchDecomposition {
  int channels{};
  std::vector<double> omega;      // the nodes of the input branch
  std::vector<Vec> density;       // density[a] = the (omega, rho_a) table, ready for TabulatedDensity
  std::vector<Matrix<S>> vectors; // vectors[k].col(a) = u_a(omega_k)
  std::vector<double> crossings;  // nodes where the tracked order differs from the sorted one

  auto nodes() const { return omega.size(); }
};

namespace detail {

// Eigen's decompositions are written for its natural storage order; Matrix<S> is row-major to match the alias of the
// nrg core, so the solver is given a column-major copy.
template<typename S> using ColumnMajor = Eigen::Matrix<S, -1, -1>;

// One node: the eigenvalues in descending order and the corresponding orthonormal eigenvectors as columns.
template<typename S>
auto diagonalize(const Matrix<S> &gamma, const double omega, const BranchOptions &options) {
  const ColumnMajor<S> input = gamma;
  Eigen::SelfAdjointEigenSolver<ColumnMajor<S>> solver(input);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Diagonalization of Gamma failed at omega=" + std::to_string(omega) + ".");

  const auto n = gamma.rows();
  std::vector<double> rho(static_cast<std::size_t>(n));
  Matrix<S> vectors = Matrix<S>::Zero(n, n);
  // Eigen returns the eigenvalues in ascending order; the branches are labelled with the dominant one first.
  for (Eigen::Index a = 0; a < n; a++) {
    rho[static_cast<std::size_t>(a)] = solver.eigenvalues()(n - 1 - a);
    vectors.col(a)                   = solver.eigenvectors().col(n - 1 - a);
  }

  // Gamma must be positive semidefinite. Rounding in the input or in the solver can produce a slightly negative
  // eigenvalue, which is clamped; anything larger is an error, reported where it happens.
  const auto scale = std::max(rho.front(), 0.0);
  for (std::size_t a = 0; a < rho.size(); a++) {
    if (rho[a] < -options.psd_tolerance * scale)
      throw std::runtime_error("Gamma is not positive semidefinite: at omega=" + std::to_string(omega)
                               + " eigenvalue " + std::to_string(a + 1) + " is " + std::to_string(rho[a]) + ".");
    if (rho[a] < 0.0) rho[a] = 0.0;
  }
  return std::make_pair(rho, vectors);
}

// Above this many channels the assignment is made greedily rather than by trying every permutation.
constexpr int max_exact_tracking_channels = 6;

// O(i,j) = |<u_i(previous)|u_j(current)>|^2. Real and nonnegative, so the arbitrary phase of each eigenvector drops
// out.
template<typename S> auto overlap_matrix(const Matrix<S> &previous, const Matrix<S> &current) {
  const Matrix<S> products = previous.adjoint() * current;
  return Matrix<double>(products.cwiseAbs2());
}

// The assignment of current eigenvectors to previous labels that maximizes the total overlap.
inline auto best_permutation(const Matrix<double> &overlap) {
  const auto n = static_cast<std::size_t>(overlap.rows());
  std::vector<std::size_t> assignment(n);

  if (n <= static_cast<std::size_t>(max_exact_tracking_channels)) {
    std::vector<std::size_t> candidate(n);
    for (std::size_t i = 0; i < n; i++) candidate[i] = i;
    auto best = -1.0;
    do {
      double total = 0.0;
      for (std::size_t i = 0; i < n; i++)
        total += overlap(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(candidate[i]));
      if (total > best) {
        best       = total;
        assignment = candidate;
      }
    } while (std::next_permutation(candidate.begin(), candidate.end()));
    return assignment;
  }

  // Greedy fallback: repeatedly take the largest remaining overlap.
  std::vector<bool> label_taken(n, false);
  std::vector<bool> vector_taken(n, false);
  for (std::size_t step = 0; step < n; step++) {
    auto best = -1.0;
    std::size_t best_label  = 0;
    std::size_t best_vector = 0;
    for (std::size_t i = 0; i < n; i++) {
      if (label_taken[i]) continue;
      for (std::size_t j = 0; j < n; j++) {
        if (vector_taken[j]) continue;
        const auto value = overlap(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
        if (value > best) {
          best        = value;
          best_label  = i;
          best_vector = j;
        }
      }
    }
    label_taken[best_label]   = true;
    vector_taken[best_vector] = true;
    assignment[best_label]    = best_vector;
  }
  return assignment;
}

// Lowdin orthonormalization W (W^dagger W)^{-1/2}: the orthonormal set closest to the columns of W. Returns false if
// W is too close to rank deficient for the inverse square root to mean anything.
template<typename S> bool lowdin(Matrix<S> &w, const double threshold = 1e-8) {
  const ColumnMajor<S> gram = (w.adjoint() * w).eval();
  Eigen::SelfAdjointEigenSolver<ColumnMajor<S>> solver(gram);
  if (solver.info() != Eigen::Success) return false;
  const auto &eigenvalues = solver.eigenvalues();
  if (eigenvalues(0) < threshold) return false;
  const ColumnMajor<S> inverse_root = solver.eigenvectors()
                                      * eigenvalues.cwiseInverse().cwiseSqrt().template cast<S>().asDiagonal()
                                      * solver.eigenvectors().adjoint();
  w = (w * Matrix<S>(inverse_root)).eval();
  return true;
}

// Within a degenerate subspace the solver returns an arbitrary basis, which no relabelling can repair. Rotate each
// such subspace to the orthonormal basis closest to the previous node's vectors before the labels are assigned.
template<typename S>
void align_degenerate_subspaces(const std::vector<double> &rho, Matrix<S> &vectors, const Matrix<S> &previous,
                                const double tolerance) {
  const auto n     = rho.size();
  const auto scale = std::max(rho.front(), 0.0);
  if (!(scale > 0.0)) return; // Gamma vanishes here, so the labelling carries no weight either way

  for (std::size_t first = 0; first < n;) {
    std::size_t last = first; // the eigenvalues are descending, so a cluster is a contiguous run
    while (last + 1 < n && std::abs(rho[last + 1] - rho[first]) <= tolerance * scale) last++;
    const auto size = last + 1 - first;
    if (size == 1) {
      first = last + 1;
      continue;
    }

    // The projector onto the subspace, and the projections of every previous vector onto it.
    const auto block = static_cast<Eigen::Index>(size);
    const auto start = static_cast<Eigen::Index>(first);
    const Matrix<S> subspace   = vectors.block(0, start, vectors.rows(), block);
    const Matrix<S> components = subspace.adjoint() * previous; // size x n

    // Keep the previous vectors that live mostly inside this subspace, one per dimension of it.
    std::vector<std::pair<double, Eigen::Index>> weight;
    for (Eigen::Index column = 0; column < components.cols(); column++)
      weight.emplace_back(components.col(column).squaredNorm(), column);
    std::sort(weight.begin(), weight.end(), [](const auto &a, const auto &b) { return a.first > b.first; });

    Matrix<S> projected = Matrix<S>::Zero(vectors.rows(), block);
    for (Eigen::Index index = 0; index < block; index++)
      projected.col(index) = subspace * components.col(weight[static_cast<std::size_t>(index)].second);

    // If the projections are (nearly) linearly dependent they carry no usable information about the previous basis,
    // and the vectors of the solver are kept.
    if (lowdin(projected)) vectors.block(0, start, vectors.rows(), block) = projected;
    first = last + 1;
  }
}

// Relabel the eigenvalues and the eigenvectors: label i takes the vector assignment[i].
template<typename S>
void apply_assignment(const std::vector<std::size_t> &assignment, std::vector<double> &rho, Matrix<S> &vectors) {
  const auto original_rho     = rho;
  const Matrix<S> original    = vectors;
  for (std::size_t i = 0; i < assignment.size(); i++) {
    rho[i]           = original_rho[assignment[i]];
    vectors.col(static_cast<Eigen::Index>(i)) = original.col(static_cast<Eigen::Index>(assignment[i]));
  }
}

} // namespace detail

// Diagonalize Gamma at a single frequency and label the eigenvectors against a reference set, exactly as
// decompose_branch() labels one node against the previous one. This is what the star needs at a representative
// energy, which lies between the input nodes; the reference is then the tracked basis at a bracketing node.
template<typename S>
auto labelled_vectors(const Matrix<S> &gamma, const double omega, const Matrix<S> &reference,
                      const BranchOptions &options) {
  auto [rho, vectors] = detail::diagonalize<S>(gamma, omega, options);
  if (options.ordering == BranchOrdering::tracked) {
    detail::align_degenerate_subspaces(rho, vectors, reference, options.degeneracy_tolerance);
    const auto assignment = detail::best_permutation(detail::overlap_matrix(reference, vectors));
    detail::apply_assignment(assignment, rho, vectors);
  }
  return std::make_pair(rho, vectors);
}

// Decompose one frequency branch of Gamma into its eigenvalue branches, labelled consistently across the nodes.
//
// The labels are assigned going outward from omega=0. Where Gamma vanishes the labelling is arbitrary, but those
// branches carry no weight; tracking becomes meaningful at the first node that does.
//
// The tracking is performed even when the sorted ordering is requested, so that the crossings can be reported in
// both modes; only the emitted labelling differs.
template<typename S> auto decompose_branch(const GammaBranch<S> &branch, const BranchOptions &options) {
  if (branch.size() < 2) throw std::runtime_error("A frequency branch must have at least two nodes.");

  BranchDecomposition<S> result;
  result.channels = static_cast<int>(branch.gamma.front().rows());
  result.omega    = branch.omega;
  result.density.assign(static_cast<std::size_t>(result.channels), Vec{});
  for (auto &table : result.density) table.reserve(branch.size());
  result.vectors.reserve(branch.size());

  const auto channels = static_cast<std::size_t>(result.channels);
  std::vector<std::size_t> identity(channels);
  for (std::size_t i = 0; i < channels; i++) identity[i] = i;

  Matrix<S> previous_tracked;
  auto previous_assignment = identity;

  for (std::size_t k = 0; k < branch.size(); k++) {
    auto [rho, vectors] = detail::diagonalize<S>(branch.gamma[k], branch.omega[k], options);
    auto tracked_rho     = rho;
    auto tracked_vectors = vectors;

    if (k > 0) {
      detail::align_degenerate_subspaces(rho, tracked_vectors, previous_tracked, options.degeneracy_tolerance);
      const auto assignment = detail::best_permutation(detail::overlap_matrix(previous_tracked, tracked_vectors));
      detail::apply_assignment(assignment, tracked_rho, tracked_vectors);
      // A crossing is where the correspondence between the tracked labels and the sorted order changes, not every
      // node at which the two happen to differ.
      if (assignment != previous_assignment) result.crossings.push_back(branch.omega[k]);
      previous_assignment = assignment;
    }
    previous_tracked = tracked_vectors;

    const auto tracked = options.ordering == BranchOrdering::tracked;
    for (std::size_t a = 0; a < channels; a++)
      result.density[a].emplace_back(branch.omega[k], tracked ? tracked_rho[a] : rho[a]);
    result.vectors.push_back(tracked ? tracked_vectors : vectors);
  }
  return result;
}

} // namespace NRG::MixChain

#endif
