#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <mixchain/branches.hpp>

using namespace NRG::MixChain;

namespace {

// A branch built from a function of omega, so that no files are involved.
template<typename S>
GammaBranch<S> make_branch(const std::vector<double> &omega, const std::function<Matrix<S>(double)> &gamma) {
  GammaBranch<S> branch;
  branch.omega = omega;
  for (const auto value : omega) branch.gamma.push_back(gamma(value));
  return branch;
}

template<typename S> Matrix<S> diagonal(const double first, const double second) {
  Matrix<S> m = Matrix<S>::Zero(2, 2);
  m(0, 0)     = first;
  m(1, 1)     = second;
  return m;
}

// A real rotation by angle theta, used to give the eigenvectors a frequency dependence.
Matrix<double> rotation(const double theta) {
  Matrix<double> r = Matrix<double>::Zero(2, 2);
  r(0, 0)          = std::cos(theta);
  r(0, 1)          = -std::sin(theta);
  r(1, 0)          = std::sin(theta);
  r(1, 1)          = std::cos(theta);
  return r;
}

BranchOptions tracked_options() {
  BranchOptions options;
  options.ordering = BranchOrdering::tracked;
  return options;
}

BranchOptions sorted_options() {
  BranchOptions options;
  options.ordering = BranchOrdering::sorted;
  return options;
}

// |<e_i|u_a(omega_k)>|, which is independent of the arbitrary phase of the eigenvector.
template<typename S>
double component(const BranchDecomposition<S> &decomposition, const std::size_t k, const int a, const int i) {
  return std::abs(decomposition.vectors[k](i, a));
}

const std::vector<double> crossing_grid{0.1, 0.2, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.8, 0.9};

} // namespace

TEST(MixChainBranches, orders_branches_by_descending_eigenvalue) { // NOLINT
  const auto branch = make_branch<double>({0.2, 0.4}, [](const double) { return diagonal<double>(0.3, 0.8); });
  const auto result = decompose_branch(branch, tracked_options());

  ASSERT_EQ(result.channels, 2);
  EXPECT_DOUBLE_EQ(result.density[0][0].second, 0.8); // the dominant branch is labelled first
  EXPECT_DOUBLE_EQ(result.density[1][0].second, 0.3);
  EXPECT_NEAR(component(result, 0, 0, 1), 1.0, 1e-15); // and it is the second channel here
  EXPECT_NEAR(component(result, 0, 1, 0), 1.0, 1e-15);
}

TEST(MixChainBranches, without_a_crossing_the_two_orderings_agree) { // NOLINT
  const auto branch = make_branch<double>(crossing_grid, [](const double omega) {
    return diagonal<double>(1.0 - 0.5 * omega, 0.1 * omega); // no crossing on this grid
  });
  const auto tracked = decompose_branch(branch, tracked_options());
  const auto sorted  = decompose_branch(branch, sorted_options());

  for (std::size_t k = 0; k < branch.size(); k++) {
    EXPECT_DOUBLE_EQ(tracked.density[0][k].second, sorted.density[0][k].second);
    EXPECT_DOUBLE_EQ(tracked.density[1][k].second, sorted.density[1][k].second);
    EXPECT_DOUBLE_EQ(tracked.density[0][k].second, 1.0 - 0.5 * branch.omega[k]);
  }
  EXPECT_TRUE(tracked.crossings.empty());
  EXPECT_TRUE(sorted.crossings.empty());
}

TEST(MixChainBranches, tracking_follows_a_branch_through_a_crossing) { // NOLINT
  // The two eigenvalues cross at omega=0.5, where the true chain is still decoupled.
  const auto branch = make_branch<double>(crossing_grid, [](const double omega) {
    return diagonal<double>(1.0 - omega, omega);
  });

  const auto tracked = decompose_branch(branch, tracked_options());
  for (std::size_t k = 0; k < branch.size(); k++) {
    const auto omega = branch.omega[k];
    EXPECT_DOUBLE_EQ(tracked.density[0][k].second, 1.0 - omega); // the label follows the channel
    EXPECT_DOUBLE_EQ(tracked.density[1][k].second, omega);
    EXPECT_NEAR(component(tracked, k, 0, 0), 1.0, 1e-14);
    EXPECT_NEAR(component(tracked, k, 1, 1), 1.0, 1e-14);
  }

  const auto sorted = decompose_branch(branch, sorted_options());
  for (std::size_t k = 0; k < branch.size(); k++) {
    const auto omega = branch.omega[k];
    // The sorted labelling jumps from one channel to the other at the crossing.
    EXPECT_DOUBLE_EQ(sorted.density[0][k].second, std::max(1.0 - omega, omega));
    EXPECT_DOUBLE_EQ(sorted.density[1][k].second, std::min(1.0 - omega, omega));
  }

  // Both modes report the crossing, at the first node past it, and only there.
  ASSERT_EQ(tracked.crossings.size(), 1U);
  EXPECT_DOUBLE_EQ(tracked.crossings.front(), 0.55);
  EXPECT_EQ(sorted.crossings, tracked.crossings);
}

TEST(MixChainBranches, tracking_follows_rotating_eigenvectors) { // NOLINT
  // An avoided crossing: the eigenvalues cross while the eigenvectors rotate slowly.
  // The grid deliberately misses omega=0.5, where the eigenvalues are exactly degenerate and the eigenvectors are
  // therefore undefined; that case is covered on its own below.
  std::vector<double> grid;
  for (double omega = 0.3025; omega < 0.70; omega += 0.005) grid.push_back(omega);
  const auto branch = make_branch<double>(grid, [](const double omega) {
    const auto r = rotation(0.5 * M_PI * omega);
    return Matrix<double>(r * diagonal<double>(1.0 - omega, omega) * r.transpose());
  });

  const auto tracked = decompose_branch(branch, tracked_options());
  for (std::size_t k = 0; k < branch.size(); k++) {
    const auto omega = branch.omega[k];
    EXPECT_NEAR(tracked.density[0][k].second, 1.0 - omega, 1e-12);
    EXPECT_NEAR(tracked.density[1][k].second, omega, 1e-12);
    // The first branch keeps the rotating eigenvector, so its overlap with the initial one decreases smoothly.
    const auto r = rotation(0.5 * M_PI * omega);
    EXPECT_NEAR(std::abs(tracked.vectors[k].col(0).dot(r.col(0))), 1.0, 1e-10);
  }
  EXPECT_EQ(tracked.crossings.size(), 1U);
}

TEST(MixChainBranches, a_degenerate_gamma_keeps_an_orthonormal_basis) { // NOLINT
  // Gamma = rho(omega) * identity: the eigenvectors are arbitrary, and the solver may return a different basis at
  // every node unless the degenerate subspace is aligned with the previous one.
  const auto branch = make_branch<double>({0.1, 0.3, 0.5, 0.7}, [](const double omega) {
    return diagonal<double>(0.5 * omega, 0.5 * omega);
  });
  const auto result = decompose_branch(branch, tracked_options());

  const Matrix<double> identity = Matrix<double>::Identity(2, 2);
  for (std::size_t k = 0; k < branch.size(); k++) {
    EXPECT_DOUBLE_EQ(result.density[0][k].second, 0.5 * branch.omega[k]);
    EXPECT_DOUBLE_EQ(result.density[1][k].second, 0.5 * branch.omega[k]);
    const Matrix<double> product = result.vectors[k].transpose() * result.vectors[k];
    EXPECT_LT((product - identity).cwiseAbs().maxCoeff(), 1e-14); // orthonormal
    // The basis does not wander from node to node.
    const Matrix<double> overlap = result.vectors[k].transpose() * result.vectors[0];
    EXPECT_NEAR(std::abs(overlap(0, 0)), 1.0, 1e-12);
    EXPECT_NEAR(std::abs(overlap(1, 1)), 1.0, 1e-12);
  }
  EXPECT_TRUE(result.crossings.empty());
}

TEST(MixChainBranches, rejects_a_gamma_that_is_not_positive_semidefinite) { // NOLINT
  const auto branch = make_branch<double>({0.2, 0.4}, [](const double) { return diagonal<double>(0.3, -0.5); });
  EXPECT_THROW(decompose_branch(branch, tracked_options()), std::runtime_error);
}

TEST(MixChainBranches, clamps_a_negative_eigenvalue_within_the_tolerance) { // NOLINT
  const auto branch = make_branch<double>({0.2, 0.4}, [](const double) { return diagonal<double>(0.3, -1e-12); });
  const auto result = decompose_branch(branch, tracked_options());
  EXPECT_DOUBLE_EQ(result.density[0][0].second, 0.3);
  EXPECT_DOUBLE_EQ(result.density[1][0].second, 0.0); // clamped, not merely small
}

TEST(MixChainBranches, decomposes_a_complex_hermitian_gamma) { // NOLINT
  using Complex = std::complex<double>;
  const double a = 0.5;
  const double b = 0.2;
  const auto branch = make_branch<Complex>({0.2, 0.4}, [a, b](const double) {
    Matrix<Complex> m = Matrix<Complex>::Zero(2, 2);
    m(0, 0)           = a;
    m(1, 1)           = a;
    m(0, 1)           = Complex(0.0, b);
    m(1, 0)           = Complex(0.0, -b);
    return m;
  });
  const auto result = decompose_branch(branch, tracked_options());

  EXPECT_NEAR(result.density[0][0].second, a + b, 1e-14);
  EXPECT_NEAR(result.density[1][0].second, a - b, 1e-14);
  const Matrix<Complex> identity = Matrix<Complex>::Identity(2, 2);
  const Matrix<Complex> product  = result.vectors[0].adjoint() * result.vectors[0];
  EXPECT_LT((product - identity).cwiseAbs().maxCoeff(), 1e-14);
}

TEST(MixChainBranches, parses_the_ordering_option) { // NOLINT
  EXPECT_EQ(branch_ordering_from_string("tracked"), BranchOrdering::tracked);
  EXPECT_EQ(branch_ordering_from_string("sorted"), BranchOrdering::sorted);
  EXPECT_EQ(branch_ordering_name(BranchOrdering::sorted), "sorted");
  EXPECT_THROW(branch_ordering_from_string("overlap"), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
