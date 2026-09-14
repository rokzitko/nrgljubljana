#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <vector>

#include <mixchain/star.hpp>

using namespace NRG::MixChain;

namespace {

using Complex = std::complex<double>;

constexpr double lambda_value = 2.0;
constexpr unsigned int mmax   = 10;

// The input grid: dense enough that the tracking has something to follow, and reaching the band edge.
std::vector<double> input_grid() {
  std::vector<double> grid;
  for (int k = 0; k <= 100; k++) grid.push_back(0.01 * k);
  return grid;
}

// Gamma(omega) is taken to be the same on both frequency branches, which are indexed by |omega|.
template<typename S> GammaInput<S> make_input(const std::function<Matrix<S>(double)> &gamma) {
  const auto grid = input_grid();
  GammaBranch<S> branch;
  branch.omega = grid;
  for (const auto omega : grid) branch.gamma.push_back(gamma(omega));
  GammaInput<S> input;
  input.channels = static_cast<int>(branch.gamma.front().rows());
  input.pos      = branch;
  input.neg      = branch;
  return input;
}

template<typename S> Matrix<S> scalar(const double value) {
  Matrix<S> m = Matrix<S>::Zero(1, 1);
  m(0, 0)     = value;
  return m;
}

template<typename S> Matrix<S> diagonal(const double first, const double second) {
  Matrix<S> m = Matrix<S>::Zero(2, 2);
  m(0, 0)     = first;
  m(1, 1)     = second;
  return m;
}

StarOptions base_options() {
  StarOptions options;
  options.Lambda = NRG::Tools::LambdaCache(lambda_value);
  options.z      = 1.0;
  options.mMAX   = mmax;
  return options;
}

// The levels of one frequency branch and one eigenvalue branch, in order of increasing m.
template<typename S>
std::vector<StarLevel<S>> select(const Star<S> &star, const Sign sign, const int branch) {
  std::vector<StarLevel<S>> levels;
  for (const auto &level : star.levels)
    if (level.sign == sign && level.branch == branch) levels.push_back(level);
  return levels;
}

} // namespace

TEST(MixChainStar, flat_band_reproduces_the_analytic_star) { // NOLINT
  const double rho = 0.3;
  const auto input = make_input<double>([rho](const double) { return scalar<double>(rho); });
  const auto star  = build_star(input, base_options());

  ASSERT_EQ(star.levels.size(), 2 * (mmax + 1));
  const auto factor = NRG::Tools::LambdaCache(lambda_value).factor(); // (1-1/Lambda)/ln(Lambda)

  for (const auto sign : {Sign::POS, Sign::NEG}) {
    const auto levels = select(star, sign, 0);
    ASSERT_EQ(levels.size(), mmax + 1);
    for (unsigned int m = 0; m <= mmax; m++) {
      // With z=1 the interval of index m is [Lambda^(-m-1), Lambda^(-m)].
      const auto power    = std::pow(lambda_value, -static_cast<double>(m));
      const auto weight   = rho * power * (1.0 - 1.0 / lambda_value);
      const auto energy   = factor * power;
      const auto coupling = levels[m].coupling(0);
      EXPECT_NEAR(std::abs(levels[m].energy), energy, 1e-13 * energy);
      EXPECT_DOUBLE_EQ(levels[m].energy, sign_value(sign) * std::abs(levels[m].energy));
      EXPECT_NEAR(coupling * coupling, weight, 1e-13 * weight);
    }
  }

  // Theta is the weight of the star; the mesh covers exactly [Lambda^(-mMAX-1), 1] on each branch.
  const auto covered = rho * (1.0 - std::pow(lambda_value, -static_cast<double>(mmax) - 1.0));
  EXPECT_NEAR(star.theta(0, 0), 2.0 * covered, 1e-13);
  EXPECT_NEAR(star.theta(0, 0), star.theta_exact(0, 0), 1e-13);
  EXPECT_LT(star.diagnostics.max_interval_deviation, 1e-13);
}

TEST(MixChainStar, is_covariant_under_a_constant_rotation) { // NOLINT
  const auto densities = [](const double omega) { return diagonal<double>(0.5 + 0.1 * omega, 0.2 + 0.05 * omega); };
  const auto plain     = build_star(make_input<double>(densities), base_options());

  Matrix<double> u = Matrix<double>::Zero(2, 2);
  const auto angle = 0.7;
  u(0, 0) = std::cos(angle);
  u(0, 1) = -std::sin(angle);
  u(1, 0) = std::sin(angle);
  u(1, 1) = std::cos(angle);
  const auto rotated = build_star(
    make_input<double>([&](const double omega) { return Matrix<double>(u * densities(omega) * u.transpose()); }),
    base_options());

  ASSERT_EQ(plain.levels.size(), rotated.levels.size());
  for (std::size_t k = 0; k < plain.levels.size(); k++) {
    EXPECT_NEAR(rotated.levels[k].energy, plain.levels[k].energy, 1e-12);
    // The overall phase of an eigenvector is arbitrary, so the projectors are compared.
    const Matrix<double> projector = plain.levels[k].coupling * plain.levels[k].coupling.transpose();
    const Matrix<double> expected  = u * projector * u.transpose();
    const Matrix<double> found = rotated.levels[k].coupling * rotated.levels[k].coupling.transpose();
    EXPECT_LT((found - expected).cwiseAbs().maxCoeff(), 1e-12);
  }
  const Matrix<double> expected_theta = u * plain.theta * u.transpose();
  EXPECT_LT((rotated.theta - expected_theta).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(MixChainStar, is_covariant_under_a_complex_rotation) { // NOLINT
  const auto densities = [](const double omega) { return diagonal<Complex>(0.5 + 0.1 * omega, 0.2 + 0.05 * omega); };
  const auto plain     = build_star(make_input<Complex>(densities), base_options());

  Matrix<Complex> u = Matrix<Complex>::Zero(2, 2);
  u(0, 0)           = Complex(0.6, 0.0);
  u(0, 1)           = Complex(0.0, 0.8);
  u(1, 0)           = Complex(0.0, 0.8);
  u(1, 1)           = Complex(0.6, 0.0);
  const auto rotated = build_star(
    make_input<Complex>([&](const double omega) { return Matrix<Complex>(u * densities(omega) * u.adjoint()); }),
    base_options());

  for (std::size_t k = 0; k < plain.levels.size(); k++) {
    EXPECT_NEAR(rotated.levels[k].energy, plain.levels[k].energy, 1e-12);
    const Matrix<Complex> projector = plain.levels[k].coupling * plain.levels[k].coupling.adjoint();
    const Matrix<Complex> expected  = u * projector * u.adjoint();
    const Matrix<Complex> found = rotated.levels[k].coupling * rotated.levels[k].coupling.adjoint();
    EXPECT_LT((found - expected).cwiseAbs().maxCoeff(), 1e-12);
  }
}

TEST(MixChainStar, scales_the_couplings_with_the_square_root_of_gamma) { // NOLINT
  const double factor = 7.5;
  const auto plain    = build_star(make_input<double>([](const double omega) { return scalar<double>(0.4 + omega); }),
                                   base_options());
  const auto scaled   = build_star(
    make_input<double>([factor](const double omega) { return scalar<double>(factor * (0.4 + omega)); }),
    base_options());

  for (std::size_t k = 0; k < plain.levels.size(); k++) {
    EXPECT_NEAR(scaled.levels[k].energy, plain.levels[k].energy, 1e-12);
    EXPECT_NEAR(scaled.levels[k].coupling(0), std::sqrt(factor) * plain.levels[k].coupling(0), 1e-12);
  }
  EXPECT_NEAR(scaled.theta(0, 0), factor * plain.theta(0, 0), 1e-12);
}

TEST(MixChainStar, degenerate_branches_give_orthogonal_couplings) { // NOLINT
  // Gamma = rho(omega) * identity: both branches have the same density, hence the same representative energies.
  const auto input = make_input<double>([](const double omega) {
    const auto rho = 0.4 + 0.2 * omega;
    return diagonal<double>(rho, rho);
  });
  const auto star = build_star(input, base_options());

  const auto first  = select(star, Sign::POS, 0);
  const auto second = select(star, Sign::POS, 1);
  ASSERT_EQ(first.size(), mmax + 1);
  for (unsigned int m = 0; m <= mmax; m++) {
    EXPECT_NEAR(first[m].energy, second[m].energy, 1e-14);
    EXPECT_NEAR(first[m].coupling.dot(second[m].coupling), 0.0, 1e-14); // one diagonalization, orthonormal vectors
    const Matrix<double> sum = first[m].coupling * first[m].coupling.transpose()
                               + second[m].coupling * second[m].coupling.transpose();
    const Matrix<double> identity = Matrix<double>::Identity(2, 2);
    EXPECT_LT((sum - first[m].coupling.squaredNorm() * identity).cwiseAbs().maxCoeff(), 1e-14);
  }
  EXPECT_LT(star.diagnostics.max_interval_deviation, 1e-12);
}

TEST(MixChainStar, a_diagonal_gamma_matches_independent_scalar_runs) { // NOLINT
  // The densities do not cross, so branch 0 is the first channel at every node: the branches are labelled by
  // descending eigenvalue at the innermost node and tracked from there.
  const auto first_density  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second_density = [](const double omega) { return 0.3 + 0.1 * omega; };

  const auto joint = build_star(make_input<double>([&](const double omega) {
                                  return diagonal<double>(first_density(omega), second_density(omega));
                                }),
                                base_options());

  const auto compare = [&joint](const std::function<double(double)> &density, const int branch) {
    const auto alone = build_star(
      make_input<double>([&](const double omega) { return scalar<double>(density(omega)); }), base_options());
    const auto joint_levels = select(joint, Sign::POS, branch);
    const auto alone_levels = select(alone, Sign::POS, 0);
    ASSERT_EQ(joint_levels.size(), alone_levels.size());
    for (std::size_t m = 0; m < joint_levels.size(); m++) {
      EXPECT_NEAR(joint_levels[m].energy, alone_levels[m].energy, 1e-13);
      EXPECT_NEAR(joint_levels[m].coupling.norm(), std::abs(alone_levels[m].coupling(0)), 1e-13);
      // The coupling points along this branch's own channel.
      EXPECT_NEAR(std::abs(joint_levels[m].coupling(branch)), joint_levels[m].coupling.norm(), 1e-13);
      EXPECT_NEAR(joint_levels[m].coupling(1 - branch), 0.0, 1e-14);
    }
  };
  compare(first_density, 0);
  compare(second_density, 1);

  EXPECT_NEAR(joint.theta(0, 1), 0.0, 1e-14); // the channels stay decoupled
  EXPECT_TRUE(joint.diagnostics.crossings_pos.empty());
}

TEST(MixChainStar, tracking_preserves_the_interval_sum_rule_at_a_crossing) { // NOLINT
  // The two densities cross at omega=0.4, inside the interval [0.25,0.5] of the Lambda=2, z=1 mesh.
  const auto input = make_input<double>(
    [](const double omega) { return diagonal<double>(0.5 + 0.25 * omega, 0.8 - 0.5 * omega); });

  auto options    = base_options();
  const auto tracked = build_star(input, options);
  options.branches.ordering = BranchOrdering::sorted;
  const auto sorted         = build_star(input, options);

  // With the branches tracked, the star reproduces the integral of Gamma over every interval exactly.
  EXPECT_LT(tracked.diagnostics.max_interval_deviation, 1e-13);
  // With the branches sorted, the interval containing the crossing receives the weight of the wrong channel. The
  // deviation there is of order a few percent for these densities.
  EXPECT_GT(sorted.diagnostics.max_interval_deviation, 1e-3);
  EXPECT_NEAR(sorted.diagnostics.max_interval_omega, 0.5, 1e-12);
  EXPECT_FALSE(tracked.diagnostics.crossings_pos.empty());
}

TEST(MixChainStar, an_empty_branch_gives_vanishing_couplings) { // NOLINT
  // Gamma = diag(rho, 0): the second channel does not hybridize at all.
  const auto input = make_input<double>([](const double omega) { return diagonal<double>(0.4 + omega, 0.0); });
  const auto star  = build_star(input, base_options());

  ASSERT_EQ(star.levels.size(), 2 * 2 * (mmax + 1));
  for (const auto &level : select(star, Sign::POS, 1)) {
    EXPECT_GT(std::abs(level.energy), 0.0);
    EXPECT_TRUE(std::isfinite(level.energy));
    EXPECT_DOUBLE_EQ(level.coupling.norm(), 0.0);
  }
  EXPECT_GT(star.theta(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(star.theta(1, 1), 0.0); // Theta is singular, which only the chain stage has to care about
}

TEST(MixChainStar, reports_a_mesh_that_reaches_below_the_input) { // NOLINT
  // The input is tabulated down to 0.01 while the mesh of mMAX=10 reaches Lambda^(-11); below the innermost node the
  // density is the constant continuation added at omega=0.
  auto input = make_input<double>([](const double) { return scalar<double>(0.3); });
  input.pos.innermost = 0.01;
  input.neg.innermost = 0.01;
  const auto star     = build_star(input, base_options());

  EXPECT_TRUE(star.diagnostics.coverage_pos.continued());
  EXPECT_TRUE(star.diagnostics.coverage_neg.continued());
  EXPECT_DOUBLE_EQ(star.diagnostics.coverage_pos.innermost_input, 0.01);
  EXPECT_NEAR(star.diagnostics.coverage_pos.lowest_mesh, std::pow(lambda_value, -11.0), 1e-15);

  // An input that is tabulated below everything the mesh reaches says nothing.
  input.pos.innermost = 1e-30;
  input.neg.innermost = 1e-30;
  EXPECT_FALSE(build_star(input, base_options()).diagnostics.coverage_pos.continued());
}

TEST(MixChainStar, reports_intervals_that_hold_no_node_of_the_input) { // NOLINT
  // The input grid of make_input() runs from 0 to 1 in steps of 0.01, so every interval below 0.01 falls inside a
  // single tabulated interval: for Lambda=2 and z=1 those are m >= 7, since Lambda^(-7) = 0.0078.
  const auto star = build_star(make_input<double>([](const double) { return scalar<double>(0.3); }), base_options());
  const auto &coverage = star.diagnostics.coverage_pos;
  EXPECT_EQ(coverage.unresolved_intervals, static_cast<int>(mmax) + 1 - 7);
  EXPECT_NEAR(coverage.unresolved_from, std::pow(lambda_value, -7.0), 1e-15);

  // A grid that resolves every interval reports nothing.
  GammaInput<double> fine;
  fine.channels = 1;
  fine.pos.omega.push_back(0.0);
  fine.pos.gamma.push_back(scalar<double>(0.3));
  for (int k = -40; k <= 0; k++) { // two nodes per interval of the mesh
    for (const double factor : {1.0, 1.5}) {
      fine.pos.omega.push_back(factor * std::pow(lambda_value, k));
      fine.pos.gamma.push_back(scalar<double>(0.3));
    }
  }
  fine.neg = fine.pos;
  EXPECT_EQ(build_star(fine, base_options()).diagnostics.coverage_pos.unresolved_intervals, 0);
}

TEST(MixChainStar, rejects_invalid_options) { // NOLINT
  const auto input = make_input<double>([](const double) { return scalar<double>(0.3); });

  auto options   = base_options();
  options.Lambda = NRG::Tools::LambdaCache(1.0);
  EXPECT_THROW(build_star(input, options), std::invalid_argument);

  options   = base_options();
  options.z = 0.0;
  EXPECT_THROW(build_star(input, options), std::invalid_argument);
  options.z = 1.5;
  EXPECT_THROW(build_star(input, options), std::invalid_argument);

  options      = base_options();
  options.mMAX = 0;
  EXPECT_THROW(build_star(input, options), std::invalid_argument);

  options               = base_options();
  options.allowed_error = 0.0;
  EXPECT_THROW(build_star(input, options), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
