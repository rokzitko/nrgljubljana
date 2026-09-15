#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <utility>
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
  EXPECT_LT(star.diagnostics[0].max_interval_deviation, 1e-13);
}

TEST(MixChainStar, is_covariant_under_a_constant_rotation) { // NOLINT
  const auto densities = [](const double omega) { return diagonal<double>(0.5 + 0.1 * omega, 0.2 + 0.05 * omega); };
  // The rotated Gamma is a single block; the plain one is kept whole too, so that both go through the same path.
  auto options         = base_options();
  options.split_blocks = false;
  const auto plain     = build_star(make_input<double>(densities), options);

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
  auto options         = base_options();
  options.split_blocks = false;
  const auto plain     = build_star(make_input<Complex>(densities), options);

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
  auto options         = base_options();
  options.split_blocks = false; // the whole matrix, whose eigenvalues are degenerate
  const auto star      = build_star(input, options);

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
  EXPECT_LT(star.diagnostics[0].max_interval_deviation, 1e-12);
}

TEST(MixChainStar, a_diagonal_gamma_matches_independent_scalar_runs) { // NOLINT
  // The densities do not cross, so branch 0 is the first channel at every node: the branches are labelled by
  // descending eigenvalue at the innermost node and tracked from there.
  const auto first_density  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second_density = [](const double omega) { return 0.3 + 0.1 * omega; };

  auto options         = base_options();
  options.split_blocks = false; // the whole matrix; split_blocks is covered by the tests below
  const auto joint     = build_star(make_input<double>([&](const double omega) {
                                  return diagonal<double>(first_density(omega), second_density(omega));
                                }),
                                options);

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
  EXPECT_TRUE(joint.diagnostics[0].crossings_pos.empty());
}

TEST(MixChainStar, tracking_preserves_the_interval_sum_rule_at_a_crossing) { // NOLINT
  // The two densities cross at omega=0.4, inside the interval [0.25,0.5] of the Lambda=2, z=1 mesh.
  const auto input = make_input<double>(
    [](const double omega) { return diagonal<double>(0.5 + 0.25 * omega, 0.8 - 0.5 * omega); });

  auto options         = base_options();
  options.split_blocks = false; // split, each channel is its own branch and there is nothing to order
  const auto tracked   = build_star(input, options);
  options.branches.ordering = BranchOrdering::sorted;
  const auto sorted         = build_star(input, options);

  // With the branches tracked, the star reproduces the integral of Gamma over every interval exactly.
  EXPECT_LT(tracked.diagnostics[0].max_interval_deviation, 1e-13);
  // With the branches sorted, the interval containing the crossing receives the weight of the wrong channel. The
  // deviation there is of order a few percent for these densities.
  EXPECT_GT(sorted.diagnostics[0].max_interval_deviation, 1e-3);
  EXPECT_NEAR(sorted.diagnostics[0].max_interval_omega, 0.5, 1e-12);
  EXPECT_FALSE(tracked.diagnostics[0].crossings_pos.empty());
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

  EXPECT_TRUE(star.diagnostics[0].coverage_pos.continued());
  EXPECT_TRUE(star.diagnostics[0].coverage_neg.continued());
  EXPECT_DOUBLE_EQ(star.diagnostics[0].coverage_pos.innermost_input, 0.01);
  EXPECT_NEAR(star.diagnostics[0].coverage_pos.lowest_mesh, std::pow(lambda_value, -11.0), 1e-15);

  // An input that is tabulated below everything the mesh reaches says nothing.
  input.pos.innermost = 1e-30;
  input.neg.innermost = 1e-30;
  EXPECT_FALSE(build_star(input, base_options()).diagnostics[0].coverage_pos.continued());
}

TEST(MixChainStar, reports_intervals_that_hold_no_node_of_the_input) { // NOLINT
  // The input grid of make_input() runs from 0 to 1 in steps of 0.01, so every interval below 0.01 falls inside a
  // single tabulated interval: for Lambda=2 and z=1 those are m >= 7, since Lambda^(-7) = 0.0078.
  const auto star = build_star(make_input<double>([](const double) { return scalar<double>(0.3); }), base_options());
  const auto &coverage = star.diagnostics[0].coverage_pos;
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
  EXPECT_EQ(build_star(fine, base_options()).diagnostics[0].coverage_pos.unresolved_intervals, 0);
}

TEST(MixChainStar, counts_levels_that_collapse_onto_the_accumulation_point) { // NOLINT
  // With hardgap the levels approach boundary=0.5 as 0.5 (1-1/Lambda)/ln(Lambda) Lambda^(-m); doubles near 0.5 are
  // 1.1e-16 apart, so beyond m of about 52 the energies can no longer be told apart from 0.5.
  const auto input = make_input<double>([](const double) { return scalar<double>(0.3); });
  auto options     = base_options();
  options.mMAX     = 80;
  options.hardgap  = true;
  options.boundary = 0.5;
  const auto star  = build_star(input, options);

  const auto &coverage = star.diagnostics[0].coverage_pos;
  EXPECT_DOUBLE_EQ(coverage.accumulation_point, 0.5);
  EXPECT_GT(coverage.collapsed_levels, 20);
  EXPECT_LT(coverage.collapsed_levels, 35);

  // The collapsed levels are exactly the inert ones, with no weight, and they sit on the accumulation point or on a
  // neighbouring double: here the density does not vanish below it, so they land just above rather than on it.
  int inert = 0;
  for (const auto &level : select(star, Sign::POS, 0)) {
    if (level.coupling.norm() != 0.0) continue;
    inert++;
    EXPECT_LE(std::abs(level.energy - 0.5), 4.0 * std::numeric_limits<double>::epsilon() * 0.5);
  }
  EXPECT_EQ(inert, coverage.collapsed_levels); // one channel, so one level per collapsed interval

  // Accumulating at zero, the energies keep their relative precision and nothing collapses.
  options.hardgap = false;
  EXPECT_EQ(build_star(input, options).diagnostics[0].coverage_pos.collapsed_levels, 0);
}

TEST(MixChainStar, one_setup_serves_every_z) { // NOLINT
  // The setup that does not depend on z is shared across evaluations. Every star must be exactly the star a fresh
  // setup gives for that z: the arithmetic is the same and in the same order, and what is shared (the caches of the
  // densities and the integration workspace) affects only the speed. The densities cross, so that the branch
  // tracking and every diagnostic are exercised.
  const auto input = make_input<double>(
    [](const double omega) { return diagonal<double>(0.5 + 0.25 * omega, 0.8 - 0.5 * omega); });
  auto options         = base_options();
  options.split_blocks = false;
  StarDiscretizer<double> shared(input, options);

  for (const double z : {0.25, 0.5, 0.75, 1.0}) {
    const auto from_shared = shared.star(z);
    options.z              = z;
    const auto fresh       = build_star(input, options);

    EXPECT_EQ(from_shared.z, z);
    ASSERT_EQ(from_shared.levels.size(), fresh.levels.size());
    for (std::size_t k = 0; k < fresh.levels.size(); k++) {
      const auto &a = from_shared.levels[k];
      const auto &b = fresh.levels[k];
      EXPECT_EQ(a.m, b.m);
      EXPECT_EQ(a.branch, b.branch);
      EXPECT_TRUE(a.sign == b.sign);
      EXPECT_EQ(a.energy, b.energy) << "z=" << z << " level " << k;
      EXPECT_TRUE((a.coupling.array() == b.coupling.array()).all()) << "z=" << z << " level " << k;
    }
    EXPECT_TRUE((from_shared.theta.array() == fresh.theta.array()).all());
    EXPECT_TRUE((from_shared.theta_exact.array() == fresh.theta_exact.array()).all());
    EXPECT_EQ(from_shared.diagnostics[0].max_interval_deviation, fresh.diagnostics[0].max_interval_deviation);
    EXPECT_EQ(from_shared.diagnostics[0].max_cquad_error, fresh.diagnostics[0].max_cquad_error);
    EXPECT_EQ(from_shared.diagnostics[0].crossings_pos, fresh.diagnostics[0].crossings_pos);
    EXPECT_EQ(from_shared.diagnostics[0].coverage_pos.unresolved_intervals,
              fresh.diagnostics[0].coverage_pos.unresolved_intervals);
  }
  EXPECT_THROW(shared.star(0.0), std::invalid_argument);
  EXPECT_THROW(shared.star(1.5), std::invalid_argument);
}

TEST(MixChainStar, a_split_diagonal_gamma_is_exactly_the_scalar_problem_per_channel) { // NOLINT
  // Densities with different structure, so that the adaptive meshes of the two channels differ.
  const auto first_density  = [](const double omega) { return 0.8 * omega * omega + 0.01; };
  const auto second_density = [](const double omega) { return omega < 0.3 ? 0.0 : 0.5; };
  const auto input          = make_input<double>(
    [&](const double omega) { return diagonal<double>(first_density(omega), second_density(omega)); });

  for (const bool adapt : {false, true}) {
    auto options  = base_options();
    options.adapt = adapt;
    const auto joint = build_star(input, options);
    EXPECT_EQ(joint.blocks, (Blocks{{0}, {1}})) << "adapt=" << adapt;
    ASSERT_EQ(joint.diagnostics.size(), 2U);

    for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first_density), 0},
                                           std::pair{std::function<double(double)>(second_density), 1}}) {
      const auto alone =
        build_star(make_input<double>([&](const double omega) { return scalar<double>(density(omega)); }), options);
      for (const auto sign : {Sign::POS, Sign::NEG}) {
        const auto joint_levels = select(joint, sign, channel);
        const auto alone_levels = select(alone, sign, 0);
        ASSERT_EQ(joint_levels.size(), alone_levels.size());
        for (std::size_t m = 0; m < joint_levels.size(); m++) {
          EXPECT_EQ(joint_levels[m].energy, alone_levels[m].energy) << "adapt=" << adapt << " m=" << m;
          EXPECT_EQ(joint_levels[m].coupling(channel), alone_levels[m].coupling(0)) << "adapt=" << adapt << " m=" << m;
          EXPECT_EQ(joint_levels[m].coupling(1 - channel), 0.0);
        }
      }
      EXPECT_EQ(joint.theta(channel, channel), alone.theta(0, 0));
    }
    EXPECT_EQ(joint.theta(0, 1), 0.0);
  }

  // Kept whole, the adaptive mesh is shared by both channels, and the energies are not those of the scalar runs.
  auto options         = base_options();
  options.adapt        = true;
  options.split_blocks = false;
  const auto whole     = build_star(input, options);
  const auto alone =
    build_star(make_input<double>([&](const double omega) { return scalar<double>(first_density(omega)); }), options);
  EXPECT_EQ(whole.blocks.size(), 1U);
  EXPECT_EQ(whole.diagnostics.size(), 1U);
  EXPECT_GT(std::abs(select(whole, Sign::POS, 0)[3].energy - select(alone, Sign::POS, 0)[3].energy), 1e-6);
}

TEST(MixChainStar, blocks_are_embedded_in_their_channels) { // NOLINT
  // Channels 1 and 3 are coupled, channel 2 is on its own.
  const auto gamma = [](const double omega) {
    Matrix<Complex> m = Matrix<Complex>::Zero(3, 3);
    m(0, 0)           = 0.5 + 0.2 * omega;
    m(2, 2)           = 0.4 - 0.1 * omega;
    m(0, 2)           = Complex(0.1 * omega, 0.05);
    m(2, 0)           = std::conj(m(0, 2));
    m(1, 1)           = 0.3 + omega * omega;
    return m;
  };
  const auto input   = make_input<Complex>(gamma);
  auto options       = base_options();
  options.adapt      = true;
  const auto star    = build_star(input, options);
  ASSERT_EQ(star.blocks, (Blocks{{0, 2}, {1}}));
  ASSERT_EQ(star.levels.size(), 2 * 3 * (mmax + 1));

  // The same blocks discretized on their own, and placed in their channels by hand.
  const auto outer = build_star(restrict_input(input, Block{0, 2}), options);
  const auto inner = build_star(restrict_input(input, Block{1}), options);
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    for (int branch = 0; branch < 3; branch++) {
      const auto levels = select(star, sign, branch);
      const auto part   = branch < 2 ? select(outer, sign, branch) : select(inner, sign, 0);
      ASSERT_EQ(levels.size(), part.size());
      for (std::size_t m = 0; m < levels.size(); m++) {
        EXPECT_EQ(levels[m].m, static_cast<int>(m));
        EXPECT_EQ(levels[m].energy, part[m].energy);
        if (branch < 2) {
          EXPECT_EQ(levels[m].coupling(0), part[m].coupling(0));
          EXPECT_EQ(levels[m].coupling(1), Complex(0.0));
          EXPECT_EQ(levels[m].coupling(2), part[m].coupling(1));
        } else {
          EXPECT_EQ(levels[m].coupling(0), Complex(0.0));
          EXPECT_EQ(levels[m].coupling(1), part[m].coupling(0));
          EXPECT_EQ(levels[m].coupling(2), Complex(0.0));
        }
      }
    }
  }
  // The levels come in the order of a single block: by sign, then interval, then branch.
  const auto half = star.levels.size() / 2;
  for (std::size_t k = 0; k < star.levels.size(); k++) {
    EXPECT_TRUE(star.levels[k].sign == (k < half ? Sign::POS : Sign::NEG));
    EXPECT_EQ(star.levels[k].m, static_cast<int>((k % half) / 3));
    EXPECT_EQ(star.levels[k].branch, static_cast<int>((k % half) % 3));
  }
  EXPECT_EQ(star.theta(0, 2), outer.theta(0, 1));
  EXPECT_EQ(star.theta(1, 1), inner.theta(0, 0));
  EXPECT_EQ(star.theta(0, 1), Complex(0.0));
  EXPECT_EQ(star.theta_exact(2, 0), outer.theta_exact(1, 0));
  EXPECT_EQ(star.diagnostics.size(), 2U);
}

TEST(MixChainStar, adapt_falls_back_to_the_fixed_mesh_where_gamma_vanishes) { // NOLINT
  // Split, the second channel is a block with no weight at all.
  auto options  = base_options();
  options.adapt = true;
  const auto star =
    build_star(make_input<double>([](const double omega) { return diagonal<double>(0.4 + omega, 0.0); }), options);
  ASSERT_EQ(star.diagnostics.size(), 2U);
  EXPECT_FALSE(star.diagnostics[0].coverage_pos.fixed_mesh_fallback);
  EXPECT_TRUE(star.diagnostics[1].coverage_pos.fixed_mesh_fallback);
  EXPECT_TRUE(star.diagnostics[1].coverage_neg.fixed_mesh_fallback);
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    const auto levels = select(star, sign, 1);
    for (unsigned int m = 0; m <= mmax; m++) {
      EXPECT_EQ(levels[m].coupling.norm(), 0.0);
      // At the centre of the interval of the fixed mesh, on the logarithmic scale.
      const auto centre = std::pow(lambda_value, -static_cast<double>(m) - 0.5);
      EXPECT_NEAR(std::abs(levels[m].energy), centre, 1e-14);
    }
  }

  // A scalar Gamma that vanishes at negative frequencies only.
  auto input = make_input<double>([](const double) { return scalar<double>(0.3); });
  for (auto &gamma : input.neg.gamma) gamma = scalar<double>(0.0);
  const auto one_sided = build_star(input, options);
  EXPECT_FALSE(one_sided.diagnostics[0].coverage_pos.fixed_mesh_fallback);
  EXPECT_TRUE(one_sided.diagnostics[0].coverage_neg.fixed_mesh_fallback);
  for (const auto &level : select(one_sided, Sign::NEG, 0)) EXPECT_EQ(level.coupling(0), 0.0);
  EXPECT_GT(one_sided.theta(0, 0), 0.0);

  // Without adapt nothing falls back.
  options.adapt = false;
  EXPECT_FALSE(build_star(input, options).diagnostics[0].coverage_neg.fixed_mesh_fallback);
}

TEST(MixChainStar, without_split_blocks_a_diagonal_gamma_is_one_block) { // NOLINT
  const auto input     = make_input<double>([](const double omega) { return diagonal<double>(0.4 + omega, 0.2); });
  auto options         = base_options();
  options.split_blocks = false;
  const auto star      = build_star(input, options);
  EXPECT_EQ(star.blocks, (Blocks{{0, 1}}));
  EXPECT_EQ(star.diagnostics.size(), 1U);
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
