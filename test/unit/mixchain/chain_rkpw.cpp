#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <mixchain/chain.hpp>
#include <mixchain/chain_lanczos.hpp>
#include <mixchain/chain_rkpw.hpp>
#include <mixchain/precision.hpp>

using namespace NRG::MixChain;

namespace {

// The reference: block Lanczos at 50 digits, on the same star.
using Real    = WideReal<50>;
using Complex = WideComplex<50>;

// The target of the scalar chain qualification, test/CHAIN_QUALIFICATION.md.
constexpr double budget = 2e-12;

constexpr double lambda_value = 2.0;

// The Wilson chain of a flat band at z=1 with the integral-method representative energies: Wilson's closed form
// divided by the Campo-Oliveira factor.
double flat_band_xi(const int n) {
  const auto L = lambda_value;
  return (1.0 - 1.0 / L) / std::log(L) * (1.0 - std::pow(L, -n - 1.0)) * std::pow(L, -n / 2.0)
         / std::sqrt((1.0 - std::pow(L, -2.0 * n - 1.0)) * (1.0 - std::pow(L, -2.0 * n - 3.0)));
}

// Gamma(omega) tabulated for the star stage, with a density of its own on each frequency branch.
GammaInput<double> input_of(const std::function<Matrix<double>(double)> &positive,
                            const std::function<Matrix<double>(double)> &negative) {
  GammaInput<double> input;
  for (int k = 0; k <= 100; k++) {
    input.pos.omega.push_back(0.01 * k);
    input.pos.gamma.push_back(positive(0.01 * k));
    input.neg.omega.push_back(0.01 * k);
    input.neg.gamma.push_back(negative(0.01 * k));
  }
  input.channels = static_cast<int>(input.pos.gamma.front().rows());
  return input;
}

Star<double> star_of(const std::function<Matrix<double>(double)> &positive,
                     const std::function<Matrix<double>(double)> &negative, const unsigned int mmax,
                     const bool split = true) {
  StarOptions options;
  options.Lambda       = NRG::Tools::LambdaCache(lambda_value);
  options.z            = 1.0;
  options.mMAX         = mmax;
  options.split_blocks = split;
  return build_star(input_of(positive, negative), options);
}

Star<double> star_of(const std::function<Matrix<double>(double)> &gamma, const unsigned int mmax,
                     const bool split = true) {
  return star_of(gamma, gamma, mmax, split);
}

template<typename S> double largest(const Matrix<S> &m) { return m.cwiseAbs().maxCoeff(); }

// The whole chain against block Lanczos at 50 digits on the same star: V and the hoppings relative to their largest
// element, the on-site blocks on the scale of their site.
template<typename S0> void expect_matches_lanczos(const Chain<S0> &chain, const Star<S0> &star, const double tolerance) {
  const auto reference = [&] {
    if constexpr (is_complex_v<S0>)
      return convert_chain<S0>(build_chain<Complex>(star, [&] {
        ChainOptions options;
        options.Nmax = chain.Nmax;
        return options;
      }()));
    else
      return convert_chain<S0>(build_chain<Real>(star, [&] {
        ChainOptions options;
        options.Nmax = chain.Nmax;
        return options;
      }()));
  }();
  EXPECT_EQ(chain.diagnostics.theta_rank, reference.diagnostics.theta_rank);
  EXPECT_EQ(chain.diagnostics.hopping_ranks, reference.diagnostics.hopping_ranks);
  EXPECT_LT(largest<S0>(chain.V - reference.V), tolerance * largest(reference.V));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    const auto hopping = largest(reference.T[n]);
    const auto scale   = std::max({largest(reference.E[n]), hopping, n > 0 ? largest(reference.T[n - 1]) : 0.0});
    EXPECT_LE(largest<S0>(chain.T[n] - reference.T[n]), tolerance * hopping) << "site " << n;
    EXPECT_LE(largest<S0>(chain.E[n] - reference.E[n]), tolerance * scale) << "site " << n;
  }
}

std::function<Matrix<double>(double)> scalar(const std::function<double(double)> &density) {
  return [density](const double omega) { return Matrix<double>::Constant(1, 1, density(omega)); };
}

Matrix<double> diagonal_of(const double first, const double second) {
  Matrix<double> m = Matrix<double>::Zero(2, 2);
  m(0, 0)          = first;
  m(1, 1)          = second;
  return m;
}

// A star with arbitrary energies and couplings, not produced by any discretization. Complex couplings when S0 is
// complex.
template<typename S0> Star<S0> arbitrary_star(const int channels, const int levels) {
  Star<S0> star;
  star.channels = channels;
  for (int k = 0; k < levels; k++) {
    StarLevel<S0> level;
    level.energy   = (k % 2 ? -1.0 : 1.0) * std::pow(2.0, -0.5 * k) * (1.0 + 0.1 * std::sin(k));
    level.coupling = Vector<S0>(channels);
    for (int i = 0; i < channels; i++)
      level.coupling(i) =
        make_scalar<S0>(0.3 * std::cos(1.3 * k + i), 0.2 * std::sin(0.7 * k - i)) * std::pow(2.0, -0.25 * k);
    star.levels.push_back(level);
  }
  return star;
}

ChainOptions chain_options(const unsigned int nmax) {
  ChainOptions options;
  options.Nmax = nmax;
  return options;
}

// One channel of the chain against the reference: the hoppings relative to themselves, the on-site energies on the
// scale of their site, since they vanish for a symmetric band.
void expect_matches_reference(const Chain<double> &chain, const Star<double> &star, const int channel = 0,
                              const double tolerance = budget) {
  const auto reference = convert_chain<double>(build_chain<Real>(star, chain_options(chain.Nmax)));
  EXPECT_NEAR(chain.V(channel, channel), reference.V(channel, channel), tolerance * reference.V(channel, channel));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    const auto xi    = reference.T[n](channel, channel);
    const auto scale = std::max({std::abs(reference.E[n](channel, channel)), xi,
                                 n > 0 ? reference.T[n - 1](channel, channel) : 0.0});
    EXPECT_NEAR(chain.T[n](channel, channel), xi, tolerance * xi) << "site " << n;
    EXPECT_NEAR(chain.E[n](channel, channel), reference.E[n](channel, channel), tolerance * scale) << "site " << n;
  }
}

} // namespace

TEST(MixChainRkpw, flat_band_gives_the_closed_form_chain) { // NOLINT
  // mMAX = 4 Nmax, so that the truncation of the star does not show at the end of the chain.
  for (const auto nmax : {20U, 60U}) {
    const auto star  = star_of(scalar([](const double) { return 0.3; }), 4 * nmax);
    const auto chain = build_chain_rkpw(star, chain_options(nmax));
    ASSERT_EQ(chain.T.size(), nmax + 1);
    for (unsigned int n = 0; n <= nmax; n++) {
      const auto xi = flat_band_xi(static_cast<int>(n));
      EXPECT_NEAR(chain.T[n](0, 0), xi, 5e-13 * xi) << "Nmax " << nmax << " site " << n;
      EXPECT_NEAR(chain.E[n](0, 0), 0.0, 5e-13 * xi) << "Nmax " << nmax << " site " << n; // particle-hole symmetry
    }
    // V^2 = Theta, the weight of the star: rho over the covered part of both frequency branches.
    EXPECT_NEAR(chain.V(0, 0) * chain.V(0, 0), 2.0 * 0.3 * (1.0 - std::pow(2.0, -(4.0 * nmax + 1.0))), 1e-14);
    EXPECT_EQ(chain.diagnostics.theta_rank, 1);
    EXPECT_FALSE(chain.diagnostics.rank_drop_site.has_value());
  }
}

TEST(MixChainRkpw, asymmetric_densities_match_the_multiprecision_chain) { // NOLINT
  // Different shapes and weights on the two frequency branches, so that the on-site energies are not zero.
  const auto positive = scalar([](const double omega) { return 0.8 - 0.3 * omega + 0.2 * omega * omega; });
  const auto negative = scalar([](const double omega) { return 0.1 + 0.4 * omega; });
  for (const auto &[mmax, nmax] : {std::pair{40U, 12U}, std::pair{120U, 50U}}) {
    const auto star = star_of(positive, negative, mmax);
    expect_matches_reference(build_chain_rkpw(star, chain_options(nmax)), star);
  }
}

TEST(MixChainRkpw, an_arbitrary_star_matches_the_multiprecision_chain) { // NOLINT
  const auto star = arbitrary_star<double>(1, 40);
  expect_matches_reference(build_chain_rkpw(star, chain_options(12)), star);
}

TEST(MixChainRkpw, a_split_diagonal_gamma_gives_exactly_the_scalar_chains) { // NOLINT
  const auto first  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second = [](const double omega) { return 0.3 + 0.1 * omega * omega; };
  const auto star   = star_of([&](const double omega) { return diagonal_of(first(omega), second(omega)); }, 40);
  const auto joint  = build_chain_rkpw(star, chain_options(8));
  EXPECT_EQ(joint.blocks, (Blocks{{0}, {1}}));
  ASSERT_EQ(joint.block_diagnostics.size(), 2U);

  for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first), 0},
                                         std::pair{std::function<double(double)>(second), 1}}) {
    const auto alone = build_chain_rkpw(star_of(scalar(density), 40), chain_options(8));
    EXPECT_EQ(joint.V(channel, channel), alone.V(0, 0));
    EXPECT_EQ(joint.V(channel, 1 - channel), 0.0);
    for (unsigned int n = 0; n <= joint.Nmax; n++) {
      EXPECT_EQ(joint.E[n](channel, channel), alone.E[n](0, 0)) << "site " << n;
      EXPECT_EQ(joint.E[n](channel, 1 - channel), 0.0);
      EXPECT_EQ(joint.T[n](channel, channel), alone.T[n](0, 0)) << "site " << n;
      EXPECT_EQ(joint.T[n](channel, 1 - channel), 0.0);
    }
    expect_matches_reference(joint, star, channel);
  }
}

TEST(MixChainRkpw, a_zero_block_has_a_zero_chain_and_does_not_spoil_the_diagnostics) { // NOLINT
  const auto chain =
    build_chain_rkpw(star_of([](const double omega) { return diagonal_of(0.4 + omega, 0.0); }, 40), chain_options(8));
  ASSERT_EQ(chain.blocks, (Blocks{{0}, {1}}));
  EXPECT_EQ(chain.diagnostics.theta_rank, 1);
  EXPECT_EQ(chain.diagnostics.theta_condition, 1.0); // from the first block alone
  EXPECT_FALSE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(chain.diagnostics.min_rank, 1);
  EXPECT_EQ(chain.block_diagnostics[1].theta_rank, 0);
  EXPECT_EQ(chain.block_diagnostics[0].levels, 82); // 2 (mMAX+1) per channel
  EXPECT_EQ(chain.block_diagnostics[0].coupled_levels, 82);
  EXPECT_EQ(chain.block_diagnostics[1].levels, 82);
  EXPECT_EQ(chain.block_diagnostics[1].coupled_levels, 0);
  EXPECT_EQ(chain.diagnostics.levels, 164);
  EXPECT_EQ(chain.diagnostics.coupled_levels, 82);
  EXPECT_EQ(chain.V(1, 1), 0.0);
  for (unsigned int n = 0; n <= chain.Nmax; n++) EXPECT_EQ(chain.E[n](1, 1), 0.0);
  for (unsigned int n = 0; n <= chain.Nmax; n++) EXPECT_EQ(chain.T[n](1, 1), 0.0);
}

TEST(MixChainRkpw, an_exhausted_channel_continues_with_zeros) { // NOLINT
  // Only 3 of the 6 levels couple, so the chain ends after 3 sites: the hopping out of site 2 is exactly zero, and so
  // is everything beyond.
  auto star = arbitrary_star<double>(1, 6);
  for (std::size_t k = 3; k < 6; k++) star.levels[k].coupling(0) = 0.0;

  const auto chain = build_chain_rkpw(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 1);
  EXPECT_EQ(chain.diagnostics.min_rank, 0);
  EXPECT_EQ(chain.diagnostics.hopping_ranks, (std::vector<int>{1, 1, 0, 0, 0}));
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 2U);
  EXPECT_EQ(chain.diagnostics.levels, 6);
  EXPECT_EQ(chain.diagnostics.coupled_levels, 3);

  const auto reference = convert_chain<double>(build_chain<Real>(star, chain_options(4)));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    if (n <= 2) {
      EXPECT_NEAR(chain.E[n](0, 0), reference.E[n](0, 0), 1e-14) << "site " << n;
    } else {
      EXPECT_EQ(chain.E[n](0, 0), 0.0) << "site " << n;
    }
    if (n <= 1) {
      EXPECT_NEAR(chain.T[n](0, 0), reference.T[n](0, 0), 1e-14) << "site " << n;
    } else {
      EXPECT_EQ(chain.T[n](0, 0), 0.0) << "site " << n;
    }
  }
}

TEST(MixChainRkpw, ranks_of_blocks_add_up_site_by_site) { // NOLINT
  // Block {1} has 12 coupled levels; block {2} has 2, padded with 10 levels without coupling, so that it passes the
  // count of levels but its chain ends after two sites.
  const auto first  = arbitrary_star<double>(1, 12);
  const auto second = arbitrary_star<double>(1, 2);
  Star<double> star;
  star.channels = 2;
  star.blocks   = {{0}, {1}};
  for (const auto &[part, channel] : {std::pair{&first, 0}, std::pair{&second, 1}})
    for (const auto &scalar_level : part->levels) {
      StarLevel<double> level;
      level.branch            = channel;
      level.energy            = scalar_level.energy;
      level.coupling          = Vector<double>::Zero(2);
      level.coupling(channel) = scalar_level.coupling(0);
      star.levels.push_back(level);
    }
  for (int k = 0; k < 10; k++) {
    StarLevel<double> level;
    level.branch   = 1;
    level.energy   = 0.9 * std::pow(3.0, -k);
    level.coupling = Vector<double>::Zero(2);
    star.levels.push_back(level);
  }

  const auto chain = build_chain_rkpw(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 2);
  EXPECT_EQ(chain.diagnostics.hopping_ranks, (std::vector<int>{2, 1, 1, 1, 1}));
  EXPECT_EQ(chain.diagnostics.min_rank, 1);
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 1U);
  EXPECT_FALSE(chain.block_diagnostics[0].rank_drop_site.has_value());
  ASSERT_TRUE(chain.block_diagnostics[1].rank_drop_site.has_value());
  EXPECT_EQ(*chain.block_diagnostics[1].rank_drop_site, 1U);
  EXPECT_EQ(chain.block_diagnostics[1].coupled_levels, 2);
  EXPECT_EQ(chain.diagnostics.coupled_levels, 14);

  const auto alone = build_chain_rkpw(first, chain_options(4));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_EQ(chain.T[n](0, 0), alone.T[n](0, 0)) << "site " << n;
    if (n > 0) {
      EXPECT_EQ(chain.T[n](1, 1), 0.0) << "site " << n;
    }
  }
}

TEST(MixChainRkpw, the_order_of_the_levels_in_the_star_does_not_matter) { // NOLINT
  // The levels are sorted by interval and frequency branch before they are handed over, so a star in any order gives
  // the same chain to the last bit.
  const auto star = star_of(scalar([](const double omega) { return 0.8 - 0.3 * omega; }),
                            scalar([](const double omega) { return 0.1 + 0.4 * omega; }), 40);
  auto reversed = star;
  std::reverse(reversed.levels.begin(), reversed.levels.end());

  const auto chain = build_chain_rkpw(star, chain_options(12));
  const auto other = build_chain_rkpw(reversed, chain_options(12));
  EXPECT_EQ(other.V(0, 0), chain.V(0, 0));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_EQ(other.E[n](0, 0), chain.E[n](0, 0)) << "site " << n;
    EXPECT_EQ(other.T[n](0, 0), chain.T[n](0, 0)) << "site " << n;
  }
}

TEST(MixChainRkpw, the_phases_of_a_single_channel_drop_out) { // NOLINT
  // One channel has no direction to rotate: only |v_k| enters, and the chain is real.
  const auto complex_star = arbitrary_star<std::complex<double>>(1, 14);
  Star<double> moduli;
  moduli.channels = 1;
  for (const auto &complex_level : complex_star.levels) {
    StarLevel<double> level;
    level.energy   = complex_level.energy;
    level.coupling = Vector<double>::Constant(1, std::abs(complex_level.coupling(0)));
    moduli.levels.push_back(level);
  }

  const auto chain = build_chain_rkpw(complex_star, chain_options(5));
  const auto real  = build_chain_rkpw(moduli, chain_options(5));
  EXPECT_EQ(chain.V(0, 0), std::complex<double>(real.V(0, 0), 0.0));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_EQ(chain.E[n](0, 0), std::complex<double>(real.E[n](0, 0), 0.0)) << "site " << n;
    EXPECT_EQ(chain.T[n](0, 0), std::complex<double>(real.T[n](0, 0), 0.0)) << "site " << n;
  }
}

TEST(MixChainRkpw, rejects_what_it_cannot_map) { // NOLINT
  // A chain of 4 sites needs 5 levels: one per site and one for the hopping out of the last one.
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(1, 4), chain_options(3)), std::invalid_argument);
  EXPECT_NO_THROW(build_chain_rkpw(arbitrary_star<double>(1, 5), chain_options(3)));
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(1, 12), chain_options(0)), std::invalid_argument);
  // With 2 channels, 2*(4+1) = 10 levels.
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(2, 9), chain_options(3)), std::invalid_argument);
  EXPECT_NO_THROW(build_chain_rkpw(arbitrary_star<double>(2, 10), chain_options(3)));
  // The nambu gauge needs blocks of two channels, as with block Lanczos.
  auto options  = chain_options(3);
  options.gauge = ChainGauge::nambu;
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(1, 12), options), std::invalid_argument);
}

// BLOCKS OF SEVERAL CHANNELS

TEST(MixChainRkpw, arbitrary_stars_of_several_channels_match_the_multiprecision_chain) { // NOLINT
  for (const int channels : {2, 3, 4}) {
    const auto real = arbitrary_star<double>(channels, 24 * channels);
    expect_matches_lanczos(build_chain_rkpw(real, chain_options(8)), real, budget);
    const auto complex = arbitrary_star<std::complex<double>>(channels, 24 * channels);
    expect_matches_lanczos(build_chain_rkpw(complex, chain_options(8)), complex, budget);
  }
}

TEST(MixChainRkpw, a_gamma_with_turning_eigenvectors_matches_the_multiprecision_chain) { // NOLINT
  // Two bands of different shape, mixed by an angle that depends on the frequency and differs between the branches.
  const auto gamma = [](const double shift) {
    return [shift](const double omega) {
      const auto angle = 0.4 + shift + 0.8 * omega;
      Matrix<double> u(2, 2);
      u << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
      return Matrix<double>(u * diagonal_of(0.7 - 0.2 * omega, 0.2 + 0.3 * omega * omega) * u.transpose());
    };
  };
  for (const auto &[mmax, nmax] : {std::pair{40U, 12U}, std::pair{100U, 40U}}) {
    const auto star = star_of(gamma(0.0), gamma(0.5), mmax);
    ASSERT_EQ(star.blocks.size(), 1U);
    expect_matches_lanczos(build_chain_rkpw(star, chain_options(nmax)), star, budget);
  }
}

TEST(MixChainRkpw, degenerate_flat_band_gives_the_scalar_chain_times_the_identity) { // NOLINT
  // Kept whole: the two branches of every interval share their energy.
  const auto star =
    star_of([](const double) { return Matrix<double>(0.3 * Matrix<double>::Identity(2, 2)); }, 80, false);
  const auto chain = build_chain_rkpw(star, chain_options(12));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    const auto xi = flat_band_xi(static_cast<int>(n));
    EXPECT_LT(largest<double>(chain.T[n] - xi * Matrix<double>::Identity(2, 2)), 1e-13 * xi) << "site " << n;
    EXPECT_LT(largest(chain.E[n]), 1e-13 * xi) << "site " << n;
  }
}

TEST(MixChainRkpw, is_covariant_under_a_constant_rotation) { // NOLINT
  const auto densities = [](const double omega) { return diagonal_of(0.5 + 0.1 * omega, 0.2 + 0.05 * omega); };
  Matrix<double> u(2, 2);
  u << std::cos(0.7), -std::sin(0.7), std::sin(0.7), std::cos(0.7);

  const auto plain   = build_chain_rkpw(star_of(densities, 40, false), chain_options(8));
  const auto rotated = build_chain_rkpw(
    star_of([&](const double omega) { return Matrix<double>(u * densities(omega) * u.transpose()); }, 40),
    chain_options(8));
  // Every block rotates with U; the polar gauge involves no preferred basis.
  EXPECT_LT(largest<double>(rotated.V - u * plain.V * u.transpose()), 1e-11);
  for (unsigned int n = 0; n <= plain.Nmax; n++) {
    EXPECT_LT(largest<double>(rotated.E[n] - u * plain.E[n] * u.transpose()), 1e-11) << "site " << n;
    EXPECT_LT(largest<double>(rotated.T[n] - u * plain.T[n] * u.transpose()), 1e-11) << "site " << n;
  }
}

TEST(MixChainRkpw, a_diagonal_gamma_kept_whole_gives_independent_scalar_chains) { // NOLINT
  const auto first  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second = [](const double omega) { return 0.3 + 0.1 * omega; };
  const auto joint  = build_chain_rkpw(
    star_of([&](const double omega) { return diagonal_of(first(omega), second(omega)); }, 40, false), chain_options(8));
  ASSERT_EQ(joint.blocks.size(), 1U);

  for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first), 0},
                                         std::pair{std::function<double(double)>(second), 1}}) {
    const auto alone = build_chain_rkpw(star_of(scalar(density), 40), chain_options(8));
    for (unsigned int n = 0; n <= joint.Nmax; n++) {
      const auto xi = alone.T[n](0, 0);
      EXPECT_NEAR(joint.T[n](channel, channel), xi, 1e-12 * xi) << "site " << n;
      EXPECT_LT(std::abs(joint.T[n](0, 1)), 1e-13 * xi) << "site " << n; // the channels stay decoupled
      EXPECT_NEAR(joint.E[n](channel, channel), alone.E[n](0, 0), 1e-12) << "site " << n;
      EXPECT_LT(std::abs(joint.E[n](0, 1)), 1e-13) << "site " << n;
    }
  }
}

TEST(MixChainRkpw, a_singular_theta_gives_a_zero_chain_for_the_decoupled_combination) { // NOLINT
  // Every coupling is c_k (1, i): the chiral case, where Theta has rank 1 of 2. Along u = (1, i)/sqrt(2) the chain is
  // that of the scalar star with couplings c_k, with V scaled by sqrt(2): every block is the scalar one times the
  // projector u u^dag.
  const auto scalar_star = arbitrary_star<double>(1, 12);
  Star<std::complex<double>> star;
  star.channels = 2;
  for (const auto &scalar_level : scalar_star.levels) {
    StarLevel<std::complex<double>> level;
    level.energy   = scalar_level.energy;
    level.coupling = Vector<std::complex<double>>(2);
    level.coupling << std::complex<double>(scalar_level.coupling(0), 0.0), std::complex<double>(0.0, scalar_level.coupling(0));
    star.levels.push_back(level);
  }
  Matrix<std::complex<double>> projector(2, 2);
  projector << 0.5, std::complex<double>(0.0, -0.5), std::complex<double>(0.0, 0.5), 0.5;

  const auto reference = build_chain_rkpw(scalar_star, chain_options(4));
  const auto chain     = build_chain_rkpw(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 1);
  EXPECT_EQ(chain.diagnostics.min_rank, 1);
  EXPECT_FALSE(chain.diagnostics.rank_drop_site.has_value());
  using M = Matrix<std::complex<double>>;
  EXPECT_LT(largest<std::complex<double>>(chain.V - M(std::sqrt(2.0) * reference.V(0, 0) * projector)), 1e-13);
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_LT(largest<std::complex<double>>(chain.E[n] - M(reference.E[n](0, 0) * projector)), 1e-13) << "site " << n;
    EXPECT_LT(largest<std::complex<double>>(chain.T[n] - M(reference.T[n](0, 0) * projector)), 1e-13) << "site " << n;
  }
}

TEST(MixChainRkpw, a_rank_drop_mid_chain_continues_with_zeros) { // NOLINT
  // A diagonal star kept whole: channel 0 has 12 levels, channel 1 only 2, so the Krylov space of channel 1 runs out
  // after two sites. The hopping T_1 has rank 1, and the chain of channel 1 is zero from site 2 on. The band matrix of
  // the rotations does not show this by itself; the second stage restores it. Also for the same star turned by a
  // constant angle, where no element is exactly zero.
  const auto first  = arbitrary_star<double>(1, 12);
  const auto second = arbitrary_star<double>(1, 2);
  Star<double> star;
  star.channels = 2;
  for (const auto &[part, channel] : {std::pair{&first, 0}, std::pair{&second, 1}})
    for (const auto &scalar_level : part->levels) {
      StarLevel<double> level;
      level.energy            = scalar_level.energy;
      level.coupling          = Vector<double>::Zero(2);
      level.coupling(channel) = scalar_level.coupling(0);
      star.levels.push_back(level);
    }
  Matrix<double> u(2, 2);
  u << std::cos(0.7), -std::sin(0.7), std::sin(0.7), std::cos(0.7);
  auto turned = star;
  for (auto &level : turned.levels) level.coupling = (u * level.coupling).eval();

  const auto alone_first = build_chain_rkpw(first, chain_options(4));
  auto second_padded     = second; // padded with levels without coupling, so that a chain can be asked for at all
  for (int k = 0; k < 2; k++) {
    StarLevel<double> level;
    level.energy   = 0.9 * std::pow(3.0, -k);
    level.coupling = Vector<double>::Zero(1);
    second_padded.levels.push_back(level);
  }
  const auto alone_second = build_chain_rkpw(second_padded, chain_options(1));

  for (const bool rotate : {false, true}) {
    const auto chain = build_chain_rkpw(rotate ? turned : star, chain_options(4));
    EXPECT_EQ(chain.diagnostics.theta_rank, 2);
    EXPECT_EQ(chain.diagnostics.hopping_ranks, (std::vector<int>{2, 1, 1, 1, 1}));
    EXPECT_EQ(chain.diagnostics.min_rank, 1);
    ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
    EXPECT_EQ(*chain.diagnostics.rank_drop_site, 1U);
    const Matrix<double> back = rotate ? Matrix<double>(u.transpose()) : Matrix<double>(Matrix<double>::Identity(2, 2));
    for (unsigned int n = 0; n <= chain.Nmax; n++) {
      const Matrix<double> onsite  = back * chain.E[n] * back.transpose();
      const Matrix<double> hopping = back * chain.T[n] * back.transpose();
      EXPECT_NEAR(onsite(0, 0), alone_first.E[n](0, 0), 1e-13) << "site " << n;
      EXPECT_NEAR(onsite(1, 1), n <= 1 ? alone_second.E[n](0, 0) : 0.0, 1e-13) << "site " << n;
      EXPECT_LT(std::abs(onsite(0, 1)), 1e-13) << "site " << n;
      EXPECT_NEAR(hopping(0, 0), alone_first.T[n](0, 0), 1e-13) << "site " << n;
      EXPECT_NEAR(hopping(1, 1), n == 0 ? alone_second.T[0](0, 0) : 0.0, 1e-13) << "site " << n;
      EXPECT_LT(std::abs(hopping(0, 1)), 1e-13) << "site " << n;
    }
    expect_matches_lanczos(chain, rotate ? turned : star, 1e-12);
  }
}

TEST(MixChainRkpw, a_star_with_too_few_coupled_levels_ends_early) { // NOLINT
  // 7 of the 12 levels couple: three full sites of two channels and one orbital of the fourth.
  auto star = arbitrary_star<double>(2, 12);
  for (std::size_t k = 7; k < 12; k++) star.levels[k].coupling.setZero();
  const auto chain = build_chain_rkpw(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.coupled_levels, 7);
  EXPECT_EQ(chain.diagnostics.hopping_ranks, (std::vector<int>{2, 2, 1, 0, 0}));
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 2U);
  EXPECT_EQ(largest(chain.T[3]), 0.0);
  EXPECT_EQ(largest(chain.E[4]), 0.0);
  expect_matches_lanczos(chain, star, 1e-12);
}

TEST(MixChainRkpw, a_weak_channel_below_the_rank_tolerance_counts_as_decoupled) { // NOLINT
  // The eigenvalues of Theta are 25 orders of magnitude apart, below the rank tolerance of 1e-20.
  const auto whole = build_chain_rkpw(star_of([](const double) { return diagonal_of(0.3, 0.3e-25); }, 40, false),
                                      chain_options(8));
  EXPECT_EQ(whole.diagnostics.theta_rank, 1);
  EXPECT_EQ(whole.V(1, 1), 0.0);
  for (unsigned int n = 0; n <= whole.Nmax; n++) EXPECT_EQ(whole.T[n](1, 1), 0.0);
}

TEST(MixChainRkpw, blocks_are_placed_in_their_channels) { // NOLINT
  // Channels 1 and 3 are coupled, channel 2 is on its own.
  GammaInput<std::complex<double>> input;
  input.channels = 3;
  for (int k = 0; k <= 100; k++) {
    const auto omega = 0.01 * k;
    Matrix<std::complex<double>> m = Matrix<std::complex<double>>::Zero(3, 3);
    m(0, 0) = 0.5 + 0.2 * omega;
    m(1, 1) = 0.3 + omega * omega;
    m(2, 2) = 0.4 - 0.1 * omega;
    m(0, 2) = std::complex<double>(0.1 * omega, 0.05);
    m(2, 0) = std::conj(m(0, 2));
    for (auto *branch : {&input.pos, &input.neg}) {
      branch->omega.push_back(omega);
      branch->gamma.push_back(m);
    }
  }
  StarOptions options;
  options.Lambda = NRG::Tools::LambdaCache(lambda_value);
  options.z      = 1.0;
  options.mMAX   = 20;
  const auto star = build_star(input, options);
  ASSERT_EQ(star.blocks, (Blocks{{0, 2}, {1}}));

  const auto chain = build_chain_rkpw(star, chain_options(6));
  EXPECT_EQ(chain.blocks, star.blocks);
  const auto zero  = std::complex<double>(0.0, 0.0);
  const auto check = [&zero](const Matrix<std::complex<double>> &m) {
    for (const int i : {0, 2}) {
      EXPECT_EQ(m(i, 1), zero);
      EXPECT_EQ(m(1, i), zero);
    }
    EXPECT_EQ(m(1, 1).imag(), 0.0);
  };
  check(chain.V);
  for (const auto &block : chain.E) check(block);
  for (const auto &block : chain.T) check(block);
  expect_matches_lanczos(chain, star, budget);
}

TEST(MixChainRkpw, the_nambu_gauge_puts_the_blocks_into_the_nambu_structure) { // NOLINT
  // A Nambu-symmetric bath: the normal part is flat and equal for the particle and the hole, and the anomalous part
  // is odd in omega.
  const double rho = 0.3, anomalous = 0.1;
  const auto block = [](const double diagonal, const double offdiagonal) {
    Matrix<double> m(2, 2);
    m << diagonal, offdiagonal, offdiagonal, diagonal;
    return m;
  };
  const auto star = star_of([&](const double) { return block(rho, anomalous); },
                            [&](const double) { return block(rho, -anomalous); }, 40);
  ASSERT_EQ(star.blocks.size(), 1U);

  auto options     = chain_options(8);
  const auto polar = build_chain_rkpw(star, options);
  options.gauge    = ChainGauge::nambu;
  const auto nambu = build_chain_rkpw(star, options);
  EXPECT_EQ(nambu.gauge, ChainGauge::nambu);
  EXPECT_LT(nambu.diagnostics.max_nambu_deviation, 1e-11);
  EXPECT_LT(std::abs(nambu.V(1, 1) + nambu.V(0, 0)), 1e-13);
  for (unsigned int n = 0; n <= nambu.Nmax; n++) {
    const auto scale = largest(polar.T[n]);
    EXPECT_LT(std::abs(nambu.E[n](1, 1) + nambu.E[n](0, 0)), 1e-11 * scale) << "site " << n;
    EXPECT_LT(std::abs(nambu.T[n](1, 1) + nambu.T[n](0, 0)), 1e-11 * scale) << "site " << n;
    EXPECT_EQ(nambu.T[n](0, 0), polar.T[n](0, 0)) << "site " << n; // the gauge only flips signs
    EXPECT_EQ(nambu.T[n](1, 1), -polar.T[n](1, 1)) << "site " << n;
  }
  expect_matches_lanczos(polar, star, budget);
  // A generic 2x2 star is one block, but its chain has no Nambu structure.
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(2, 12), [] {
    auto refused  = chain_options(3);
    refused.gauge = ChainGauge::nambu;
    return refused;
  }()), std::runtime_error);
}

TEST(MixChainRkpw, real_data_in_complex_arithmetic_stays_real) { // NOLINT
  const auto real = arbitrary_star<double>(2, 24);
  Star<std::complex<double>> star;
  star.channels = 2;
  for (const auto &real_level : real.levels) {
    StarLevel<std::complex<double>> level;
    level.energy   = real_level.energy;
    level.coupling = real_level.coupling.cast<std::complex<double>>();
    star.levels.push_back(level);
  }
  const auto chain    = build_chain_rkpw(star, chain_options(5));
  const auto expected = build_chain_rkpw(real, chain_options(5));
  const auto check    = [](const Matrix<std::complex<double>> &m, const Matrix<double> &reference) {
    EXPECT_LT(largest<std::complex<double>>(m - reference.cast<std::complex<double>>()), 1e-13);
  };
  check(chain.V, expected.V);
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    check(chain.E[n], expected.E[n]);
    check(chain.T[n], expected.T[n]);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
