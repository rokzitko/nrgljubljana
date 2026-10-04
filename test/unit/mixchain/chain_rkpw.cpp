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
using Real = WideReal<50>;

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
                     const std::function<Matrix<double>(double)> &negative, const unsigned int mmax) {
  StarOptions options;
  options.Lambda = NRG::Tools::LambdaCache(lambda_value);
  options.z      = 1.0;
  options.mMAX   = mmax;
  return build_star(input_of(positive, negative), options);
}

Star<double> star_of(const std::function<Matrix<double>(double)> &gamma, const unsigned int mmax) {
  return star_of(gamma, gamma, mmax);
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
  // Blocks of several channels are not handled yet.
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(2, 12), chain_options(3)), std::invalid_argument);
  // The nambu gauge needs blocks of two channels, as with block Lanczos.
  auto options  = chain_options(3);
  options.gauge = ChainGauge::nambu;
  EXPECT_THROW(build_chain_rkpw(arbitrary_star<double>(1, 12), options), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
