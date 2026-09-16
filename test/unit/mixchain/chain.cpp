#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <functional>
#include <optional>
#include <stdexcept>
#include <vector>

#include <mixchain/chain.hpp>
#include <mixchain/precision.hpp>

using namespace NRG::MixChain;

namespace {

using Real    = WideReal<50>;
using Complex = WideComplex<50>;

constexpr double lambda_value = 2.0;

// The Wilson chain of a flat band at z=1 with the integral-method representative energies: Wilson's closed form
// divided by the Campo-Oliveira factor.
double flat_band_xi(const int n) {
  const auto L = lambda_value;
  return (1.0 - 1.0 / L) / std::log(L) * (1.0 - std::pow(L, -n - 1.0)) * std::pow(L, -n / 2.0)
         / std::sqrt((1.0 - std::pow(L, -2.0 * n - 1.0)) * (1.0 - std::pow(L, -2.0 * n - 3.0)));
}

// Gamma(omega) tabulated for the star stage, the same on both frequency branches.
template<typename S0> GammaInput<S0> input_of(const std::function<Matrix<S0>(double)> &gamma) {
  GammaBranch<S0> branch;
  for (int k = 0; k <= 100; k++) {
    branch.omega.push_back(0.01 * k);
    branch.gamma.push_back(gamma(0.01 * k));
  }
  GammaInput<S0> input;
  input.channels = static_cast<int>(branch.gamma.front().rows());
  input.pos      = branch;
  input.neg      = branch;
  return input;
}

template<typename S0>
Star<S0> star_from(const GammaInput<S0> &input, const unsigned int mmax, const bool split = true) {
  StarOptions options;
  options.Lambda       = NRG::Tools::LambdaCache(lambda_value);
  options.z            = 1.0;
  options.mMAX         = mmax;
  options.split_blocks = split;
  return build_star(input, options);
}

// A star built by the star stage from Gamma(omega).
template<typename S0>
Star<S0> star_of(const std::function<Matrix<S0>(double)> &gamma, const unsigned int mmax, const bool split = true) {
  return star_from(input_of(gamma), mmax, split);
}

Matrix<double> diagonal_of(const double first, const double second) {
  Matrix<double> m = Matrix<double>::Zero(2, 2);
  m(0, 0)          = first;
  m(1, 1)          = second;
  return m;
}

// A star with arbitrary energies and couplings, not produced by any discretization, so that nothing about it is
// special. Complex couplings when S0 is complex.
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

template<typename S> double largest(const Matrix<S> &m) { return static_cast<double>(m.cwiseAbs().maxCoeff()); }

// The chain Hamiltonian as one block tridiagonal matrix: E_n on the diagonal, T_n below it and T_n^dag above.
template<typename S> Matrix<S> chain_hamiltonian(const Chain<S> &chain) {
  const auto n     = static_cast<Eigen::Index>(chain.channels);
  const auto sites = static_cast<Eigen::Index>(chain.Nmax + 1);
  Matrix<S> h      = Matrix<S>::Zero(n * sites, n * sites);
  for (Eigen::Index s = 0; s < sites; s++) {
    h.block(s * n, s * n, n, n) = chain.E[static_cast<std::size_t>(s)];
    if (s + 1 < sites) {
      h.block((s + 1) * n, s * n, n, n) = chain.T[static_cast<std::size_t>(s)];
      h.block(s * n, (s + 1) * n, n, n) = chain.T[static_cast<std::size_t>(s)].adjoint();
    }
  }
  return h;
}

} // namespace

TEST(MixChainChain, flat_band_gives_the_closed_form_chain) { // NOLINT
  // mMAX = 4 Nmax: with the usual 2 Nmax the truncation of the star shows at the end of the chain, at 1e-10 for Nmax=20.
  const auto star  = star_of<double>([](const double) { return Matrix<double>::Constant(1, 1, 0.3); }, 80);
  const auto chain = build_chain<Real>(star, chain_options(20));

  for (unsigned int n = 0; n <= chain.Nmax; n++) { // one hopping per site, the last one out of the chain
    const auto xi = flat_band_xi(static_cast<int>(n));
    EXPECT_NEAR(static_cast<double>(chain.T[n](0, 0)), xi, 1e-14 * xi) << "site " << n;
  }
  for (const auto &onsite : chain.E) EXPECT_LT(largest(onsite), 1e-40); // particle-hole symmetry
  // V^2 = Theta, the weight of the star: rho over the covered part of both frequency branches.
  EXPECT_NEAR(static_cast<double>(chain.V(0, 0) * chain.V(0, 0)), 2.0 * 0.3 * (1.0 - std::pow(2.0, -81.0)), 1e-14);
}

TEST(MixChainChain, degenerate_flat_band_gives_the_scalar_chain_times_the_identity) { // NOLINT
  const auto star  = star_of<double>([](const double) { return Matrix<double>(0.3 * Matrix<double>::Identity(2, 2)); }, 80);
  const auto chain = build_chain<Real>(star, chain_options(12));

  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    const auto xi = flat_band_xi(static_cast<int>(n));
    const Matrix<Real> expected = Matrix<Real>::Identity(2, 2) * Real(xi);
    EXPECT_LT(largest<Real>(chain.T[n] - expected), 1e-14 * xi) << "site " << n;
  }
  for (const auto &onsite : chain.E) EXPECT_LT(largest(onsite), 1e-40);
}

TEST(MixChainChain, moments_of_the_star_are_reproduced) { // NOLINT
  // Block Lanczos with Nmax+1 blocks reproduces sum_k v_k E_k^p v_k^dag = V (H_chain^p)_00 V for p <= 2 Nmax + 1.
  for (const bool complex_data : {false, true}) {
    const auto check = [&]<typename S, typename S0>() {
      const auto star  = arbitrary_star<S0>(2, 12);
      const auto chain = build_chain<S>(star, chain_options(3));
      const auto wide  = to_wide<S>(star);
      const auto h     = chain_hamiltonian(chain);

      Matrix<S> power = Matrix<S>::Identity(h.rows(), h.cols());
      for (unsigned int p = 0; p <= 2 * chain.Nmax + 1; p++) {
        Matrix<S> exact = Matrix<S>::Zero(2, 2);
        for (Eigen::Index k = 0; k < wide.start.rows(); k++) {
          auto energy_power = make_scalar<S>(1, 0);
          for (unsigned int q = 0; q < p; q++) energy_power *= make_scalar<S>(wide.energies[static_cast<std::size_t>(k)], 0);
          const Vector<S> a = wide.start.row(k).transpose(); // conj(v_k)
          exact += energy_power * (a.conjugate() * a.transpose());
        }
        const Matrix<S> from_chain = chain.V * power.topLeftCorner(2, 2) * chain.V;
        EXPECT_LT(largest<S>(from_chain - exact), 1e-45 * largest(exact)) << "moment " << p;
        power = (power * h).eval();
      }
    };
    if (complex_data)
      check.template operator()<Complex, std::complex<double>>();
    else
      check.template operator()<Real, double>();
  }
}

TEST(MixChainChain, the_lanczos_blocks_span_the_star_and_obey_the_recursion) { // NOLINT
  // With exactly channels*(Nmax+2) levels the blocks Q_0..Q_{Nmax+1} span the whole space, so stacked side by side
  // they are unitary. They must also satisfy the recursion they came from,
  //
  //   H_star Q_n = Q_{n-1} T_{n-1}^dag + Q_n E_n + Q_{n+1} T_n,
  //
  // which is what fixes the coefficients: it checks every block including the last hopping T_Nmax, and it catches what
  // equality of spectra would not, such as a transposed coupling.
  const auto star = arbitrary_star<std::complex<double>>(2, 10);
  std::vector<Matrix<Complex>> blocks;
  const auto chain = build_chain(to_wide<Complex>(star), chain_options(3), &blocks);
  ASSERT_EQ(blocks.size(), chain.Nmax + 2);
  ASSERT_EQ(chain.T.size(), chain.Nmax + 1);

  Matrix<Complex> q(10, 10);
  for (std::size_t s = 0; s < blocks.size(); s++) q.block(0, static_cast<Eigen::Index>(2 * s), 10, 2) = blocks[s];
  EXPECT_LT(largest<Complex>(q.adjoint() * q - Matrix<Complex>::Identity(10, 10)), 1e-45);

  Matrix<Complex> star_hamiltonian = Matrix<Complex>::Zero(10, 10);
  for (std::size_t k = 0; k < star.levels.size(); k++)
    star_hamiltonian(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(k)) =
      make_scalar<Complex>(star.levels[k].energy, 0);

  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    Matrix<Complex> recursion = blocks[n] * chain.E[n] + blocks[n + 1] * chain.T[n];
    if (n > 0) recursion += blocks[n - 1] * chain.T[n - 1].adjoint();
    EXPECT_LT(largest<Complex>(star_hamiltonian * blocks[n] - recursion), 1e-45) << "site " << n;
  }
  // The first Nmax+1 blocks carry the chain itself: Q^dag H_star Q is its block tridiagonal Hamiltonian.
  const Matrix<Complex> spanned = q.leftCols(8);
  EXPECT_LT(largest<Complex>(spanned.adjoint() * star_hamiltonian * spanned - chain_hamiltonian(chain)), 1e-45);
}

TEST(MixChainChain, is_covariant_under_a_constant_rotation) { // NOLINT
  const auto densities = [](const double omega) {
    Matrix<double> m = Matrix<double>::Zero(2, 2);
    m(0, 0)          = 0.5 + 0.1 * omega;
    m(1, 1)          = 0.2 + 0.05 * omega;
    return m;
  };
  Matrix<double> u = Matrix<double>::Zero(2, 2);
  u(0, 0) = std::cos(0.7);
  u(0, 1) = -std::sin(0.7);
  u(1, 0) = std::sin(0.7);
  u(1, 1) = std::cos(0.7);

  const auto plain   = build_chain<Real>(star_of<double>(densities, 40), chain_options(8));
  const auto rotated = build_chain<Real>(
    star_of<double>([&](const double omega) { return Matrix<double>(u * densities(omega) * u.transpose()); }, 40),
    chain_options(8));

  // Every block rotates with U; the polar gauge involves no preferred basis.
  const Matrix<Real> U = u.cast<Real>();
  EXPECT_LT(largest<Real>(rotated.V - U * plain.V * U.transpose()), 1e-11);
  for (unsigned int n = 0; n <= plain.Nmax; n++)
    EXPECT_LT(largest<Real>(rotated.E[n] - U * plain.E[n] * U.transpose()), 1e-11) << "site " << n;
  for (unsigned int n = 0; n <= plain.Nmax; n++)
    EXPECT_LT(largest<Real>(rotated.T[n] - U * plain.T[n] * U.transpose()), 1e-11) << "site " << n;
}

TEST(MixChainChain, a_diagonal_gamma_gives_independent_scalar_chains) { // NOLINT
  // Kept whole, so that the block Lanczos of the whole matrix is what is tested; split_blocks is tested below.
  const auto first  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second = [](const double omega) { return 0.3 + 0.1 * omega; };
  const auto joint  = build_chain<Real>(
    star_of<double>([&](const double omega) { return diagonal_of(first(omega), second(omega)); }, 40, false),
    chain_options(8));

  for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first), 0},
                                         std::pair{std::function<double(double)>(second), 1}}) {
    const auto alone = build_chain<Real>(
      star_of<double>([&](const double omega) { return Matrix<double>::Constant(1, 1, density(omega)); }, 40),
      chain_options(8));
    for (unsigned int n = 0; n <= joint.Nmax; n++) {
      EXPECT_NEAR(static_cast<double>(joint.T[n](channel, channel)), static_cast<double>(alone.T[n](0, 0)),
                  1e-12 * static_cast<double>(alone.T[n](0, 0)));
      EXPECT_LT(static_cast<double>(abs(joint.T[n](0, 1))), 1e-40); // the channels stay decoupled
      EXPECT_NEAR(static_cast<double>(joint.E[n](channel, channel)), static_cast<double>(alone.E[n](0, 0)), 1e-12);
    }
  }
}

TEST(MixChainChain, real_data_in_complex_arithmetic_stays_real) { // NOLINT
  const auto star  = arbitrary_star<double>(2, 12);
  const auto chain = build_chain<Complex>(star, chain_options(3));
  const auto imaginary = [](const Matrix<Complex> &m) {
    double worst = 0;
    for (Eigen::Index i = 0; i < m.rows(); i++)
      for (Eigen::Index j = 0; j < m.cols(); j++) worst = std::max(worst, static_cast<double>(abs(m(i, j).imag())));
    return worst;
  };
  EXPECT_LT(imaginary(chain.V), 1e-45);
  for (const auto &block : chain.E) EXPECT_LT(imaginary(block), 1e-45);
  for (const auto &block : chain.T) EXPECT_LT(imaginary(block), 1e-45);
}

TEST(MixChainChain, a_singular_theta_gives_a_zero_chain_for_the_decoupled_combination) { // NOLINT
  // Every coupling is c_k (1, i) = sqrt(2) c_k u with u = (1, i)/sqrt(2): the chiral case, where Theta has rank 1 of 2.
  // The combination orthogonal to u does not couple to the bath. Along u the chain is that of the scalar star with
  // couplings c_k, with V scaled by sqrt(2): every block is the scalar one times the projector u u^dag.
  const auto scalar = arbitrary_star<double>(1, 12);
  Star<std::complex<double>> star;
  star.channels = 2;
  for (const auto &scalar_level : scalar.levels) {
    StarLevel<std::complex<double>> level;
    level.energy   = scalar_level.energy;
    level.coupling = Vector<std::complex<double>>(2);
    level.coupling << std::complex<double>(scalar_level.coupling(0), 0.0), std::complex<double>(0.0, scalar_level.coupling(0));
    star.levels.push_back(level);
  }
  Matrix<std::complex<double>> projector(2, 2);
  projector << 0.5, std::complex<double>(0.0, -0.5), std::complex<double>(0.0, 0.5), 0.5;

  const auto reference = build_chain<Real>(scalar, chain_options(4));
  const auto check     = [&]<typename S>(const double tolerance) {
    const auto chain       = build_chain<S>(star, chain_options(4));
    const Matrix<S> u      = projector.cast<S>();
    const auto scalar_of   = [](const Real &x) { return make_scalar<S>(static_cast<real_type<S>>(x), 0); };
    EXPECT_EQ(chain.diagnostics.theta_rank, 1);
    EXPECT_EQ(chain.diagnostics.min_rank, 1);
    EXPECT_FALSE(chain.diagnostics.rank_drop_site.has_value());
    EXPECT_LT(largest<S>(chain.V - scalar_of(sqrt(Real(2)) * reference.V(0, 0)) * u), tolerance);
    for (unsigned int n = 0; n <= chain.Nmax; n++)
      EXPECT_LT(largest<S>(chain.E[n] - scalar_of(reference.E[n](0, 0)) * u), tolerance) << "site " << n;
    for (unsigned int n = 0; n <= chain.Nmax; n++)
      EXPECT_LT(largest<S>(chain.T[n] - scalar_of(reference.T[n](0, 0)) * u), tolerance) << "site " << n;
  };
  check.template operator()<Complex>(1e-40);
  // In double precision rounding leaves a small but nonzero eigenvalue of Theta; the tolerance floor still catches it.
  check.template operator()<std::complex<double>>(1e-13);
}

TEST(MixChainChain, a_rank_drop_mid_chain_continues_with_zeros) { // NOLINT
  // A diagonal star: channel 0 has 12 levels, channel 1 only 2, so the Krylov space of channel 1 runs out after two
  // sites. The hopping T_1 has rank 1, and the chain of channel 1 is zero from site 2 on.
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

  const auto chain = build_chain<Real>(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 2);
  EXPECT_EQ(chain.diagnostics.min_rank, 1);
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 1U);

  const auto alone_first  = build_chain<Real>(first, chain_options(4));
  // The two levels of channel 2 alone, padded with levels without coupling so that the star is large enough to ask
  // for a chain at all; they add nothing to the Krylov space.
  auto second_padded = second;
  for (int k = 0; k < 2; k++) {
    StarLevel<double> level;
    level.energy   = 0.9 * std::pow(3.0, -k);
    level.coupling = Vector<double>::Zero(1);
    second_padded.levels.push_back(level);
  }
  const auto alone_second = build_chain<Real>(second_padded, chain_options(1));
  const auto near         = [](const Real &a, const Real &b) { return static_cast<double>(abs(a - b)); };
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_LT(near(chain.E[n](0, 0), alone_first.E[n](0, 0)), 1e-40) << "site " << n;
    EXPECT_LT(near(chain.E[n](1, 1), n <= 1 ? alone_second.E[n](0, 0) : Real(0)), 1e-40) << "site " << n;
    EXPECT_LT(static_cast<double>(abs(chain.E[n](0, 1))), 1e-40) << "site " << n;
  }
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_LT(near(chain.T[n](0, 0), alone_first.T[n](0, 0)), 1e-40) << "site " << n;
    EXPECT_LT(near(chain.T[n](1, 1), n == 0 ? alone_second.T[0](0, 0) : Real(0)), 1e-40) << "site " << n;
    EXPECT_LT(static_cast<double>(abs(chain.T[n](0, 1))), 1e-40) << "site " << n;
  }
}

TEST(MixChainChain, an_exhausted_single_channel_gives_rank_zero) { // NOLINT
  // Only 3 of the 6 levels couple, so a scalar chain ends after 3 sites. Every direction of the residual at site 2 is
  // rounding, which the relative test alone would not see.
  auto star = arbitrary_star<double>(1, 6);
  for (std::size_t k = 3; k < 6; k++) star.levels[k].coupling(0) = 0.0;
  auto coupled = star;
  coupled.levels.resize(4); // the 3 that couple, and one that does not, so that a chain can be asked for at all

  const auto chain = build_chain<Real>(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 1);
  EXPECT_EQ(chain.diagnostics.min_rank, 0);
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 2U);

  const auto alone = build_chain<Real>(coupled, chain_options(2));
  for (unsigned int n = 0; n <= chain.Nmax; n++)
    EXPECT_LT(static_cast<double>(abs(chain.E[n](0, 0) - (n <= 2 ? alone.E[n](0, 0) : Real(0)))), 1e-40) << "site " << n;
  for (unsigned int n = 0; n <= chain.Nmax; n++)
    EXPECT_LT(static_cast<double>(abs(chain.T[n](0, 0) - (n <= 2 ? alone.T[n](0, 0) : Real(0)))), 1e-40) << "site " << n;
}

TEST(MixChainChain, a_split_diagonal_gamma_gives_exactly_the_scalar_chains) { // NOLINT
  const auto first  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second = [](const double omega) { return 0.3 + 0.1 * omega * omega; };
  const auto joint  = build_chain<Real>(
    star_of<double>([&](const double omega) { return diagonal_of(first(omega), second(omega)); }, 40),
    chain_options(8));
  EXPECT_EQ(joint.blocks, (Blocks{{0}, {1}}));
  ASSERT_EQ(joint.block_diagnostics.size(), 2U);

  for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first), 0},
                                         std::pair{std::function<double(double)>(second), 1}}) {
    const auto alone = build_chain<Real>(
      star_of<double>([&](const double omega) { return Matrix<double>::Constant(1, 1, density(omega)); }, 40),
      chain_options(8));
    EXPECT_EQ(joint.V(channel, channel), alone.V(0, 0));
    EXPECT_EQ(joint.V(channel, 1 - channel), Real(0));
    for (unsigned int n = 0; n <= joint.Nmax; n++) {
      EXPECT_EQ(joint.E[n](channel, channel), alone.E[n](0, 0)) << "site " << n;
      EXPECT_EQ(joint.E[n](channel, 1 - channel), Real(0));
    }
    for (unsigned int n = 0; n <= joint.Nmax; n++) {
      EXPECT_EQ(joint.T[n](channel, channel), alone.T[n](0, 0)) << "site " << n;
      EXPECT_EQ(joint.T[n](channel, 1 - channel), Real(0));
    }
  }
}

TEST(MixChainChain, blocks_are_compared_only_within_themselves) { // NOLINT
  // The eigenvalues of Theta are 25 orders of magnitude apart, below the rank tolerance of 1e-20.
  const auto gamma = [](const double) { return diagonal_of(0.3, 0.3e-25); };

  const auto split = build_chain<Real>(star_of<double>(gamma, 40), chain_options(8));
  EXPECT_EQ(split.diagnostics.theta_rank, 2);
  EXPECT_FALSE(split.diagnostics.rank_drop_site.has_value());
  EXPECT_NEAR(static_cast<double>(split.V(1, 1) / split.V(0, 0)), std::sqrt(1e-25), 1e-14 * std::sqrt(1e-25));
  for (unsigned int n = 0; n <= split.Nmax; n++) // the hoppings do not depend on the normalization of Gamma
    EXPECT_NEAR(static_cast<double>(split.T[n](1, 1)), static_cast<double>(split.T[n](0, 0)),
                1e-14 * static_cast<double>(split.T[n](0, 0)));

  // Kept whole, the weak channel is taken for a combination that does not couple.
  const auto whole = build_chain<Real>(star_of<double>(gamma, 40, false), chain_options(8));
  EXPECT_EQ(whole.diagnostics.theta_rank, 1);
  EXPECT_EQ(whole.V(1, 1), Real(0));
  for (unsigned int n = 0; n <= whole.Nmax; n++) EXPECT_LT(static_cast<double>(abs(whole.T[n](1, 1))), 1e-40);
}

TEST(MixChainChain, blocks_are_placed_in_their_channels_and_the_map_stays_unitary) { // NOLINT
  // Channels 1 and 3 are coupled, channel 2 is on its own.
  const auto gamma = [](const double omega) {
    Matrix<std::complex<double>> m = Matrix<std::complex<double>>::Zero(3, 3);
    m(0, 0)                        = 0.5 + 0.2 * omega;
    m(1, 1)                        = 0.3 + omega * omega;
    m(2, 2)                        = 0.4 - 0.1 * omega;
    m(0, 2)                        = std::complex<double>(0.1 * omega, 0.05);
    m(2, 0)                        = std::conj(m(0, 2));
    return m;
  };
  // mMAX=3 gives 2*3*4 = 24 levels, exactly the 3*(Nmax+1) that a chain with Nmax=7 spans.
  const auto input = input_of<std::complex<double>>(gamma);
  const auto star  = star_from(input, 3);
  ASSERT_EQ(star.blocks, (Blocks{{0, 2}, {1}}));
  std::vector<Matrix<Complex>> lanczos;
  const auto chain = build_chain(to_wide<Complex>(star), chain_options(6), &lanczos);
  EXPECT_EQ(chain.blocks, star.blocks);

  const auto outer = build_chain<Complex>(star_from(restrict_input(input, Block{0, 2}), 3), chain_options(6));
  const auto inner = build_chain<Complex>(star_from(restrict_input(input, Block{1}), 3), chain_options(6));
  const auto check = [](const Matrix<Complex> &whole, const Matrix<Complex> &outer_part,
                        const Matrix<Complex> &inner_part) {
    const int outer_channel[] = {0, 2};
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) EXPECT_EQ(whole(outer_channel[i], outer_channel[j]), outer_part(i, j));
    EXPECT_EQ(whole(1, 1), inner_part(0, 0));
    for (const int i : {0, 2}) {
      EXPECT_EQ(whole(i, 1), Complex(0));
      EXPECT_EQ(whole(1, i), Complex(0));
    }
  };
  check(chain.V, outer.V, inner.V);
  for (unsigned int n = 0; n <= chain.Nmax; n++) check(chain.E[n], outer.E[n], inner.E[n]);
  for (unsigned int n = 0; n <= chain.Nmax; n++) check(chain.T[n], outer.T[n], inner.T[n]);

  // With Nmax=6 the blocks Q_0..Q_7 of the two parts fill the 24 levels: stacked side by side they are unitary, and
  // they satisfy the recursion of the whole star, blocks and all.
  Matrix<Complex> q(24, 24);
  ASSERT_EQ(lanczos.size(), chain.Nmax + 2);
  for (std::size_t s = 0; s < lanczos.size(); s++) q.block(0, static_cast<Eigen::Index>(3 * s), 24, 3) = lanczos[s];
  EXPECT_LT(largest<Complex>(q.adjoint() * q - Matrix<Complex>::Identity(24, 24)), 1e-40);
  Matrix<Complex> star_hamiltonian = Matrix<Complex>::Zero(24, 24);
  for (std::size_t k = 0; k < star.levels.size(); k++)
    star_hamiltonian(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(k)) =
      make_scalar<Complex>(star.levels[k].energy, 0);
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    Matrix<Complex> recursion = lanczos[n] * chain.E[n] + lanczos[n + 1] * chain.T[n];
    if (n > 0) recursion += lanczos[n - 1] * chain.T[n - 1].adjoint();
    EXPECT_LT(largest<Complex>(star_hamiltonian * lanczos[n] - recursion), 1e-40) << "site " << n;
  }
}

TEST(MixChainChain, a_zero_block_has_a_zero_chain_and_does_not_spoil_the_diagnostics) { // NOLINT
  const auto chain =
    build_chain<Real>(star_of<double>([](const double omega) { return diagonal_of(0.4 + omega, 0.0); }, 40),
                      chain_options(8));
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
  EXPECT_EQ(chain.V(1, 1), Real(0));
  for (unsigned int n = 0; n <= chain.Nmax; n++) EXPECT_EQ(chain.E[n](1, 1), Real(0));
  for (unsigned int n = 0; n <= chain.Nmax; n++) EXPECT_EQ(chain.T[n](1, 1), Real(0));
}

TEST(MixChainChain, ranks_of_blocks_add_up_site_by_site) { // NOLINT
  // Block {1} has 12 coupled levels; block {2} has 2, padded with 10 levels without coupling, so that it passes the
  // count of levels but its Krylov space runs out after two sites.
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

  const auto chain = build_chain<Real>(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 2);
  EXPECT_EQ(chain.diagnostics.hopping_ranks, (std::vector<int>{2, 1, 1, 1, 1})); // one hopping per site, 0..Nmax
  EXPECT_EQ(chain.diagnostics.min_rank, 1);
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 1U);
  EXPECT_FALSE(chain.block_diagnostics[0].rank_drop_site.has_value());
  ASSERT_TRUE(chain.block_diagnostics[1].rank_drop_site.has_value());
  EXPECT_EQ(*chain.block_diagnostics[1].rank_drop_site, 1U);
  EXPECT_EQ(chain.block_diagnostics[1].min_rank, 0);
  EXPECT_EQ(chain.block_diagnostics[0].levels, 12);
  EXPECT_EQ(chain.block_diagnostics[0].coupled_levels, 12);
  EXPECT_EQ(chain.block_diagnostics[1].levels, 12);
  EXPECT_EQ(chain.block_diagnostics[1].coupled_levels, 2);
  EXPECT_EQ(chain.diagnostics.coupled_levels, 14);

  const auto alone_first = build_chain<Real>(first, chain_options(4));
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_EQ(chain.T[n](0, 0), alone_first.T[n](0, 0)) << "site " << n;
    if (n > 0) {
      EXPECT_EQ(chain.T[n](1, 1), Real(0)) << "site " << n;
    }
  }
}

TEST(MixChainChain, the_nambu_gauge_puts_the_blocks_into_the_nambu_structure) { // NOLINT
  // A Nambu-symmetric bath: the normal part is flat and equal for the particle and the hole, and the anomalous part
  // is odd in omega, as it is for a superconductor. The chain then has the structure that a consumer of xi, zeta,
  // scdelta and sckappa relies on, but only in the nambu gauge: the polar gauge, which makes every T_n positive
  // semidefinite, hides the sign of the hole component in the Lanczos blocks instead.
  const double rho = 0.3, anomalous = 0.1;
  const auto block = [](const double diagonal, const double offdiagonal) {
    Matrix<double> m = Matrix<double>::Zero(2, 2);
    m(0, 0)          = diagonal;
    m(1, 1)          = diagonal;
    m(0, 1)          = offdiagonal;
    m(1, 0)          = offdiagonal;
    return m;
  };
  GammaInput<double> input;
  input.channels = 2;
  for (int k = 0; k <= 100; k++) {
    input.pos.omega.push_back(0.01 * k);
    input.pos.gamma.push_back(block(rho, anomalous));
    input.neg.omega.push_back(0.01 * k);
    input.neg.gamma.push_back(block(rho, -anomalous)); // odd in omega
  }
  const auto star = star_from(input, 40);
  ASSERT_EQ(star.blocks.size(), 1U); // particle and hole are coupled, hence one block

  auto options     = chain_options(8);
  const auto polar = build_chain<Real>(star, options);
  options.gauge    = ChainGauge::nambu;
  const auto nambu = build_chain<Real>(star, options);
  EXPECT_EQ(polar.gauge, ChainGauge::polar);
  EXPECT_EQ(nambu.gauge, ChainGauge::nambu);
  EXPECT_LT(nambu.diagnostics.max_nambu_deviation, 1e-30);

  // V(2,2) = -V(1,1): the impurity index of V is physical, so the hole row must carry the Nambu sign, which is why
  // U_0 = diag(1,-1) rather than the identity.
  EXPECT_LT(static_cast<double>(abs(nambu.V(1, 1) + nambu.V(0, 0))), 1e-30);
  EXPECT_EQ(nambu.V(0, 0), polar.V(0, 0));
  EXPECT_EQ(nambu.V(1, 1), -polar.V(1, 1));
  EXPECT_EQ(nambu.V(0, 1), -polar.V(0, 1));

  for (unsigned int n = 0; n <= nambu.Nmax; n++) {
    // E(2,2) = -E(1,1) and T(2,2) = -T(1,1), which is what the four stored numbers per site rely on.
    EXPECT_LT(static_cast<double>(abs(nambu.E[n](1, 1) + nambu.E[n](0, 0))), 1e-30) << "site " << n;
    EXPECT_LT(static_cast<double>(abs(nambu.T[n](1, 1) + nambu.T[n](0, 0))), 1e-30) << "site " << n;
    // The gauge only flips signs: the hole component of the even sites.
    const auto sign = Real(n % 2 == 0 ? -1 : 1);
    EXPECT_EQ(nambu.E[n](0, 1), sign * polar.E[n](0, 1)) << "site " << n;
    EXPECT_EQ(nambu.E[n](0, 0), polar.E[n](0, 0)) << "site " << n;
    EXPECT_EQ(nambu.T[n](1, 1), -polar.T[n](1, 1)) << "site " << n;
    EXPECT_EQ(nambu.T[n](0, 0), polar.T[n](0, 0)) << "site " << n;
  }
}

TEST(MixChainChain, the_nambu_gauge_refuses_a_chain_without_the_structure) { // NOLINT
  auto options  = chain_options(3);
  options.gauge = ChainGauge::nambu;
  // Two independent channels: every block has one channel, so there is no particle-hole pair to flip.
  EXPECT_THROW(build_chain<Real>(star_of<double>([](const double omega) { return diagonal_of(0.4 + omega, 0.2); }, 40),
                                 options),
               std::invalid_argument);
  // A generic 2x2 star is one block, but its chain has no Nambu structure.
  EXPECT_THROW(build_chain<Real>(arbitrary_star<double>(2, 12), options), std::runtime_error);
}

TEST(MixChainChain, rejects_a_star_too_small_for_the_chain) { // NOLINT
  // A chain of 4 sites with 2 channels needs 2*(4+1) = 10 levels: one block per site and one for the hopping out of
  // the last one.
  const auto star = arbitrary_star<double>(2, 9);
  EXPECT_THROW(build_chain<Real>(star, chain_options(3)), std::invalid_argument);
  EXPECT_NO_THROW(build_chain<Real>(arbitrary_star<double>(2, 10), chain_options(3)));
  EXPECT_THROW(build_chain<Real>(arbitrary_star<double>(2, 12), chain_options(0)), std::invalid_argument);
}

TEST(MixChainChain, works_in_double_precision_too) { // NOLINT
  // The unit tests mostly use 50 digits; the double instantiation must give the same chain to rounding.
  const auto star   = arbitrary_star<double>(2, 12);
  const auto narrow = build_chain<double>(star, chain_options(3));
  const auto wide   = convert_chain<double>(build_chain<Real>(star, chain_options(3)));
  EXPECT_LT(largest<double>(narrow.V - wide.V), 1e-14);
  for (unsigned int n = 0; n <= 3; n++) EXPECT_LT(largest<double>(narrow.T[n] - wide.T[n]), 1e-13);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
