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

// A star built by the star stage from Gamma(omega), the same on both frequency branches.
template<typename S0>
Star<S0> star_of(const std::function<Matrix<S0>(double)> &gamma, const unsigned int mmax) {
  GammaBranch<S0> branch;
  for (int k = 0; k <= 100; k++) {
    branch.omega.push_back(0.01 * k);
    branch.gamma.push_back(gamma(0.01 * k));
  }
  GammaInput<S0> input;
  input.channels = static_cast<int>(branch.gamma.front().rows());
  input.pos      = branch;
  input.neg      = branch;
  StarOptions options;
  options.Lambda = NRG::Tools::LambdaCache(lambda_value);
  options.z      = 1.0;
  options.mMAX   = mmax;
  return build_star(input, options);
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

  for (unsigned int n = 0; n < chain.Nmax; n++) {
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

  for (unsigned int n = 0; n < chain.Nmax; n++) {
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

TEST(MixChainChain, full_length_lanczos_is_a_unitary_map_of_the_star) { // NOLINT
  // With exactly channels*(Nmax+1) levels the blocks span the whole space: stacked side by side they form a unitary Q
  // with Q H_chain Q^dag = H_star. This catches what equality of spectra would not, such as a transposed coupling.
  const auto star = arbitrary_star<std::complex<double>>(2, 8);
  std::vector<Matrix<Complex>> blocks;
  const auto chain = build_chain(to_wide<Complex>(star), chain_options(3), &blocks);

  Matrix<Complex> q(8, 8);
  for (std::size_t s = 0; s < blocks.size(); s++) q.block(0, static_cast<Eigen::Index>(2 * s), 8, 2) = blocks[s];
  EXPECT_LT(largest<Complex>(q.adjoint() * q - Matrix<Complex>::Identity(8, 8)), 1e-45);

  Matrix<Complex> star_hamiltonian = Matrix<Complex>::Zero(8, 8);
  for (std::size_t k = 0; k < star.levels.size(); k++)
    star_hamiltonian(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(k)) =
      make_scalar<Complex>(star.levels[k].energy, 0);
  EXPECT_LT(largest<Complex>(q * chain_hamiltonian(chain) * q.adjoint() - star_hamiltonian), 1e-45);
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
  for (unsigned int n = 0; n < plain.Nmax; n++)
    EXPECT_LT(largest<Real>(rotated.T[n] - U * plain.T[n] * U.transpose()), 1e-11) << "site " << n;
}

TEST(MixChainChain, a_diagonal_gamma_gives_independent_scalar_chains) { // NOLINT
  const auto first  = [](const double omega) { return 0.8 - 0.1 * omega; };
  const auto second = [](const double omega) { return 0.3 + 0.1 * omega; };
  const auto joint  = build_chain<Real>(star_of<double>([&](const double omega) {
                                           Matrix<double> m = Matrix<double>::Zero(2, 2);
                                           m(0, 0)          = first(omega);
                                           m(1, 1)          = second(omega);
                                           return m;
                                         }, 40),
                                       chain_options(8));

  for (const auto &[density, channel] : {std::pair{std::function<double(double)>(first), 0},
                                         std::pair{std::function<double(double)>(second), 1}}) {
    const auto alone = build_chain<Real>(
      star_of<double>([&](const double omega) { return Matrix<double>::Constant(1, 1, density(omega)); }, 40),
      chain_options(8));
    for (unsigned int n = 0; n < joint.Nmax; n++) {
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
    for (unsigned int n = 0; n < chain.Nmax; n++)
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
  const auto alone_second = build_chain<Real>(second, chain_options(1));
  const auto near         = [](const Real &a, const Real &b) { return static_cast<double>(abs(a - b)); };
  for (unsigned int n = 0; n <= chain.Nmax; n++) {
    EXPECT_LT(near(chain.E[n](0, 0), alone_first.E[n](0, 0)), 1e-40) << "site " << n;
    EXPECT_LT(near(chain.E[n](1, 1), n <= 1 ? alone_second.E[n](0, 0) : Real(0)), 1e-40) << "site " << n;
    EXPECT_LT(static_cast<double>(abs(chain.E[n](0, 1))), 1e-40) << "site " << n;
  }
  for (unsigned int n = 0; n < chain.Nmax; n++) {
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
  coupled.levels.resize(3);

  const auto chain = build_chain<Real>(star, chain_options(4));
  EXPECT_EQ(chain.diagnostics.theta_rank, 1);
  EXPECT_EQ(chain.diagnostics.min_rank, 0);
  ASSERT_TRUE(chain.diagnostics.rank_drop_site.has_value());
  EXPECT_EQ(*chain.diagnostics.rank_drop_site, 2U);

  const auto alone = build_chain<Real>(coupled, chain_options(2));
  for (unsigned int n = 0; n <= chain.Nmax; n++)
    EXPECT_LT(static_cast<double>(abs(chain.E[n](0, 0) - (n <= 2 ? alone.E[n](0, 0) : Real(0)))), 1e-40) << "site " << n;
  for (unsigned int n = 0; n < chain.Nmax; n++)
    EXPECT_LT(static_cast<double>(abs(chain.T[n](0, 0) - (n <= 1 ? alone.T[n](0, 0) : Real(0)))), 1e-40) << "site " << n;
}

TEST(MixChainChain, rejects_a_star_too_small_for_the_chain) { // NOLINT
  const auto star = arbitrary_star<double>(2, 7); // a chain of 4 sites with 2 channels needs 8 levels
  EXPECT_THROW(build_chain<Real>(star, chain_options(3)), std::invalid_argument);
  EXPECT_THROW(build_chain<Real>(arbitrary_star<double>(2, 12), chain_options(0)), std::invalid_argument);
}

TEST(MixChainChain, works_in_double_precision_too) { // NOLINT
  // The unit tests mostly use 50 digits; the double instantiation must give the same chain to rounding.
  const auto star   = arbitrary_star<double>(2, 12);
  const auto narrow = build_chain<double>(star, chain_options(3));
  const auto wide   = convert_chain<double>(build_chain<Real>(star, chain_options(3)));
  EXPECT_LT(largest<double>(narrow.V - wide.V), 1e-14);
  for (unsigned int n = 0; n < 3; n++) EXPECT_LT(largest<double>(narrow.T[n] - wide.T[n]), 1e-13);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
