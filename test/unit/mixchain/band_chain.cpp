#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <star-to-chain.hpp>

#include <mixchain/band_chain.hpp>

using namespace NRG::MixChain;

namespace {

// The levels of a star: energies falling off as Lambda^(-m), of both signs, with couplings of matching size whose
// direction turns from level to level. Complex couplings when S is complex.
template<typename S> struct Levels {
  std::vector<double> energies;
  Matrix<S> start;
};

template<typename S> Levels<S> graded_levels(const int channels, const int count, const double lambda = 2.0) {
  Levels<S> levels;
  levels.start = Matrix<S>(count, channels);
  for (int k = 0; k < count; k++) {
    const auto m = k / 2; // two levels per interval, one of each sign
    levels.energies.push_back((k % 2 ? -1.0 : 1.0) * std::pow(lambda, -m) * (0.7 + 0.1 * std::sin(1.7 * k)));
    for (int i = 0; i < channels; i++)
      levels.start(k, i) = make_scalar<S>(std::cos(0.9 * k + 1.3 * i) + 0.2, 0.6 * std::sin(0.5 * k - i))
                           * std::pow(lambda, -0.5 * m) * (k % 2 ? 0.8 : 1.0);
  }
  return levels;
}

template<typename S> double largest(const Matrix<S> &m) { return m.size() ? m.cwiseAbs().maxCoeff() : 0.0; }

// The kept part of the band as one matrix: E_n on the diagonal, T_n below it and T_n^dag above.
template<typename S> Matrix<S> band_matrix(const BandChain<S> &chain) {
  const auto p      = chain.R.rows();
  const auto blocks = static_cast<Eigen::Index>(chain.E.size());
  Matrix<S> h       = Matrix<S>::Zero(p * blocks, p * blocks);
  for (Eigen::Index s = 0; s < blocks; s++) {
    h.block(s * p, s * p, p, p) = chain.E[static_cast<std::size_t>(s)];
    if (s + 1 < blocks) {
      h.block((s + 1) * p, s * p, p, p) = chain.T[static_cast<std::size_t>(s)];
      h.block(s * p, (s + 1) * p, p, p) = chain.T[static_cast<std::size_t>(s)].adjoint();
    }
  }
  return h;
}

} // namespace

TEST(MixChainBand, one_channel_is_the_scalar_chain) { // NOLINT
  const auto levels = graded_levels<double>(1, 100);
  const std::size_t blocks = 31;
  const auto band = band_star_to_chain(levels.energies, levels.start, blocks);

  std::vector<NRG::StarPoint> points;
  double weight = 0.0;
  for (std::size_t k = 0; k < levels.energies.size(); k++) {
    const auto amplitude = std::abs(levels.start(static_cast<Eigen::Index>(k), 0));
    points.push_back({levels.energies[k], amplitude});
    weight += amplitude * amplitude;
  }
  const auto scalar = NRG::scalar_star_to_chain(points, blocks - 1);

  EXPECT_EQ(band.rows, blocks);
  EXPECT_NEAR(band.R(0, 0), std::sqrt(weight), 1e-14 * std::sqrt(weight));
  for (std::size_t n = 0; n + 1 < blocks; n++) {
    const auto xi = scalar.xi[n];
    EXPECT_NEAR(band.T[n](0, 0), xi, 1e-13 * xi) << "site " << n;
    EXPECT_NEAR(band.E[n](0, 0), scalar.zeta[n], 1e-13 * std::max(xi, n > 0 ? scalar.xi[n - 1] : 0.0)) << "site " << n;
  }
}

TEST(MixChainBand, the_kept_blocks_are_those_of_the_full_reduction) { // NOLINT
  for (const bool complex_data : {false, true}) {
    const auto check = [&]<typename S>() {
      const auto levels = graded_levels<S>(2, 40);
      const auto full   = band_star_to_chain(levels.energies, levels.start, 20);
      const auto kept   = band_star_to_chain(levels.energies, levels.start, 5);
      EXPECT_EQ(full.rows, 40U);
      EXPECT_EQ(kept.rows, 10U);
      EXPECT_LT(largest<S>(kept.R - full.R), 1e-14);
      for (std::size_t n = 0; n < 5; n++) EXPECT_LT(largest<S>(kept.E[n] - full.E[n]), 1e-14) << "site " << n;
      for (std::size_t n = 0; n < 4; n++) EXPECT_LT(largest<S>(kept.T[n] - full.T[n]), 1e-14) << "site " << n;
    };
    if (complex_data)
      check.template operator()<std::complex<double>>();
    else
      check.template operator()<double>();
  }
}

TEST(MixChainBand, the_full_reduction_is_a_unitary_transformation_of_the_star) { // NOLINT
  // Kept whole, the band matrix has the spectrum of the star, and R^dag R = A^dag A.
  const auto levels = graded_levels<std::complex<double>>(2, 12);
  const auto band   = band_star_to_chain(levels.energies, levels.start, 6);
  ASSERT_EQ(band.rows, 12U);

  const Eigen::MatrixXcd h = band_matrix(band);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(h, Eigen::EigenvaluesOnly);
  ASSERT_EQ(solver.info(), Eigen::Success);
  auto expected = levels.energies;
  std::sort(expected.begin(), expected.end());
  for (std::size_t k = 0; k < expected.size(); k++)
    EXPECT_NEAR(solver.eigenvalues()(static_cast<Eigen::Index>(k)), expected[k], 1e-13) << "level " << k;

  const Matrix<std::complex<double>> theta = levels.start.adjoint() * levels.start;
  EXPECT_LT(largest<std::complex<double>>(band.R.adjoint() * band.R - theta), 1e-13 * largest(theta));
}

TEST(MixChainBand, moments_of_the_star_are_reproduced) { // NOLINT
  // n+1 blocks reproduce A^dag D^q A = R^dag (H^q)_00 R for q <= 2n+1, also when later levels were only chased
  // through the band.
  for (const bool complex_data : {false, true}) {
    const auto check = [&]<typename S>() {
      const auto levels = graded_levels<S>(2, 24);
      const auto band   = band_star_to_chain(levels.energies, levels.start, 4);
      const auto h      = band_matrix(band);
      Matrix<S> power   = Matrix<S>::Identity(h.rows(), h.cols());
      Vector<S> weights = Vector<S>::Ones(levels.start.rows());
      for (int q = 0; q <= 7; q++) {
        const Matrix<S> exact      = levels.start.adjoint() * weights.asDiagonal() * levels.start;
        const Matrix<S> from_chain = band.R.adjoint() * power.topLeftCorner(2, 2) * band.R;
        EXPECT_LT(largest<S>(from_chain - exact), 1e-12 * largest(exact)) << "moment " << q;
        power = (power * h).eval();
        for (Eigen::Index k = 0; k < weights.size(); k++)
          weights(k) *= make_scalar<S>(levels.energies[static_cast<std::size_t>(k)], 0);
      }
    };
    if (complex_data)
      check.template operator()<std::complex<double>>();
    else
      check.template operator()<double>();
  }
}

TEST(MixChainBand, the_blocks_are_triangular_with_a_real_nonnegative_diagonal) { // NOLINT
  // More levels than rows of the band, so that every kept row has served as a pivot.
  const auto levels = graded_levels<std::complex<double>>(3, 60);
  const auto band   = band_star_to_chain(levels.energies, levels.start, 6);
  const auto check  = [](const Matrix<std::complex<double>> &m) {
    for (Eigen::Index i = 0; i < m.rows(); i++) {
      for (Eigen::Index j = 0; j < i; j++) EXPECT_EQ(m(i, j), std::complex<double>(0.0, 0.0));
      EXPECT_EQ(m(i, i).imag(), 0.0);
      EXPECT_GT(m(i, i).real(), 0.0);
    }
  };
  check(band.R);
  for (const auto &hopping : band.T) check(hopping);
  for (const auto &onsite : band.E) EXPECT_EQ(largest<std::complex<double>>(onsite - onsite.adjoint()), 0.0);
}

TEST(MixChainBand, real_data_in_complex_arithmetic_stays_real) { // NOLINT
  const auto real_levels = graded_levels<double>(2, 40);
  const Matrix<std::complex<double>> start = real_levels.start.cast<std::complex<double>>();
  const auto band      = band_star_to_chain(real_levels.energies, start, 6);
  const auto reference = band_star_to_chain(real_levels.energies, real_levels.start, 6);
  const auto check     = [](const Matrix<std::complex<double>> &m, const Matrix<double> &expected) {
    for (Eigen::Index i = 0; i < m.rows(); i++)
      for (Eigen::Index j = 0; j < m.cols(); j++) {
        EXPECT_EQ(m(i, j).imag(), 0.0);
        EXPECT_NEAR(m(i, j).real(), expected(i, j), 1e-14);
      }
  };
  check(band.R, reference.R);
  for (std::size_t n = 0; n < band.E.size(); n++) check(band.E[n], reference.E[n]);
  for (std::size_t n = 0; n < band.T.size(); n++) check(band.T[n], reference.T[n]);
}

TEST(MixChainBand, a_diagonal_star_gives_two_scalar_chains_and_exact_zeros_between_them) { // NOLINT
  // Every level couples to one channel only, the channels alternating: the zeros are never touched by a rotation.
  const auto first  = graded_levels<double>(1, 30);
  const auto second = graded_levels<double>(1, 30, 3.0);
  std::vector<double> energies;
  Matrix<double> start = Matrix<double>::Zero(60, 2);
  for (int k = 0; k < 30; k++) {
    energies.push_back(first.energies[static_cast<std::size_t>(k)]);
    start(2 * k, 0) = first.start(k, 0);
    energies.push_back(second.energies[static_cast<std::size_t>(k)]);
    start(2 * k + 1, 1) = second.start(k, 0);
  }
  const auto band = band_star_to_chain(energies, start, 8);
  const std::vector<BandChain<double>> alone{band_star_to_chain(first.energies, first.start, 8),
                                             band_star_to_chain(second.energies, second.start, 8)};
  const auto check = [&](const Matrix<double> &joint, const auto &part_of) {
    EXPECT_EQ(joint(0, 1), 0.0);
    EXPECT_EQ(joint(1, 0), 0.0);
    for (int channel = 0; channel < 2; channel++) {
      const auto expected = part_of(alone[static_cast<std::size_t>(channel)]);
      EXPECT_NEAR(joint(channel, channel), expected, 1e-14 * std::max(std::abs(expected), 1e-3));
    }
  };
  check(band.R, [](const BandChain<double> &part) { return part.R(0, 0); });
  for (std::size_t n = 0; n < band.E.size(); n++)
    check(band.E[n], [n](const BandChain<double> &part) { return part.E[n](0, 0); });
  for (std::size_t n = 0; n < band.T.size(); n++)
    check(band.T[n], [n](const BandChain<double> &part) { return part.T[n](0, 0); });
}

TEST(MixChainBand, scales_with_the_energies_and_with_the_couplings) { // NOLINT
  const auto levels    = graded_levels<std::complex<double>>(2, 40);
  const auto reference = band_star_to_chain(levels.energies, levels.start, 6);
  using M = Matrix<std::complex<double>>;
  for (const double energy_scale : {1e-100, 1.0, 1e100}) {
    for (const double coupling_scale : {1e-200, 1.0, 1e200}) {
      auto energies = levels.energies;
      for (auto &energy : energies) energy *= energy_scale;
      const M start   = levels.start * coupling_scale;
      const auto band = band_star_to_chain(energies, start, 6);
      EXPECT_LT(largest<std::complex<double>>(M(band.R / coupling_scale) - reference.R), 1e-13);
      for (std::size_t n = 0; n < band.E.size(); n++)
        EXPECT_LT(largest<std::complex<double>>(M(band.E[n] / energy_scale) - reference.E[n]), 1e-13) << "site " << n;
      for (std::size_t n = 0; n < band.T.size(); n++)
        EXPECT_LT(largest<std::complex<double>>(M(band.T[n] / energy_scale) - reference.T[n]),
                  1e-13 * largest(reference.T[n]))
          << "site " << n;
    }
  }
}

TEST(MixChainBand, a_short_star_fills_only_part_of_the_band) { // NOLINT
  const auto levels = graded_levels<double>(2, 5);
  const auto band   = band_star_to_chain(levels.energies, levels.start, 4);
  EXPECT_EQ(band.rows, 5U);
  ASSERT_EQ(band.E.size(), 4U);
  ASSERT_EQ(band.T.size(), 3U);
  // Rows 0..4 hold the star: two full sites and the first orbital of the third.
  EXPECT_GT(std::abs(band.E[2](0, 0)), 0.0);
  EXPECT_EQ(band.E[2](1, 1), 0.0);
  EXPECT_EQ(band.T[1].row(1).cwiseAbs().maxCoeff(), 0.0);
  EXPECT_EQ(largest(band.E[3]), 0.0);
  EXPECT_EQ(largest(band.T[2]), 0.0);
  // What is there is still the whole star.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Eigen::MatrixXd(band_matrix(band).topLeftCorner(5, 5)),
                                                        Eigen::EigenvaluesOnly);
  auto expected = levels.energies;
  std::sort(expected.begin(), expected.end());
  for (std::size_t k = 0; k < 5; k++) EXPECT_NEAR(solver.eigenvalues()(static_cast<Eigen::Index>(k)), expected[k], 1e-14);
}

TEST(MixChainBand, rejects_invalid_input) { // NOLINT
  const auto levels = graded_levels<double>(2, 10);
  EXPECT_THROW(band_star_to_chain(levels.energies, levels.start, 0), std::invalid_argument);
  auto short_energies = levels.energies;
  short_energies.pop_back();
  EXPECT_THROW(band_star_to_chain(short_energies, levels.start, 3), std::invalid_argument);
  EXPECT_THROW(band_star_to_chain(levels.energies, Matrix<double>(10, 0), 3), std::invalid_argument);
  auto energies = levels.energies;
  energies[3]   = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(band_star_to_chain(energies, levels.start, 3), std::invalid_argument);
  auto start  = levels.start;
  start(4, 1) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(band_star_to_chain(levels.energies, start, 3), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
