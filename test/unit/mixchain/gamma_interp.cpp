#include <gtest/gtest.h>

#include <complex>
#include <functional>
#include <vector>

#include <mixchain/gamma_interp.hpp>

using namespace NRG::MixChain;

namespace {

using Complex = std::complex<double>;

const std::vector<double> grid{0.0, 0.25, 0.5, 0.75, 1.0};

template<typename S>
GammaBranch<S> make_branch(const std::function<Matrix<S>(double)> &gamma) {
  GammaBranch<S> branch;
  branch.omega = grid;
  for (const auto omega : grid) branch.gamma.push_back(gamma(omega));
  return branch;
}

// A Gamma linear in omega, so that the linear interpolation is exact everywhere.
Matrix<double> real_gamma(const double omega) {
  Matrix<double> m = Matrix<double>::Zero(2, 2);
  m(0, 0)          = 1.0 + omega;
  m(1, 1)          = 2.0 - 0.5 * omega;
  m(0, 1)          = 0.5 - omega;
  m(1, 0)          = 0.5 - omega;
  return m;
}

Matrix<Complex> complex_gamma(const double omega) {
  Matrix<Complex> m = Matrix<Complex>::Zero(2, 2);
  m(0, 0)           = 1.0 + omega;
  m(1, 1)           = 2.0 - 0.5 * omega;
  m(0, 1)           = Complex(0.5 - omega, 0.3 * omega);
  m(1, 0)           = Complex(0.5 - omega, -0.3 * omega);
  return m;
}

template<typename S> double distance(const Matrix<S> &a, const Matrix<S> &b) {
  return (a - b).cwiseAbs().maxCoeff();
}

} // namespace

TEST(MixChainGammaInterp, reproduces_the_input_at_the_nodes) { // NOLINT
  const auto branch = make_branch<double>(real_gamma);
  for (const auto method : {NRG::Tools::InterpolationMethod::linear, NRG::Tools::InterpolationMethod::steffen}) {
    GammaInterpolation<double> interpolation(branch, method);
    EXPECT_EQ(interpolation.channels(), 2);
    for (std::size_t k = 0; k < grid.size(); k++) EXPECT_LT(distance(interpolation(grid[k]), branch.gamma[k]), 1e-14);
  }
}

TEST(MixChainGammaInterp, reproduces_a_complex_input_at_the_nodes) { // NOLINT
  const auto branch = make_branch<Complex>(complex_gamma);
  for (const auto method : {NRG::Tools::InterpolationMethod::linear, NRG::Tools::InterpolationMethod::steffen}) {
    GammaInterpolation<Complex> interpolation(branch, method);
    for (std::size_t k = 0; k < grid.size(); k++) EXPECT_LT(distance(interpolation(grid[k]), branch.gamma[k]), 1e-14);
  }
}

TEST(MixChainGammaInterp, is_hermitian_between_the_nodes) { // NOLINT
  const auto branch = make_branch<Complex>(complex_gamma);
  for (const auto method : {NRG::Tools::InterpolationMethod::linear, NRG::Tools::InterpolationMethod::steffen}) {
    GammaInterpolation<Complex> interpolation(branch, method);
    for (const double omega : {0.1, 0.37, 0.6, 0.99}) {
      const auto m = interpolation(omega);
      // Only the upper triangle is interpolated, so this holds exactly rather than to rounding.
      EXPECT_EQ(m(0, 1), std::conj(m(1, 0)));
      EXPECT_EQ(m(0, 0).imag(), 0.0);
      EXPECT_EQ(m(1, 1).imag(), 0.0);
    }
  }
}

TEST(MixChainGammaInterp, linear_interpolation_is_exact_for_a_linear_gamma) { // NOLINT
  const auto branch = make_branch<double>(real_gamma);
  GammaInterpolation<double> interpolation(branch, NRG::Tools::InterpolationMethod::linear);
  for (const double omega : {0.125, 0.375, 0.625, 0.875}) {
    EXPECT_LT(distance(interpolation(omega), real_gamma(omega)), 1e-15);
  }
}

TEST(MixChainGammaInterp, interpolates_the_off_diagonal_imaginary_part) { // NOLINT
  const auto branch = make_branch<Complex>(complex_gamma);
  GammaInterpolation<Complex> interpolation(branch, NRG::Tools::InterpolationMethod::linear);
  const auto m = interpolation(0.125); // midway between the first two nodes
  EXPECT_NEAR(m(0, 1).imag(), 0.3 * 0.125, 1e-15);
  EXPECT_NEAR(m(0, 1).real(), 0.5 - 0.125, 1e-15);
}

TEST(MixChainGammaInterp, holds_the_values_constant_beyond_the_ends) { // NOLINT
  const auto branch = make_branch<Complex>(complex_gamma);
  for (const auto method : {NRG::Tools::InterpolationMethod::linear, NRG::Tools::InterpolationMethod::steffen}) {
    GammaInterpolation<Complex> interpolation(branch, method);
    // PiecewisePolynomial::evaluate() would throw outside its domain; the argument is clamped instead.
    EXPECT_LT(distance(interpolation(-0.5), branch.gamma.front()), 1e-14);
    EXPECT_LT(distance(interpolation(1.5), branch.gamma.back()), 1e-14);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
