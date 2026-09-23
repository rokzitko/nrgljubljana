#include <gtest/gtest.h>

#include <complex>
#include <stdexcept>
#include <type_traits>
#include <string>

#include <Eigen/Dense>

#include <mixchain/precision.hpp>
#include <mixchain/types.hpp>

using namespace NRG::MixChain;

namespace {

// The Hermitian square root of [[2, c], [conj(c), 3]], which is what the polar gauge of the block Lanczos needs, and
// the deepest thing the multiprecision types have to support.
template<typename S> double square_root_residual(const S &offdiagonal) {
  using Matrix = Eigen::Matrix<S, -1, -1>;
  Matrix a(2, 2);
  a(0, 0) = S(2);
  a(1, 1) = S(3);
  a(0, 1) = offdiagonal;
  a(1, 0) = Eigen::numext::conj(offdiagonal);

  Eigen::SelfAdjointEigenSolver<Matrix> solver(a);
  EXPECT_EQ(solver.info(), Eigen::Success);
  const Matrix root = solver.eigenvectors()
                      * solver.eigenvalues().cwiseSqrt().template cast<S>().asDiagonal()
                      * solver.eigenvectors().adjoint();
  return static_cast<double>((root * root - a).cwiseAbs().maxCoeff());
}

} // namespace

TEST(MixChainPrecision, resolves_a_request_to_the_smallest_rung_that_covers_it) { // NOLINT
  EXPECT_EQ(resolve_precision(100), 50U);  // 31 digits
  EXPECT_EQ(resolve_precision(166), 50U);  // 50 digits exactly
  EXPECT_EQ(resolve_precision(167), 200U); // just past the first rung
  EXPECT_EQ(resolve_precision(600), 200U); // 181 digits
  EXPECT_EQ(resolve_precision(2000), 800U); // the default preccpp of nrgchain, 603 digits
}

TEST(MixChainPrecision, converts_between_bits_and_digits) { // NOLINT
  EXPECT_EQ(digits_for_bits(2000), 603U);
  EXPECT_GE(bits_for_digits(800), 2657U);
  // Every rung must be reachable by some request, or it would be dead code.
  for (const auto rung : precision_ladder) EXPECT_EQ(resolve_precision(bits_for_digits(rung)), rung);
}

TEST(MixChainPrecision, rejects_a_request_it_cannot_serve) { // NOLINT
  EXPECT_THROW(resolve_precision(10), std::invalid_argument);  // the guard of nrgchain
  EXPECT_THROW(resolve_precision(0), std::invalid_argument);
  EXPECT_THROW(resolve_precision(4000), std::invalid_argument); // beyond the top rung
  try {
    resolve_precision(4000);
  } catch (const std::invalid_argument &error) {
    const std::string message = error.what();
    EXPECT_NE(message.find("800"), std::string::npos); // it says what the ceiling is
  }
}

TEST(MixChainPrecision, dispatches_to_a_real_and_to_a_complex_type) { // NOLINT
  for (const bool complex_data : {false, true}) {
    const auto is_complex = with_precision(2000, complex_data, []<typename S>() {
      return static_cast<bool>(Eigen::NumTraits<S>::IsComplex);
    });
    EXPECT_EQ(is_complex, complex_data);
  }
}

TEST(MixChainPrecision, carries_far_more_precision_than_double) { // NOLINT
  // The dispatch compiles both branches whatever the flag says at run time, so the body must be valid for a real and
  // for a complex scalar alike: the off-diagonal element is assembled with make_scalar rather than written as S(re, im).
  const auto residual =
    with_precision(100, false, []<typename S>() { return square_root_residual<S>(make_scalar<S>(1, 0)); });
  EXPECT_LT(residual, 1e-40); // in double precision this would be around 1e-16

  const auto complex_residual =
    with_precision(100, true, []<typename S>() { return square_root_residual<S>(make_scalar<S>(0, 1)); });
  EXPECT_LT(complex_residual, 1e-40);

  // The top rung is far beyond the smallest chain coefficient of any reasonable run, Lambda^(-Nmax/2).
  const auto deep =
    with_precision(2000, true, []<typename S>() { return square_root_residual<S>(make_scalar<S>(0.5, 0.5)); });
  EXPECT_LT(deep, 1e-300);
}

TEST(MixChainPrecision, the_scalar_traits_cover_the_multiprecision_types) { // NOLINT
  static_assert(!is_complex_v<WideReal<50>>);
  static_assert(is_complex_v<WideComplex<50>>);
  static_assert(!is_complex_v<double>);
  static_assert(is_complex_v<std::complex<double>>);
  static_assert(std::is_same_v<real_type<WideComplex<50>>, WideReal<50>>);
  static_assert(std::is_same_v<real_type<double>, double>);
  EXPECT_EQ(make_scalar<WideComplex<50>>(1, 2).imag(), 2);
  EXPECT_EQ(make_scalar<WideReal<50>>(1, 2), 1); // the imaginary part is dropped, as for any real scalar
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
