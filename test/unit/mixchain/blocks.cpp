#include <gtest/gtest.h>

#include <complex>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <mixchain/blocks.hpp>

using namespace NRG::MixChain;

namespace {

template<typename S> GammaBranch<S> make_branch(const std::function<Matrix<S>(double)> &gamma) {
  GammaBranch<S> branch;
  for (const auto omega : {1e-99, 0.1, 0.5, 1.0}) {
    branch.omega.push_back(omega);
    branch.gamma.push_back(gamma(omega));
  }
  branch.innermost = 0.1;
  return branch;
}

template<typename S>
GammaInput<S> make_input(const int channels, const std::function<Matrix<S>(double)> &pos,
                         const std::function<Matrix<S>(double)> &neg) {
  GammaInput<S> input;
  input.channels = channels;
  input.pos      = make_branch(pos);
  input.neg      = make_branch(neg);
  return input;
}

template<typename S>
GammaInput<S> make_input(const int channels, const std::function<Matrix<S>(double)> &gamma) {
  return make_input(channels, gamma, gamma);
}

// A positive diagonal, with the given off-diagonal element (i,j) and its conjugate set to 'value' at omega=0.5 only.
template<typename S>
auto with_elements(const int channels, const std::vector<std::pair<int, int>> &elements, const S value) {
  return [=](const double omega) {
    Matrix<S> m = Matrix<S>::Identity(channels, channels);
    if (omega == 0.5)
      for (const auto &[i, j] : elements) {
        m(i, j) = value;
        m(j, i) = Eigen::numext::conj(value);
      }
    return m;
  };
}

} // namespace

TEST(MixChainBlocks, diagonal_gamma_has_one_block_per_channel) { // NOLINT
  const auto input = make_input<double>(3, [](const double) { return Matrix<double>(Matrix<double>::Identity(3, 3)); });
  EXPECT_EQ(gamma_blocks(input), (Blocks{{0}, {1}, {2}}));
}

TEST(MixChainBlocks, blocks_need_not_be_contiguous) { // NOLINT
  const auto input = make_input<double>(4, with_elements<double>(4, {{0, 2}, {1, 3}}, 0.3));
  EXPECT_EQ(gamma_blocks(input), (Blocks{{0, 2}, {1, 3}}));
}

TEST(MixChainBlocks, links_are_transitive) { // NOLINT
  const auto input = make_input<double>(4, with_elements<double>(4, {{0, 1}, {1, 2}}, 0.3));
  EXPECT_EQ(gamma_blocks(input), (Blocks{{0, 1, 2}, {3}}));
}

TEST(MixChainBlocks, a_link_on_one_frequency_branch_is_enough) { // NOLINT
  const auto diagonal = [](const double) { return Matrix<double>(Matrix<double>::Identity(2, 2)); };
  const auto input    = make_input<double>(2, diagonal, with_elements<double>(2, {{0, 1}}, 0.3));
  EXPECT_EQ(gamma_blocks(input), (Blocks{{0, 1}}));
}

TEST(MixChainBlocks, only_exact_zeros_separate_blocks) { // NOLINT
  const auto tiny = make_input<double>(2, with_elements<double>(2, {{0, 1}}, 1e-300));
  EXPECT_EQ(gamma_blocks(tiny), (Blocks{{0, 1}}));
  // A purely imaginary element links as well.
  using Complex         = std::complex<double>;
  const auto imaginary = make_input<Complex>(2, with_elements<Complex>(2, {{0, 1}}, Complex(0.0, 0.2)));
  EXPECT_EQ(gamma_blocks(imaginary), (Blocks{{0, 1}}));
}

TEST(MixChainBlocks, restrict_input_takes_the_submatrix) { // NOLINT
  const auto gamma = [](const double omega) {
    Matrix<std::complex<double>> m(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) m(i, j) = std::complex<double>(10.0 * i + j, omega);
    return m;
  };
  const auto neg      = [&](const double omega) { return Matrix<std::complex<double>>(2.0 * gamma(omega)); };
  const auto input    = make_input<std::complex<double>>(3, gamma, neg);
  const auto restricted = restrict_input(input, Block{0, 2});

  EXPECT_EQ(restricted.channels, 2);
  ASSERT_EQ(restricted.pos.size(), input.pos.size());
  EXPECT_EQ(restricted.pos.omega, input.pos.omega);
  EXPECT_EQ(restricted.pos.innermost, input.pos.innermost);
  for (std::size_t k = 0; k < input.pos.size(); k++) {
    const auto omega = input.pos.omega[k];
    EXPECT_EQ(restricted.pos.gamma[k](0, 0), std::complex<double>(0.0, omega));
    EXPECT_EQ(restricted.pos.gamma[k](0, 1), std::complex<double>(2.0, omega));
    EXPECT_EQ(restricted.pos.gamma[k](1, 0), std::complex<double>(20.0, omega));
    EXPECT_EQ(restricted.pos.gamma[k](1, 1), std::complex<double>(22.0, omega));
    EXPECT_EQ(restricted.neg.gamma[k](1, 0), std::complex<double>(40.0, 2.0 * omega));
  }
}

TEST(MixChainBlocks, names_round_trip) { // NOLINT
  const Blocks blocks{{0, 2}, {1}, {3, 4}};
  EXPECT_EQ(blocks_name(blocks), "{1,3} {2} {4,5}");
  EXPECT_EQ(parse_blocks(blocks_name(blocks), 5), blocks);
  // Any order and spacing, returned in the canonical order.
  EXPECT_EQ(parse_blocks("  { 5 , 4 }{2}   {3,1} ", 5), blocks);
}

TEST(MixChainBlocks, parse_rejects_what_is_not_a_partition) { // NOLINT
  EXPECT_THROW(parse_blocks("{1,2} {3}", 4), std::invalid_argument);   // channel 4 missing
  EXPECT_THROW(parse_blocks("{1,2} {2,3}", 3), std::invalid_argument); // channel 2 twice
  EXPECT_THROW(parse_blocks("{1,4}", 3), std::invalid_argument);       // out of range
  EXPECT_THROW(parse_blocks("{0,1}", 2), std::invalid_argument);       // 1-based
  EXPECT_THROW(parse_blocks("{1,2} {}", 2), std::invalid_argument);    // empty block
  EXPECT_THROW(parse_blocks("{1,2", 2), std::invalid_argument);        // unterminated
  EXPECT_THROW(parse_blocks("1,2", 2), std::invalid_argument);         // no braces
  EXPECT_THROW(parse_blocks("", 2), std::invalid_argument);            // nothing
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
