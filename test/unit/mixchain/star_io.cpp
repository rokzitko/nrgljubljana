#include <gtest/gtest.h>

#include <complex>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mixchain/star_io.hpp>

using namespace NRG::MixChain;

namespace {

using Complex = std::complex<double>;

constexpr unsigned int mmax = 4;

template<typename S> GammaInput<S> make_input(const std::function<Matrix<S>(double)> &gamma) {
  GammaBranch<S> branch;
  for (int k = 0; k <= 20; k++) {
    const auto omega = 0.05 * k;
    branch.omega.push_back(omega);
    branch.gamma.push_back(gamma(omega));
  }
  GammaInput<S> input;
  input.channels = static_cast<int>(branch.gamma.front().rows());
  input.pos      = branch;
  input.neg      = branch;
  return input;
}

StarOptions options_for() {
  StarOptions options;
  options.Lambda      = NRG::Tools::LambdaCache(2.0);
  options.z           = 1.0;
  options.mMAX        = mmax;
  options.bandrescale = 2.5;
  return options;
}

// A 2x2 Gamma with a nonzero off-diagonal element, so that the coupling vectors have every component filled.
Matrix<double> real_gamma(const double omega) {
  Matrix<double> m = Matrix<double>::Zero(2, 2);
  m(0, 0)          = 0.8 - 0.1 * omega;
  m(1, 1)          = 0.3 + 0.1 * omega;
  m(0, 1)          = 0.1 * omega;
  m(1, 0)          = 0.1 * omega;
  return m;
}

Matrix<Complex> complex_gamma(const double omega) {
  Matrix<Complex> m = Matrix<Complex>::Zero(2, 2);
  m(0, 0)           = 0.8 - 0.1 * omega;
  m(1, 1)           = 0.3 + 0.1 * omega;
  m(0, 1)           = Complex(0.1 * omega, 0.05 * omega);
  m(1, 0)           = Complex(0.1 * omega, -0.05 * omega);
  return m;
}

void write_file(const std::string &filename, const std::string &contents) {
  std::ofstream file(filename);
  file << contents;
}

// A minimal valid star: one channel, mMAX=1, hence four levels.
const char *minimal_star = "# mixchain star\n"
                           "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0\n"
                           "# m sign a E v1\n"
                           "0 + 0 0.7 0.5\n"
                           "1 + 0 0.35 0.3\n"
                           "0 - 0 -0.7 0.5\n"
                           "1 - 0 -0.35 0.3\n";

template<typename S> void expect_same_star(const Star<S> &saved, const Star<S> &loaded) {
  EXPECT_EQ(loaded.channels, saved.channels);
  EXPECT_EQ(loaded.mMAX, saved.mMAX);
  EXPECT_DOUBLE_EQ(loaded.z, saved.z);
  EXPECT_DOUBLE_EQ(loaded.Lambda, saved.Lambda);
  EXPECT_DOUBLE_EQ(loaded.bandrescale, saved.bandrescale);
  ASSERT_EQ(loaded.levels.size(), saved.levels.size());
  for (std::size_t k = 0; k < saved.levels.size(); k++) {
    EXPECT_EQ(loaded.levels[k].m, saved.levels[k].m);
    EXPECT_EQ(loaded.levels[k].branch, saved.levels[k].branch);
    EXPECT_TRUE(loaded.levels[k].sign == saved.levels[k].sign);
    EXPECT_DOUBLE_EQ(loaded.levels[k].energy, saved.levels[k].energy);
    for (int i = 0; i < saved.channels; i++)
      EXPECT_EQ(loaded.levels[k].coupling(i), saved.levels[k].coupling(i));
  }
  // Theta is recomputed on loading, in the same order, so it is identical rather than merely close.
  EXPECT_EQ((loaded.theta - saved.theta).cwiseAbs().maxCoeff(), 0.0);
}

} // namespace

TEST(MixChainStarIO, round_trip_preserves_a_real_star) { // NOLINT
  const auto filename = "star_io_real.dat";
  const auto saved    = build_star(make_input<double>(real_gamma), options_for());
  save_star(saved, filename);

  const auto header = read_star_header(filename);
  EXPECT_EQ(header.channels, 2);
  EXPECT_EQ(header.mMAX, mmax);
  EXPECT_DOUBLE_EQ(header.bandrescale, 2.5);
  EXPECT_FALSE(header.complex_data);
  EXPECT_FALSE(star_is_complex(filename));

  expect_same_star(saved, load_star<double>(filename));
  std::remove(filename);
}

TEST(MixChainStarIO, round_trip_preserves_a_complex_star) { // NOLINT
  const auto filename = "star_io_complex.dat";
  const auto saved    = build_star(make_input<Complex>(complex_gamma), options_for());
  save_star(saved, filename);

  EXPECT_TRUE(star_is_complex(filename));
  expect_same_star(saved, load_star<Complex>(filename));
  std::remove(filename);
}

TEST(MixChainStarIO, writes_one_row_per_level_with_the_expected_columns) { // NOLINT
  const auto saved = build_star(make_input<double>(real_gamma), options_for());
  std::ostringstream out;
  save_star(saved, out);

  std::istringstream in(out.str());
  std::string line;
  std::size_t rows = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line.front() == '#') continue;
    rows++;
    std::istringstream fields(line);
    std::string token;
    std::size_t columns = 0;
    while (fields >> token) columns++;
    EXPECT_EQ(columns, 4U + 2U); // m, sign, a, E and two real components
  }
  EXPECT_EQ(rows, 2U * 2U * (mmax + 1));
}

TEST(MixChainStarIO, rejects_loading_with_the_wrong_scalar_type) { // NOLINT
  const auto filename = "star_io_scalar.dat";
  save_star(build_star(make_input<Complex>(complex_gamma), options_for()), filename);
  EXPECT_THROW(load_star<double>(filename), std::logic_error);
  std::remove(filename);
}

TEST(MixChainStarIO, reads_a_minimal_star) { // NOLINT
  const auto filename = "star_io_minimal.dat";
  write_file(filename, minimal_star);
  const auto star = load_star<double>(filename);
  EXPECT_EQ(star.levels.size(), 4U);
  EXPECT_DOUBLE_EQ(star.levels[0].energy, 0.7);
  EXPECT_DOUBLE_EQ(star.levels[2].energy, -0.7);
  EXPECT_DOUBLE_EQ(star.theta(0, 0), 2.0 * (0.25 + 0.09));
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_a_file_without_a_header) { // NOLINT
  const auto filename = "star_io_noheader.dat";
  write_file(filename, "# mixchain star\n0 + 0 0.7 0.5\n");
  EXPECT_THROW(read_star_header(filename), std::runtime_error);
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_an_incomplete_header) { // NOLINT
  const auto filename = "star_io_partial_header.dat";
  write_file(filename, "# channels=1 mMAX=1 z=1\n0 + 0 0.7 0.5\n");
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_a_truncated_file) { // NOLINT
  const auto filename = "star_io_truncated.dat";
  write_file(filename, "# mixchain star\n"
                       "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0\n"
                       "0 + 0 0.7 0.5\n"
                       "1 + 0 0.35 0.3\n");
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_a_row_with_the_wrong_number_of_columns) { // NOLINT
  const auto filename = "star_io_columns.dat";
  write_file(filename, "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0\n"
                       "0 + 0 0.7 0.5 0.1\n");
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_inconsistent_rows) { // NOLINT
  const auto base = std::string("# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0\n");
  const auto filename = "star_io_rows.dat";

  write_file(filename, base + "0 * 0 0.7 0.5\n"); // the sign column
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  write_file(filename, base + "0 + 0 -0.7 0.5\n"); // the energy contradicts the frequency branch
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  write_file(filename, base + "0 + 0 0.0 0.5\n"); // a vanishing energy
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  write_file(filename, base + "7 + 0 0.7 0.5\n"); // the interval index
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  write_file(filename, base + "0 + 3 0.7 0.5\n"); // the branch index
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  std::remove(filename);
}

TEST(MixChainStarIO, rejects_an_invalid_header) { // NOLINT
  const auto filename = "star_io_bad_header.dat";
  for (const auto *header : {"# channels=0 mMAX=1 z=1 Lambda=2 complex=0\n",
                             "# channels=1 mMAX=0 z=1 Lambda=2 complex=0\n",
                             "# channels=1 mMAX=1 z=0 Lambda=2 complex=0\n",
                             "# channels=1 mMAX=1 z=1 Lambda=1 complex=0\n",
                             "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=0 complex=0\n"}) {
    write_file(filename, header);
    EXPECT_THROW(read_star_header(filename), std::runtime_error);
  }
  std::remove(filename);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
