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

TEST(MixChainStarIO, round_trip_preserves_the_blocks) { // NOLINT
  const auto filename = "star_io_blocks.dat";
  // A diagonal Gamma splits into one block per channel.
  const auto diagonal = [](const double omega) {
    Matrix<double> m = Matrix<double>::Zero(2, 2);
    m(0, 0)          = 0.8 - 0.1 * omega;
    m(1, 1)          = 0.3 + 0.1 * omega;
    return m;
  };
  const auto star = build_star(make_input<double>(diagonal), options_for());
  ASSERT_EQ(star.blocks, (Blocks{{0}, {1}}));
  std::ostringstream out;
  save_star(star, out);
  EXPECT_NE(out.str().find("\n# blocks= {1} {2}\n"), std::string::npos);
  save_star(star, filename);
  const auto loaded = load_star<double>(filename);
  expect_same_star(star, loaded);
  EXPECT_EQ(loaded.blocks, star.blocks);

  // Blocks that are not contiguous: channels 1 and 3 are coupled.
  const auto three = [](const double omega) {
    Matrix<Complex> m = Matrix<Complex>::Zero(3, 3);
    m(0, 0)           = 0.5 + 0.2 * omega;
    m(1, 1)           = 0.3 + omega * omega;
    m(2, 2)           = 0.4 - 0.1 * omega;
    m(0, 2)           = Complex(0.1 * omega, 0.05);
    m(2, 0)           = std::conj(m(0, 2));
    return m;
  };
  const auto split = build_star(make_input<Complex>(three), options_for());
  ASSERT_EQ(split.blocks, (Blocks{{0, 2}, {1}}));
  save_star(split, filename);
  const auto loaded_split = load_star<Complex>(filename);
  expect_same_star(split, loaded_split);
  EXPECT_EQ(loaded_split.blocks, split.blocks);
  std::remove(filename);
}

TEST(MixChainStarIO, the_untabulated_region_round_trips) { // NOLINT
  const auto filename = "star_io_untabulated.dat";
  auto input          = make_input<double>(real_gamma);
  input.pos.innermost = 0.2;  // as the branches record it, in the rescaled band
  input.neg.innermost = 0.15; // the wider of the two regions is what the star keeps
  const auto star     = build_star(input, options_for());
  ASSERT_TRUE(star.untabulated_known);
  EXPECT_EQ(star.untabulated_from, 0.0); // the plain mesh accumulates at 0
  EXPECT_DOUBLE_EQ(star.untabulated_to, 0.2);
  std::ostringstream out;
  save_star(star, out);
  EXPECT_NE(out.str().find("untabulated=0,0.2"), std::string::npos) << out.str().substr(0, 200);
  save_star(star, filename);
  const auto loaded = load_star<double>(filename);
  EXPECT_TRUE(loaded.untabulated_known);
  EXPECT_EQ(loaded.untabulated_from, 0.0);
  EXPECT_DOUBLE_EQ(loaded.untabulated_to, 0.2);

  // Tabulated down to omega = 0: no region.
  input.pos.innermost = 0.0;
  input.neg.innermost = 0.0;
  std::ostringstream none;
  save_star(build_star(input, options_for()), none);
  EXPECT_NE(none.str().find("untabulated=none"), std::string::npos);

  // A file that does not record it, as those written before it existed.
  write_file(filename, minimal_star);
  EXPECT_FALSE(load_star<double>(filename).untabulated_known);

  // A malformed value.
  write_file(filename, "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0 untabulated=0.2\n");
  EXPECT_THROW(read_star_header(filename), std::runtime_error);
  write_file(filename, "# channels=1 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0 untabulated=0.3,0.2\n");
  EXPECT_THROW(read_star_header(filename), std::runtime_error);
  std::remove(filename);
}

TEST(MixChainStarIO, a_single_block_writes_no_blocks_line) { // NOLINT
  auto options         = options_for();
  const auto whole     = build_star(make_input<double>(real_gamma), options);
  options.split_blocks = false;
  const auto unsplit   = build_star(make_input<double>([](const double omega) {
                                    Matrix<double> m = Matrix<double>::Identity(2, 2);
                                    m(0, 0)          = 0.5 + omega;
                                    return m;
                                  }),
                                  options);
  const auto filename  = "star_io_single_block.dat";
  for (const auto *star : {&whole, &unsplit}) {
    ASSERT_EQ(star->blocks, (Blocks{{0, 1}}));
    std::ostringstream out;
    save_star(*star, out);
    EXPECT_EQ(out.str().find("blocks="), std::string::npos);
    save_star(*star, filename);
    EXPECT_EQ(load_star<double>(filename).blocks, (Blocks{{0, 1}}));
  }

  // A file written before blocks existed is a single block.
  write_file(filename, minimal_star);
  EXPECT_EQ(load_star<double>(filename).blocks, (Blocks{{0}}));
  std::remove(filename);
}

TEST(MixChainStarIO, rejects_inconsistent_blocks) { // NOLINT
  const auto filename = "star_io_bad_blocks.dat";
  // mMAX=1 and 2 channels: the file must hold 2*2*(1+1) = 8 rows.
  const auto header = std::string("# channels=2 mMAX=1 z=1 Lambda=2 bandrescale=1 complex=0\n");
  // Two blocks of one channel each: branch 0 belongs to channel 1, branch 1 to channel 2.
  const auto rows = [](const std::string &first) {
    return first + "0 + 1 0.7 0 0.4\n"
                   "1 + 0 0.35 0.3 0\n1 + 1 0.35 0 0.2\n"
                   "0 - 0 -0.7 0.5 0\n0 - 1 -0.7 0 0.4\n"
                   "1 - 0 -0.35 0.3 0\n1 - 1 -0.35 0 0.2\n";
  };

  write_file(filename, header + "# blocks= {1} {2}\n" + rows("0 + 0 0.7 0.5 0\n"));
  EXPECT_EQ(load_star<double>(filename).blocks, (Blocks{{0}, {1}}));
  write_file(filename, header + "# blocks= {1} {2}\n" + rows("0 + 0 0.7 0 0\n")); // a level with no coupling
  EXPECT_NO_THROW(load_star<double>(filename));

  write_file(filename, header + "# blocks= {1} {2}\n" + rows("0 + 0 0.7 0 0.5\n")); // in the other block
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);
  write_file(filename, header + "# blocks= {1} {2}\n" + rows("0 + 0 0.7 0.5 0.5\n")); // across both blocks
  EXPECT_THROW(load_star<double>(filename), std::runtime_error);

  for (const auto *line : {"# blocks= {1,1}\n", "# blocks= {1}\n", "# blocks= {1} {2}\n# blocks= {1} {2}\n"}) {
    write_file(filename, header + line + rows("0 + 0 0.7 0.5 0\n"));
    EXPECT_THROW(read_star_header(filename), std::runtime_error) << line;
  }
  std::remove(filename);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
