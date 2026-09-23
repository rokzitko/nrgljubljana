#include <gtest/gtest.h>

#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <mixchain/load.hpp>

using namespace NRG::MixChain;

namespace {

using Table = std::vector<std::pair<double, double>>;

// The grid used by most of the tests: two nodes on each side, reaching the band edge.
const std::vector<double> standard_grid{-1.0, -0.5, 0.5, 1.0};

Table table(const std::vector<double> &grid, const std::vector<double> &values) {
  Table result;
  for (std::size_t k = 0; k < grid.size(); k++) result.emplace_back(grid[k], values[k]);
  return result;
}

Table constant_table(const std::vector<double> &grid, const double value) {
  return table(grid, std::vector<double>(grid.size(), value));
}

void write_component(const std::string &prefix, const int i, const int j, const bool imaginary, const Table &values) {
  std::ofstream file(gamma_filename(prefix, i, j, imaginary));
  file << std::setprecision(17);
  for (const auto &[omega, value] : values) file << omega << " " << value << "\n";
}

// Each test uses its own prefix, so no test can be affected by the files of another.
GammaOptions options_for(const std::string &prefix, const int channels) {
  GammaOptions options;
  options.prefix   = prefix;
  options.channels = channels;
  return options;
}

void remove_components(const std::string &prefix, const int channels) {
  for (int i = 1; i <= channels; i++) {
    for (int j = 1; j <= channels; j++) {
      std::remove(gamma_filename(prefix, i, j, false).c_str());
      std::remove(gamma_filename(prefix, i, j, true).c_str());
    }
  }
}

// A 2x2 real symmetric input with distinct entries.
void write_real_2x2(const std::string &prefix) {
  write_component(prefix, 1, 1, false, constant_table(standard_grid, 0.3));
  write_component(prefix, 1, 2, false, constant_table(standard_grid, 0.1));
  write_component(prefix, 2, 1, false, constant_table(standard_grid, 0.1));
  write_component(prefix, 2, 2, false, constant_table(standard_grid, 0.4));
}

} // namespace

TEST(MixChainGammaInput, reads_a_single_channel_real_input) { // NOLINT
  const std::string prefix = "gi_scalar";
  // Distinct values everywhere, so that a mirrored or misordered branch cannot pass.
  write_component(prefix, 1, 1, false, table(standard_grid, {0.1, 0.2, 0.3, 0.4}));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  EXPECT_FALSE(gamma_is_complex(options));
  const auto input = load_gamma<double>(options, out);

  ASSERT_EQ(input.pos.size(), 3U); // two nodes plus the zero point
  ASSERT_EQ(input.neg.size(), 3U);
  EXPECT_DOUBLE_EQ(input.pos.omega[0], zero_point);
  EXPECT_DOUBLE_EQ(input.pos.omega[1], 0.5);
  EXPECT_DOUBLE_EQ(input.pos.omega[2], 1.0);
  EXPECT_DOUBLE_EQ(input.pos.gamma[0](0, 0), 0.3); // the zero point continues the innermost node
  EXPECT_DOUBLE_EQ(input.pos.gamma[1](0, 0), 0.3);
  EXPECT_DOUBLE_EQ(input.pos.gamma[2](0, 0), 0.4);

  // The negative branch is indexed by |omega|, so it is the input read inwards.
  EXPECT_DOUBLE_EQ(input.neg.omega[1], 0.5);
  EXPECT_DOUBLE_EQ(input.neg.omega[2], 1.0);
  EXPECT_DOUBLE_EQ(input.neg.gamma[0](0, 0), 0.2);
  EXPECT_DOUBLE_EQ(input.neg.gamma[1](0, 0), 0.2);
  EXPECT_DOUBLE_EQ(input.neg.gamma[2](0, 0), 0.1);

  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, only_the_offdiagonal_imaginary_parts_decide_the_arithmetic) { // NOLINT
  const std::string prefix = "gi_complexity";
  write_real_2x2(prefix);
  auto options = options_for(prefix, 2);
  EXPECT_FALSE(gamma_is_complex(options));

  // A diagonal imaginary part must vanish, so it does not make the problem complex.
  write_component(prefix, 1, 1, true, constant_table(standard_grid, 0.0));
  EXPECT_FALSE(gamma_is_complex(options));

  write_component(prefix, 1, 2, true, constant_table(standard_grid, 0.05));
  write_component(prefix, 2, 1, true, constant_table(standard_grid, -0.05));
  EXPECT_TRUE(gamma_is_complex(options));

  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, reads_a_complex_hermitian_input) { // NOLINT
  const std::string prefix = "gi_complex";
  write_real_2x2(prefix);
  write_component(prefix, 1, 2, true, constant_table(standard_grid, 0.05));
  write_component(prefix, 2, 1, true, constant_table(standard_grid, -0.05));
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  ASSERT_TRUE(gamma_is_complex(options));
  const auto input = load_gamma<std::complex<double>>(options, out);

  const auto &m = input.pos.gamma[1];
  EXPECT_DOUBLE_EQ(m(0, 0).real(), 0.3);
  EXPECT_DOUBLE_EQ(m(0, 0).imag(), 0.0);
  EXPECT_DOUBLE_EQ(m(0, 1).real(), 0.1);
  EXPECT_DOUBLE_EQ(m(0, 1).imag(), 0.05);
  EXPECT_DOUBLE_EQ(m(1, 0).real(), 0.1);
  EXPECT_DOUBLE_EQ(m(1, 0).imag(), -0.05);
  EXPECT_DOUBLE_EQ(m(1, 1).real(), 0.4);

  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_a_partial_set_of_offdiagonal_imaginary_parts) { // NOLINT
  const std::string prefix = "gi_partial_im";
  write_real_2x2(prefix);
  write_component(prefix, 1, 2, true, constant_table(standard_grid, 0.05));
  auto options = options_for(prefix, 2);

  EXPECT_THROW(gamma_is_complex(options), std::runtime_error);
  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_a_non_hermitian_complex_input) { // NOLINT
  const std::string prefix = "gi_non_hermitian";
  write_real_2x2(prefix);
  write_component(prefix, 1, 2, true, constant_table(standard_grid, 0.05));
  write_component(prefix, 2, 1, true, constant_table(standard_grid, 0.05)); // should be -0.05
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  EXPECT_THROW(load_gamma<std::complex<double>>(options, out), std::runtime_error);
  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_a_non_symmetric_real_input) { // NOLINT
  const std::string prefix = "gi_non_symmetric";
  write_real_2x2(prefix);
  write_component(prefix, 2, 1, false, constant_table(standard_grid, 0.2)); // should be 0.1
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);
  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, symmetrizes_a_deviation_within_the_tolerance) { // NOLINT
  const std::string prefix = "gi_symmetrize";
  write_real_2x2(prefix);
  write_component(prefix, 2, 1, false, constant_table(standard_grid, 0.1 + 1e-10));
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  const auto input = load_gamma<double>(options, out);
  const auto &m    = input.pos.gamma[1];
  EXPECT_DOUBLE_EQ(m(0, 1), m(1, 0));
  EXPECT_NEAR(m(0, 1), 0.1 + 0.5e-10, 1e-16);

  // The same deviation scaled up past the tolerance is an error.
  options.hermiticity_tolerance = 1e-12;
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);

  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_a_nonzero_diagonal_imaginary_part_in_a_real_run) { // NOLINT
  const std::string prefix = "gi_diagonal_im";
  write_component(prefix, 1, 1, false, constant_table(standard_grid, 0.3));
  write_component(prefix, 1, 1, true, constant_table(standard_grid, 0.01));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  ASSERT_FALSE(gamma_is_complex(options));
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);

  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, rejects_complex_input_in_real_arithmetic) { // NOLINT
  const std::string prefix = "gi_wrong_scalar";
  write_real_2x2(prefix);
  write_component(prefix, 1, 2, true, constant_table(standard_grid, 0.05));
  write_component(prefix, 2, 1, true, constant_table(standard_grid, -0.05));
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  EXPECT_THROW(load_gamma<double>(options, out), std::logic_error);
  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_grids_that_differ) { // NOLINT
  const std::string prefix = "gi_grid";
  write_real_2x2(prefix);
  auto options = options_for(prefix, 2);
  std::ostringstream out;

  // A different number of nodes.
  write_component(prefix, 2, 2, false, constant_table({-1.0, -0.5, 0.0, 0.5, 1.0}, 0.4));
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);

  // The same number of nodes, one of them displaced.
  write_component(prefix, 2, 2, false, constant_table({-1.0, -0.5, 0.6, 1.0}, 0.4));
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);

  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_a_non_monotonic_grid) { // NOLINT
  const std::string prefix = "gi_monotonic";
  write_component(prefix, 1, 1, false, constant_table({-1.0, 0.5, -0.5, 1.0}, 0.3));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);
  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, rejects_a_missing_component_file) { // NOLINT
  const std::string prefix = "gi_missing";
  write_component(prefix, 1, 1, false, constant_table(standard_grid, 0.3));
  write_component(prefix, 1, 2, false, constant_table(standard_grid, 0.1));
  write_component(prefix, 2, 1, false, constant_table(standard_grid, 0.1));
  auto options = options_for(prefix, 2); // Gamma_22-re.dat was never written

  std::ostringstream out;
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);
  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, rejects_one_sided_data) { // NOLINT
  const std::string prefix = "gi_one_sided";
  write_component(prefix, 1, 1, false, constant_table({0.25, 0.5, 1.0}, 0.3));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  EXPECT_THROW(load_gamma<double>(options, out), std::runtime_error);
  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, rejects_invalid_options) { // NOLINT
  const std::string prefix = "gi_options";
  write_component(prefix, 1, 1, false, constant_table(standard_grid, 0.3));

  auto options     = options_for(prefix, 1);
  options.channels = 0;
  EXPECT_THROW(gamma_is_complex(options), std::invalid_argument);
  options.channels = max_channels + 1;
  EXPECT_THROW(gamma_is_complex(options), std::invalid_argument);

  options             = options_for(prefix, 1);
  options.bandrescale = 0.0;
  EXPECT_THROW(gamma_is_complex(options), std::invalid_argument);

  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, applies_bandrescale) { // NOLINT
  const std::string prefix = "gi_bandrescale";
  write_component(prefix, 1, 1, false, constant_table({-2.0, -1.0, 1.0, 2.0}, 0.25));
  auto options        = options_for(prefix, 1);
  options.bandrescale = 2.0;

  std::ostringstream out;
  const auto input = load_gamma<double>(options, out);

  // The band edge is mapped to 1 and the weight is preserved: Gamma is scaled by bandrescale.
  EXPECT_DOUBLE_EQ(input.pos.omega.back(), 1.0);
  EXPECT_DOUBLE_EQ(input.pos.omega[1], 0.5);
  EXPECT_DOUBLE_EQ(input.pos.gamma.back()(0, 0), 0.5);
  EXPECT_EQ(out.str().find("discarded"), std::string::npos);

  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, reports_weight_discarded_beyond_the_band_edge) { // NOLINT
  const std::string prefix = "gi_discarded";
  const std::vector<double> wide_grid{-2.0, -1.0, 1.0, 2.0};
  write_component(prefix, 1, 1, false, constant_table(wide_grid, 0.25));
  write_component(prefix, 1, 2, false, constant_table(wide_grid, 0.0));
  write_component(prefix, 2, 1, false, constant_table(wide_grid, 0.0));
  write_component(prefix, 2, 2, false, constant_table(wide_grid, 0.0));
  auto options = options_for(prefix, 2);

  std::ostringstream out;
  load_gamma<double>(options, out);
  const auto report = out.str();

  // Half of the weight of Gamma_11 is beyond |omega|=1. Gamma_22 carries no weight at all, so it has none to lose.
  EXPECT_NE(report.find("# Gamma_11: 0.5 of 1"), std::string::npos);
  EXPECT_NE(report.find("discarded"), std::string::npos);
  EXPECT_EQ(report.find("# Gamma_22: "), std::string::npos);

  remove_components(prefix, 2);
}

TEST(MixChainGammaInput, reports_weight_added_by_extrapolation) { // NOLINT
  const std::string prefix = "gi_extrapolated";
  write_component(prefix, 1, 1, false, constant_table({-0.5, -0.25, 0.25, 0.5}, 0.4));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  load_gamma<double>(options, out);
  const auto report = out.str();

  // The density is continued at 0.4 over (0.5,1] on each side: 0.4 of weight added to 0.4 tabulated.
  EXPECT_NE(report.find("# Gamma_11: 0.4 of 0.4"), std::string::npos);
  EXPECT_NE(report.find("extrapolation"), std::string::npos);

  remove_components(prefix, 1);
}

TEST(MixChainGammaInput, reports_nothing_when_the_band_is_covered) { // NOLINT
  const std::string prefix = "gi_covered";
  write_component(prefix, 1, 1, false, constant_table(standard_grid, 0.3));
  auto options = options_for(prefix, 1);

  std::ostringstream out;
  load_gamma<double>(options, out);
  const auto report = out.str();

  EXPECT_EQ(report.find("discarded"), std::string::npos);
  EXPECT_EQ(report.find("extrapolation"), std::string::npos);

  remove_components(prefix, 1);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
