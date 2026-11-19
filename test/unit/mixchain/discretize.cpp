#include <gtest/gtest.h>

#include <complex>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mixchain/configuration.hpp>
#include <mixchain/discretize.hpp>
#include <mixchain/load.hpp>
#include <mixchain/star_io.hpp>

using namespace NRG::MixChain;

namespace {

// A 2x2 Gamma that mixes the channels, so that it does not split into blocks; complex adds an imaginary part to the
// off-diagonal elements.
GammaTable mixing_table(const bool complex) {
  GammaTable table;
  table.channels = 2;
  table.complex  = complex;
  for (int k = -50; k <= 50; k++) {
    const double w = k / 50.0;
    table.omega.push_back(w);
    const std::complex<double> offdiagonal{0.1 * (1.0 - w * w), complex ? 0.05 * w : 0.0};
    for (const auto value : {std::complex<double>{0.5, 0.0}, offdiagonal, std::conj(offdiagonal),
                             std::complex<double>{0.25 * (1.0 - w * w / 2.0), 0.0}})
      table.gamma.push_back(value);
  }
  return table;
}

void write_table_files(const GammaTable &table, const std::string &prefix) {
  for (int i = 0; i < table.channels; i++)
    for (int j = 0; j < table.channels; j++)
      for (const bool imaginary : {false, true}) {
        if (imaginary && (!table.complex || i == j)) continue;
        std::ofstream file(gamma_filename(prefix, i + 1, j + 1, imaginary));
        file.precision(17);
        for (std::size_t k = 0; k < table.omega.size(); k++)
          file << table.omega[k] << ' ' << (imaginary ? table.at(k, i, j).imag() : table.at(k, i, j).real()) << '\n';
      }
}

void remove_table_files(const GammaTable &table, const std::string &prefix) {
  for (int i = 1; i <= table.channels; i++)
    for (int j = 1; j <= table.channels; j++)
      for (const bool imaginary : {false, true}) std::remove(gamma_filename(prefix, i, j, imaginary).c_str());
}

template<typename S> void expect_same_branch(const GammaBranch<S> &a, const GammaBranch<S> &b) {
  EXPECT_EQ(a.omega, b.omega);
  EXPECT_EQ(a.innermost, b.innermost);
  ASSERT_EQ(a.gamma.size(), b.gamma.size());
  for (std::size_t k = 0; k < a.gamma.size(); k++) EXPECT_TRUE(a.gamma[k] == b.gamma[k]) << k;
}

template<typename S> void expect_same_chain(const Chain<S> &a, const Chain<S> &b) {
  EXPECT_EQ(a.Nmax, b.Nmax);
  EXPECT_EQ(a.blocks, b.blocks);
  EXPECT_TRUE(a.V == b.V);
  ASSERT_EQ(a.E.size(), b.E.size());
  ASSERT_EQ(a.T.size(), b.T.size());
  for (std::size_t n = 0; n < a.E.size(); n++) {
    EXPECT_TRUE(a.E[n] == b.E[n]) << n;
    EXPECT_TRUE(a.T[n] == b.T[n]) << n;
  }
  EXPECT_EQ(a.diagnostics.theta_rank, b.diagnostics.theta_rank);
  EXPECT_EQ(a.diagnostics.max_star_sensitivity, b.diagnostics.max_star_sensitivity);
  EXPECT_EQ(a.diagnostics.sensitive_from_site, b.diagnostics.sensitive_from_site);
}

const std::map<std::string, std::string> entries{{"channels", "2"}, {"Lambda", "2"}, {"Nmax", "12"}};

// The chain of each z in memory equals the one built, as the executable does, from the star reloaded from star.dat.
template<typename S> void in_memory_equals_file_round_trip(const bool complex) {
  const auto configuration = read_configuration(entries, true, true);
  std::ostringstream log;
  const auto input  = gamma_input<S>(mixing_table(complex), configuration.gamma, log);
  const auto result = discretize(input, configuration, {0.5, 1.0}, log);
  ASSERT_EQ(result.size(), 2u);

  StarDiscretizer<S> discretizer(input, configuration.star);
  for (const auto &d : result) {
    const auto filename = "discretize_star.dat";
    save_star(discretizer.star(d.star.z), filename);
    ASSERT_EQ(star_is_complex(filename), complex);
    const auto star = load_star<S>(filename);
    std::remove(filename);
    expect_same_chain(make_chain(star, configuration.chain), d.chain);
    EXPECT_EQ(d.chain.Nmax, 12u);
  }
}

} // namespace

TEST(MixChainDiscretize, table_equals_files) { // NOLINT
  auto options   = read_configuration(entries, true, false).gamma;
  options.prefix = "discretize_gamma";
  for (const bool complex : {false, true}) {
    const auto table = mixing_table(complex);
    write_table_files(table, options.prefix);
    std::ostringstream from_files, from_table;
    if (complex) {
      const auto a = load_gamma<std::complex<double>>(options, from_files);
      const auto b = gamma_input<std::complex<double>>(table, options, from_table);
      expect_same_branch(a.pos, b.pos);
      expect_same_branch(a.neg, b.neg);
    } else {
      const auto a = load_gamma<double>(options, from_files);
      const auto b = gamma_input<double>(table, options, from_table);
      expect_same_branch(a.pos, b.pos);
      expect_same_branch(a.neg, b.neg);
    }
    EXPECT_EQ(from_files.str(), from_table.str());
    remove_table_files(table, options.prefix);
  }
}

TEST(MixChainDiscretize, table_errors) { // NOLINT
  const GammaOptions options{.prefix = "unused", .channels = 2};
  std::ostringstream log;
  auto table = mixing_table(false);
  table.gamma.pop_back();
  EXPECT_THROW(gamma_input<double>(table, options, log), std::invalid_argument); // size

  table = mixing_table(false);
  table.omega[3] = table.omega[2];
  EXPECT_THROW(gamma_input<double>(table, options, log), std::runtime_error); // not increasing

  EXPECT_THROW(gamma_input<double>(mixing_table(true), options, log), std::logic_error); // complex in real arithmetic

  table = mixing_table(false);
  table.gamma[0] = {0.5, 0.1}; // imaginary diagonal element in real arithmetic
  EXPECT_THROW(gamma_input<double>(table, options, log), std::runtime_error);

  table = mixing_table(false);
  table.gamma[1] = 0.3; // Gamma_12 != Gamma_21
  EXPECT_THROW(gamma_input<double>(table, options, log), std::runtime_error);

  for (auto &value : table.gamma) value = 0.0;
  EXPECT_THROW(gamma_input<double>(table, options, log), std::runtime_error); // vanishes
}

TEST(MixChainDiscretize, configuration_from_entries) { // NOLINT
  const auto c = read_configuration(entries, true, true);
  EXPECT_EQ(c.gamma.channels, 2);
  EXPECT_EQ(c.star.mMAX, 24u);
  EXPECT_TRUE(c.mmax_from_nmax);
  EXPECT_EQ(c.chain.Nmax, 12u);
  EXPECT_EQ(c.chain.gauge, ChainGauge::polar);
  EXPECT_DOUBLE_EQ(c.star.z, 1.0);
  EXPECT_FALSE(c.discretization_files);

  auto e = entries;
  e["mMAX"] = "30";
  EXPECT_EQ(read_configuration(e, true, false).star.mMAX, 30u);
  EXPECT_THROW(read_configuration({{"Lambda", "2"}}, true, false), std::invalid_argument);  // neither mMAX nor Nmax
  EXPECT_THROW(read_configuration({{"mMAX", "20"}}, false, true), std::invalid_argument);   // the chain needs Nmax
  EXPECT_THROW(read_configuration({{"Nmax", "0"}}, true, true), std::invalid_argument);
  e["rank_tolerance"] = "-1";
  EXPECT_THROW(read_configuration(e, true, true), std::invalid_argument);
}

TEST(MixChainDiscretize, in_memory_equals_file_round_trip_real) { in_memory_equals_file_round_trip<double>(false); } // NOLINT

TEST(MixChainDiscretize, in_memory_equals_file_round_trip_complex) { // NOLINT
  in_memory_equals_file_round_trip<std::complex<double>>(true);
}
