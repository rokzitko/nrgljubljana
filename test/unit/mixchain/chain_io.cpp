#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <mixchain/chain_io.hpp>
#include <mixchain/precision.hpp>

using namespace NRG::MixChain;

namespace {

using Real = WideReal<50>;

// A star with arbitrary energies and couplings; nothing here depends on it being a discretization of anything.
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

// The elements of one written chain, by block name, site and indices.
struct Written {
  std::vector<std::string> comments;
  std::map<std::tuple<std::string, int, int, int>, double> value;
};

Written parse(const std::string &text) {
  Written written;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind('#', 0) == 0) {
      written.comments.push_back(line);
      continue;
    }
    std::istringstream fields(line);
    std::string name;
    int n = 0, i = 0, j = 0;
    double value = 0;
    fields >> name >> n >> i >> j >> value;
    written.value[{name, n, i, j}] = value;
  }
  return written;
}

Written write(const Chain<Real> &chain, const double bandrescale) {
  std::ostringstream out;
  save_chain(chain, ChainFileHeader{1.0, 2.0, bandrescale, 50}, out);
  return parse(out.str());
}

} // namespace

TEST(MixChainChainIO, the_chain_is_written_in_the_units_of_the_input) { // NOLINT
  // E_n and T_n carry bandrescale back, as nrgchain applies it to xi.dat and zeta.dat; V is invariant under the
  // rescaling of omega and Gamma and is written as it is.
  const auto chain    = build_chain<Real>(arbitrary_star<double>(2, 12), [] {
    ChainOptions options;
    options.Nmax = 3;
    return options;
  }());
  const auto bare     = write(chain, 1.0);
  const auto rescaled = write(chain, 2.5);

  for (int i = 1; i <= 2; i++)
    for (int j = 1; j <= 2; j++) {
      const auto v = static_cast<double>(chain.V(i - 1, j - 1));
      EXPECT_EQ(bare.value.at({"V", 0, i, j}), v);
      EXPECT_EQ(rescaled.value.at({"V", 0, i, j}), v) << "V must not be rescaled";
      for (unsigned int n = 0; n <= chain.Nmax; n++) {
        const auto e = static_cast<double>(chain.E[n](i - 1, j - 1));
        EXPECT_EQ(bare.value.at({"E", static_cast<int>(n), i, j}), e);
        EXPECT_NEAR(rescaled.value.at({"E", static_cast<int>(n), i, j}), 2.5 * e, 1e-15 * std::abs(2.5 * e) + 1e-300);
      }
      for (unsigned int n = 0; n < chain.Nmax; n++) {
        const auto t = static_cast<double>(chain.T[n](i - 1, j - 1));
        EXPECT_EQ(bare.value.at({"T", static_cast<int>(n), i, j}), t);
        EXPECT_NEAR(rescaled.value.at({"T", static_cast<int>(n), i, j}), 2.5 * t, 1e-15 * std::abs(2.5 * t));
      }
    }
}

TEST(MixChainChainIO, the_header_records_the_run) { // NOLINT
  const auto chain = build_chain<Real>(arbitrary_star<double>(2, 12), [] {
    ChainOptions options;
    options.Nmax = 3;
    return options;
  }());
  const auto written = write(chain, 2.5);
  ASSERT_GE(written.comments.size(), 4U);
  EXPECT_EQ(written.comments[0], "# mixchain Wilson chain");
  EXPECT_NE(written.comments[1].find("channels=2 Nmax=3 z=1 Lambda=2 bandrescale=2.5 complex=0 digits=50"),
            std::string::npos)
    << written.comments[1];
  EXPECT_NE(written.comments[2].find("levels=12 coupled_levels=12 theta_rank=2"), std::string::npos)
    << written.comments[2];
  // A single block writes no blocks line.
  for (const auto &comment : written.comments) EXPECT_EQ(comment.find("blocks="), std::string::npos);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
