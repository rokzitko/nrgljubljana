#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>

using namespace std;

#include <nrgchain/io.h>
#include <nrgchain/parser.h>

TEST(NRGChainParser, skips_blank_lines) { // NOLINT
  const auto filename = "nrgchain_blank_lines.param";
  {
    ofstream file(filename);
    file << "[param]\n\nLambda=4\n\n# comment\n\nNmax=8\n";
  }

  params.clear();
  parser(filename);

  EXPECT_EQ(P("Lambda", 0.0), 4.0);
  EXPECT_EQ(Pint("Nmax", 0), 8);

  std::remove(filename);
}

TEST(NRGChainParser, trims_keys_and_values) { // NOLINT
  const auto filename = "nrgchain_whitespace.param";
  {
    ofstream file(filename);
    file << "[param]\nLambda = 2.0\nadapt = true\nband = custom.dat\nNmax\t=\t8\n";
  }

  params.clear();
  parser(filename);

  EXPECT_EQ(P("Lambda", 0.0), 2.0);
  EXPECT_TRUE(Pbool("adapt", false));
  EXPECT_EQ(Pstr("band", ""), "custom.dat");
  EXPECT_EQ(Pint("Nmax", 0), 8);

  std::remove(filename);
}

TEST(NRGChainParser, accepts_extended_bool_values_and_rejects_invalid_values) { // NOLINT
  const auto filename = "nrgchain_bool.param";
  {
    ofstream file(filename);
    file << "[param]\n"
         << "yes = YeS\n"
         << "true = tRuE\n"
         << "one = 1\n"
         << "no = nO\n"
         << "false = FaLsE\n"
         << "zero = 0\n"
         << "invalid = maybe\n";
  }

  params.clear();
  parser(filename);

  EXPECT_TRUE(Pbool("yes", false));
  EXPECT_TRUE(Pbool("true", false));
  EXPECT_TRUE(Pbool("one", false));
  EXPECT_FALSE(Pbool("no", true));
  EXPECT_FALSE(Pbool("false", true));
  EXPECT_FALSE(Pbool("zero", true));
  EXPECT_TRUE(Pbool("missing", true));
  EXPECT_THROW(Pbool("invalid", false), std::runtime_error);

  std::remove(filename);
}

TEST(NRGChainParser, rejects_partial_nonfinite_and_out_of_range_numbers) { // NOLINT
  const auto filename = "nrgchain_invalid_numbers.param";
  {
    ofstream file(filename);
    file << "[param]\n"
         << "partial = 2junk\n"
         << "infinite = inf\n"
         << "integer = 999999999999999999999\n";
  }

  params.clear();
  parser(filename);
  EXPECT_THROW(P("partial", 0.0), std::invalid_argument);
  EXPECT_THROW(P("infinite", 0.0), std::invalid_argument);
  EXPECT_THROW(Pint("integer", 0), std::invalid_argument);

  std::remove(filename);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
