#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
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

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
