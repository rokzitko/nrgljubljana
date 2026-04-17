#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>

#include <adapt/adapt.hpp>

using namespace NRG::Adapt;

namespace {

void write_file(const std::string &filename, const std::string &contents) {
  std::ofstream file(filename);
  file << contents;
}

} // namespace

TEST(Adapt, parser_skips_blank_lines) { // NOLINT
  const auto filename = "adapt_blank_lines.param";
  write_file(filename, "[param]\n\nLambda=3\n\n# comment\n\nxmax=7\n");

  Params P(filename);
  EXPECT_EQ(P.P("Lambda", 0.0), 3.0);
  EXPECT_EQ(P.P("xmax", 0.0), 7.0);

  std::remove(filename);
}

TEST(Adapt, linint_requires_two_points) { // NOLINT
  EXPECT_THROW(LinInt(Vec{}), std::runtime_error);
  EXPECT_THROW(LinInt(Vec{{1.0, 2.0}}), std::runtime_error);
}

TEST(Adapt, int_with_to_throws_when_step_limit_is_exceeded) { // NOLINT
  const auto filename = "adapt_int_with_to.param";
  write_file(filename, "[param]\nLambda=2\nadapt=false\nxmax=2\nxfine=2\n");

  Params P(filename);
  Adapt calc(P, Sign::POS);
  calc.x = 0.0;
  calc.y = 0.0;

  EXPECT_THROW(calc.int_with_to(0.0, 1.0, []([[maybe_unused]] const auto x, [[maybe_unused]] const auto y) { return 0.0; }, false, 1e-10, 0), std::runtime_error);

  std::remove(filename);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
