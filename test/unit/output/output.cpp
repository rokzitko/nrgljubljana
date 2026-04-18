#include <gtest/gtest.h>

#include <cstdio>
#include <map>
#include <string>

#include <output.hpp>

using namespace std::string_literals;
using namespace NRG;

TEST(Output, AnnotatedHandlesEmptyDiag) { // NOLINT
  Params P;
  P.dumpannotated = 1;
  Stats<double> stats(P, {}, 0.0);
  Step step(P);
  DiagInfo<double> diag;
  Annotated annotated(P);
  const auto filename = "annotated-empty.dat"s;

  EXPECT_NO_THROW(annotated.dump(step, diag, stats, []([[maybe_unused]] const auto &I) { return 1; }, filename));

  std::remove(filename.c_str());
}

TEST(Output, ExpvOutputThrowsOnBadPath) { // NOLINT
  Params P;
  std::map<std::string, double> values;
  EXPECT_THROW(ExpvOutput<double>("/no/such/dir/custom.dat", values, {}, P), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
