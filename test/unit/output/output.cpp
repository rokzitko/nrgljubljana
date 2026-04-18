#include <gtest/gtest.h>

#include <cstdio>
#include <map>
#include <string>

#include <h5.hpp>
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

TEST(Output, StatsH5SaveStoresAbsoluteGroundStateEnergy) { // NOLINT
  Params P;
  P.Nmax = 2;

  Stats<double> stats(P, {}, 0.0);
  Step step(P);

  stats.Egs = 2.0;
  stats.update(step);
  step.next();
  stats.Egs = 3.0;
  stats.update(step);

  H5Easy::File file("stats.h5", H5Easy::File::Overwrite);
  stats.h5save_nrg(file);

  EXPECT_DOUBLE_EQ(H5Easy::load<std::vector<double>>(file, "stats/1/rel_Egs").front(), stats.rel_Egs[0]);
  EXPECT_DOUBLE_EQ(H5Easy::load<std::vector<double>>(file, "stats/1/abs_Egs").front(), stats.abs_Egs[0]);
  EXPECT_DOUBLE_EQ(H5Easy::load<std::vector<double>>(file, "stats/2/abs_Egs").front(), stats.abs_Egs[1]);

  std::remove("stats.h5");
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
