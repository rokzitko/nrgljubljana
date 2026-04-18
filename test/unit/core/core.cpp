#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <core.hpp>
#include <truncation.hpp>

using namespace NRG;

TEST(Core, HighestRetainedUsesSafeguardWithoutFullSort) { // NOLINT
  Params P;
  P.keep = 2;
  P.keepmin = 1;
  P.safeguard = 0.1;
  P.safeguardmax = 2;

  DiagInfo<double> diag;
  diag[Invar()] = NRG::Eigen<double>(std::vector<double>{2.0, 0.0, 1.0, 1.05}, 1.0);

  Step step(P);
  EXPECT_DOUBLE_EQ(highest_retained(step, diag, P), 1.05);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
