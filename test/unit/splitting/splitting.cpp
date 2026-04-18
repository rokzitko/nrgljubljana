#include <gtest/gtest.h>

#include <splitting.hpp>

using namespace NRG;

TEST(Splitting, FixesTrailingCluster) { // NOLINT
  Params P;
  DiagInfo<double> diag;
  diag[Invar()] = NRG::Eigen<double>(std::vector<double>{0.0, 1.0, 1.0 + 1e-8}, 1.0);

  Clusters<double> clusters(diag, 1e-6, P);
  const auto corrected = diag.begin()->second.values.all_corr();

  EXPECT_DOUBLE_EQ(corrected[1], corrected[2]);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
