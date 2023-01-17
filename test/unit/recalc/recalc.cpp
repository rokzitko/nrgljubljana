#include <gtest/gtest.h>

#include <recalc.hpp>
#include "test_common.hpp"

TEST(recalc, split_in_blocks_Eigen) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag(P, Sym);
#ifdef DISABLED_DUE_TO_ISSUES
  SubspaceStructure substruct{diag, Sym};
  substruct.dump();
#endif
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
