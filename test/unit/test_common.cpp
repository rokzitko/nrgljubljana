#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <invar.hpp>
#include <subspaces.hpp>
#include <core.hpp>

using namespace NRG;

TEST(test_common, basic) {
  Params P;
  P.symtype.setvalue("QS");
  P.set_channels(1);
  EXPECT_EQ(P.combs, 4);
  auto Sym = set_symmetry<double>(P);
  EXPECT_EQ(Sym->nr_combs(), 4);
}

TEST(test_common, basic2) {
//  auto [P, Sym] = test_setup_basic<double>();
//  EXPECT_EQ(P.combs, 4);
//  EXPECT_EQ(Sym->nr_combs(), 4);
}

TEST(test_common, diag) {
//  auto [P, Sym] = test_setup_basic<double>();
//  auto diag = test_setup_diag(P, Sym.get());
// EXPECT_EQ(Sym->nr_combs(), 4);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
