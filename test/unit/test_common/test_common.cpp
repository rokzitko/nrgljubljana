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
  auto Sym = set_symmetry<double>(P, "QS", 1);
  EXPECT_EQ(P.combs, 4);
  EXPECT_EQ(Sym->nr_combs(), 4);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
