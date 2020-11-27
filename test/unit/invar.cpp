#include <gtest/gtest.h>

#include <invar.hpp>
#include <mk_sym.hpp>
#include "test_common.hpp"

TEST(invar, initInvar) { // NOLINT
  auto [P, sym] = test_setup_basic();
  // should have Q and SS degrees of freedom
  EXPECT_EQ(Invar::invdim, 2);
  EXPECT_EQ(Invar::qntype.size(), 2);
  EXPECT_EQ(Invar::qntype[0], additive);
  EXPECT_EQ(Invar::qntype[1], additive);
  EXPECT_EQ(Invar::names.size(), 2);
  EXPECT_EQ(Invar::names["Q"], 0);
  EXPECT_EQ(Invar::names["SS"], 1);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
