#include <gtest/gtest.h>

#include <recalc.hpp>
#include "test_common.hpp"

TEST(recalc, split_in_blocks_Eigen) { // NOLINT
   int a{0};
   int b{2};
   auto c = a + b;
   EXPECT_EQ(c, b); // NOLINT
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
