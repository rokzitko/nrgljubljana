#include <gtest/gtest.h>
#include <nrg-general.hpp> // common
#include <nrg-lib.hpp>     // exposed in library

TEST(int, Add) { // NOLINT
   int a{0};
   int b{2};
   auto c = a + b;
   EXPECT_EQ(c, b); // NOLINT
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
