#include <gtest/gtest.h>
#include <numerics.hpp>

using namespace NRG;

TEST(reim, traits) {
  auto & [r, i] = reim(std::complex(1.0,2.0));
  EXPECT_EQ(r, 1.0);
  EXPECT_EQ(i, 2.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
