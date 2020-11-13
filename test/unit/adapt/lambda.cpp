#include <gtest/gtest.h>

#include <cmath>
#include <adapt/lambda.hpp>

using namespace NRG::Adapt;

TEST(Lambda, basic) { // NOLINT
  LAMBDA Lambda(2.0);
  EXPECT_EQ(Lambda, 2.0); // NOLINT
  EXPECT_EQ(Lambda.logL(), log(2.0));
  EXPECT_EQ(Lambda.factor(), (1.0-1.0/2.0)/log(2.0));
  EXPECT_EQ(Lambda.power(2.0), pow(2.0, 2.0));
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
