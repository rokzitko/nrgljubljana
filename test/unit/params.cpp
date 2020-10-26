#include <gtest/gtest.h>
#include <params.hpp>

using namespace NRG;

TEST(parser, params) {
  {
    std::list<parambase *> all;
    param<double> p{"p", "Testing paraemter", 1.0, all};
    EXPECT_EQ(p.value(), 1.0);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
