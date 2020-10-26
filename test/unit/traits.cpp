#include <gtest/gtest.h>
#include <traits.hpp>

using namespace NRG;

TEST(double, traits) {
  {
    EXPECT_EQ(std::is_same_v(typename traits<double>::t_matel, double), true);
  }
}

TEST(complex, traits) {
  {
    EXPECT_EQ(std::is_same_v(typename traits<std::complex<double>>::t_matel, std::complex<double>), true);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
