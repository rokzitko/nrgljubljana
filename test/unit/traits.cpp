#include <gtest/gtest.h>
#include <type_traits>
#include <traits.hpp>

using namespace NRG;

TEST(traits, double) {
  {
    bool t = std::is_same_v<typename traits<double>::t_matel, double>;
    EXPECT_EQ(t, true);
  }
}

TEST(traits, complex) {
  {
    bool t = std::is_same_v<typename traits<std::complex<double>>::t_matel, std::complex<double>>;
    EXPECT_EQ(t, true);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
