#include <gtest/gtest.h>
#include <nrgljubljana/nrgljubljana.hpp>

using namespace nrgljubljana;

TEST(Toto, Add) { // NOLINT

  toto a(0);
  toto b(2);

  auto c = a + b;
  EXPECT_EQ(c, b); // NOLINT
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS(); // NOLINT
}
