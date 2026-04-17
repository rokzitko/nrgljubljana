#include <gtest/gtest.h>

#include <kk/kk.hpp>

TEST(KK, empty_input_throws) { // NOLINT
  auto construct = [] { return NRG::KK::KK(NRG::KK::XYFUNC{}); };
  EXPECT_THROW(construct(), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
