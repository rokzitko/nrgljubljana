#include <gtest/gtest.h>

#include <kk/kk.hpp>

TEST(KK, empty_input_throws) { // NOLINT
  auto construct = [] { return NRG::KK::KK(NRG::KK::XYFUNC{}); };
  EXPECT_THROW(construct(), std::runtime_error);
}

TEST(KK, invalid_akima_input_throws) { // NOLINT
  auto construct = [] {
    return NRG::KK::KK(NRG::KK::XYFUNC{{-1.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}, {1.0, 0.0}});
  };
  EXPECT_THROW(construct(), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
