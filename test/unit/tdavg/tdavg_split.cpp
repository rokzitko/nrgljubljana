#include <gtest/gtest.h>

#define main tdavg_program_main
#include <tdavg/tdavg.cc>
#undef main

TEST(TDAvg, split_keeps_high_bit_bytes) { // NOLINT
  const std::string input_ = std::string(1, static_cast<char>(0x80)) + "1 2";
  const auto columns_ = split(input_);

  ASSERT_EQ(columns_.size(), 2U);
  EXPECT_DOUBLE_EQ(columns_[0], 0.0);
  EXPECT_DOUBLE_EQ(columns_[1], 2.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
