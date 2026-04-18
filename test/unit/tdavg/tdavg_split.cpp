#include <gtest/gtest.h>

#define main tdavg_program_main
#include <tdavg/tdavg.cc>
#undef main

TEST(TDAvg, split_keeps_high_bit_bytes) { // NOLINT
  const std::string input = std::string(1, static_cast<char>(0x80)) + "1 2";
  const auto columns = split(input);

  ASSERT_EQ(columns.size(), 2U);
  EXPECT_DOUBLE_EQ(columns[0], 0.0);
  EXPECT_DOUBLE_EQ(columns[1], 2.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
