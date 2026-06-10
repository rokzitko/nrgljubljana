#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include <time_mem.hpp>

using namespace NRG;

TEST(Time_Mem, FormatMemoryUsage) { // NOLINT
  EXPECT_EQ(format_memory_usage(0), "0");
  EXPECT_EQ(format_memory_usage(12), "12");
  EXPECT_EQ(format_memory_usage(1234), "1'234");
  EXPECT_EQ(format_memory_usage(12434123), "12'434'123");
}

TEST(Time_Mem, MemTime) { // NOLINT
  MemTime mt;
  mt.brief_report();
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
