#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include <time_mem.hpp>

using namespace NRG;

TEST(Time_Mem, MemTime) { // NOLINT
  MemTime mt;
  mt.brief_report();
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
