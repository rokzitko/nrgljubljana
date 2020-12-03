#include <gtest/gtest.h>
#include <openmp.hpp>

using namespace NRG;

TEST(openMP, report) {
  report_openMP();
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
