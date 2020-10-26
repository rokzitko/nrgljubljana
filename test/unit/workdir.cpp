#include <gtest/gtest.h>
#include <workdir.hpp>

using namespace NRG;

TEST(workdir, workdir) {
  {
    Workdir workdir("testdir");
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
