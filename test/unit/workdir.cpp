#include <string>
using namespace std::string_literals;
#include <gtest/gtest.h>
#include <workdir.hpp>

using namespace NRG;

TEST(workdir, workdir) {
  Workdir workdir("testdir");
  EXPECT_EQ(workdir.get(), "."); // because no testdir/, this defaulted to .
  EXPECT_EQ(workdir.rhofn(1, "rho"), "./rho1"s);
  EXPECT_EQ(workdir.unitaryfn(1), "./unitary1"s);
}

TEST(workdir, dtemp) {
  Workdir workdir(".");
  EXPECT_EQ(workdir.get().size(), 8); // ./XXXXXX
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
