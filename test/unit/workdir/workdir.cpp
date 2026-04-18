#include <string>
using namespace std::string_literals;
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <workdir.hpp>

using namespace NRG;

TEST(workdir, workdir) {
  Workdir workdir(".", true); // true=quiet
  EXPECT_EQ(workdir.rhofn(1, "rho"), workdir.get() + "/rho1"s);
  EXPECT_EQ(workdir.unitaryfn(1), workdir.get() + "/unitary1"s);
}

TEST(workdir, dtemp) {
  Workdir workdir(".", true); // true=quiet
  EXPECT_EQ(workdir.get().size(), 8); // ./XXXXXX
}

TEST(workdir, invalid_parent_throws) {
  EXPECT_THROW(Workdir("testdir", true), std::runtime_error);
}

TEST(workdir, remove_workdir_removes_files) {
  std::string path;
  {
    Workdir workdir(".", true);
    path = workdir.get();
    std::ofstream(workdir.rhofn(1, "rho")) << "data";
    EXPECT_TRUE(std::filesystem::exists(path));
  }
  EXPECT_FALSE(std::filesystem::exists(path));
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
