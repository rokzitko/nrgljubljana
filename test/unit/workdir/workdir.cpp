#include <string>
using namespace std::string_literals;
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <ostream>
#include <stdexcept>
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
    EXPECT_TRUE(workdir.remove_workdir());
    EXPECT_FALSE(std::filesystem::exists(path));
  }
  EXPECT_FALSE(std::filesystem::exists(path));
}

TEST(workdir, temporary_workdir_is_removed_at_destruction) {
  std::string path;
  {
    Workdir workdir(".", true);
    path = workdir.get();
    std::ofstream(workdir.rhofn(1, "rho")) << "data";
    ASSERT_TRUE(std::filesystem::exists(path));
  }
  EXPECT_FALSE(std::filesystem::exists(path));
}

TEST(workdir, persistent_exact_reuses_path_and_survives) {
  Workdir test_parent(".", true);
  const auto path = test_parent.get() + "/checkpoint";
  const auto marker = path + "/marker";
  ASSERT_FALSE(std::filesystem::exists(path));

  {
    Workdir workdir(path, WorkdirMode::persistent_exact, true);
    EXPECT_EQ(workdir.get(), path);
    ASSERT_TRUE(std::filesystem::is_directory(path));
    std::ofstream(marker) << "checkpoint data";
  }
  ASSERT_TRUE(std::filesystem::exists(marker));

  {
    Workdir workdir(path, WorkdirMode::persistent_exact, true);
    EXPECT_EQ(workdir.get(), path);
    EXPECT_TRUE(std::filesystem::exists(marker));
  }
  EXPECT_TRUE(std::filesystem::exists(path));
}

TEST(workdir, persistent_exact_rejects_concurrent_owner) {
  Workdir test_parent(".", true);
  const auto path = test_parent.get() + "/checkpoint";
  Workdir owner(path, WorkdirMode::persistent_exact, true);

  EXPECT_THROW(Workdir(path, WorkdirMode::persistent_exact, true), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
