#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <unistd.h>

#include <hilb/hilb.hpp>

TEST(Hilb, empty_dos_throws) { // NOLINT
  const auto filename = "hilb_empty_dos.dat";
  std::ofstream file(filename);

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_empty_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_THROW(NRG::Hilb::Hilb(5, argv), std::runtime_error);

  std::remove(filename);
}

TEST(Hilb, unsorted_dos_is_sorted_before_use) { // NOLINT
  const auto filename = "hilb_unsorted_dos.dat";
  {
    std::ofstream file(filename);
    file << "1 1\n-1 1\n0 1\n";
  }

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_unsorted_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_NO_THROW(NRG::Hilb::Hilb(5, argv));

  std::remove(filename);
}

TEST(Hilb, interpolator_move_assignment_releases_old_state) { // NOLINT
  auto a = NRG::Hilb::interpolator(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{1.0, 2.0, 3.0});
  auto b = NRG::Hilb::interpolator(std::vector<double>{-2.0, 0.0, 2.0}, std::vector<double>{4.0, 5.0, 6.0});
  a = std::move(b);
  EXPECT_DOUBLE_EQ(a(0.0), 5.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
