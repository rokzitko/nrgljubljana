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

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
