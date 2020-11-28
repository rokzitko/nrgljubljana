#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <params.hpp>
#include <diag.hpp>

using namespace NRG;

TEST(Diag, diagonalise) { // NOLINT
  auto P = test_setup_P();
  auto DP = DiagParams(P, 0.5);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
