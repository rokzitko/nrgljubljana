#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <core.hpp>

using namespace NRG;

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
