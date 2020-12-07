#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "compare.hpp"
#include "test_common.hpp"
#include <traits.hpp>
#include <eigen.hpp>

using namespace NRG;

TEST(Values, access) { // NOLINT
  Values<double> values;
  std::vector v = {1.0, 2.0, 3.0, 4.0, 5.0};
  values.copy(v);
  EXPECT_EQ(values.size(), 5);
  EXPECT_DOUBLE_EQ(values.lowest_rel(), 1.0);
  EXPECT_DOUBLE_EQ(values.rel(0), 1.0);
  EXPECT_DOUBLE_EQ(values.rel(4), 5.0);
  values.set_shift(1.0);
  EXPECT_DOUBLE_EQ(values.rel_zero(0), 0.0);
  EXPECT_DOUBLE_EQ(values.rel_zero(4), 4.0);
  values.set_scale(2.0);
  EXPECT_DOUBLE_EQ(values.abs(0), 2.0);
  EXPECT_DOUBLE_EQ(values.abs(4), 10.0);
  EXPECT_DOUBLE_EQ(values.abs_zero(0),  0.0);
  EXPECT_DOUBLE_EQ(values.abs_zero(4), 8.0);
  values.set_GS_energy(0.5);
  EXPECT_DOUBLE_EQ(values.absG(0), 1.5);
  EXPECT_DOUBLE_EQ(values.absG(4), 9.5);
  const auto vv = values.all_rel();
  VECTOR_DOUBLE_EQ(v, vv);
  const auto vz = values.all_rel_zero();
  EXPECT_DOUBLE_EQ(vz[0], 0.0);
  EXPECT_DOUBLE_EQ(vz[1], 1.0);
  EXPECT_DOUBLE_EQ(vz[2], 2.0);
  EXPECT_DOUBLE_EQ(vz[3], 3.0);
  EXPECT_DOUBLE_EQ(vz[4], 4.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
