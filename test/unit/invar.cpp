#include <gtest/gtest.h>

#include <invar.hpp>
#include <mk_sym.hpp>
#include "test_common.hpp"

TEST(Invar, initInvar) { // NOLINT
  Params P;
  auto Sym = setup_Sym<double>(P);
  // should have Q and SS degrees of freedom
  EXPECT_EQ(Invar::invdim, 2);
  EXPECT_EQ(Invar::qntype.size(), 2);
  EXPECT_EQ(Invar::qntype[0], additive);
  EXPECT_EQ(Invar::qntype[1], additive);
  EXPECT_EQ(Invar::names.size(), 2);
  EXPECT_EQ(Invar::names["Q"], 0);
  EXPECT_EQ(Invar::names["SS"], 1);
}

TEST(Invar, InvarQS) { // NOLINT
  Params P;
  auto Sym = setup_Sym<double>(P);
  {
    Invar I(1,2);
  }
  {
    std::string str = "1 2";
    std::istringstream ss(str);
    Invar I;
    ss >> I;
    EXPECT_EQ(I, Invar(1,2));
  }
  {
    std::ostringstream ss;
    Invar I(1,2);
    ss << I;
    EXPECT_EQ(ss.str(), "1 2"s);
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.str(), "1 2"s);
    EXPECT_EQ(I.name(), "1_2"s);
  }
  {
    Invar I1(1);
    EXPECT_EQ(I1.str(), "1"s);
    Invar I2(1,2);
    EXPECT_EQ(I2.str(), "1 2"s);
    Invar I3(1,2,3);
    EXPECT_EQ(I3.str(), "1 2 3"s);
  }
  {
    Invar I1(1,1);
    Invar I2(1,2);
    EXPECT_EQ(I1 == I1, true);
    EXPECT_EQ(I1 != I2, true);
    EXPECT_EQ(I1 < I2, true);
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.getqn(0), 1);
    EXPECT_EQ(I.getqn(1), 2);
  }
  {
    Invar I1(2,2);
    Invar I2(3,2);
    I1.combine(I2);
    EXPECT_EQ(I1, Invar(5, 4));
  }
  {
    Invar I(1,2);
    I.inverse();
    EXPECT_EQ(I, Invar(-1,-2)); // both additive!
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.get("Q"), 1);
    EXPECT_EQ(I.get("SS"), 2);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
