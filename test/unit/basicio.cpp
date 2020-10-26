#include <gtest/gtest.h>
#include <basicio.hpp>
#include <sstream>

using namespace NRG;

TEST(to_string, basicio) {
  {
    auto res = to_string(std::complex(1.0,2.0));
    EXPECT_EQ(res, "(1.0,2.0)");
  }
}

TEST(inserters, basicio) {
  {
    auto res = to_string(std::make_pair(1,2));
    EXPECT_EQ(res, "1 2");
  }
  {
    auto res = to_string(std::vector<int>({1,2,3,4}));
    EXPECT_EQ(res, "1 2 3 4");
  }
  {
    ublas::vector<int> a(2);
    a(0) = 1;
    a(1) = 2;
    auto res = to_string(a);
    EXPECT_EQ(res, "1 2");
  }
  {
    ublas::vector<int> m(2,2);
    m(0,0) = 1;
    m(0,1) = 2;
    m(1,0) = 3;
    m(1,1) = 4;
    auto res = to_string(m);
    EXPECT_EQ(res, "1 2\n3 4\n");
  }
}

TEST(from_string, basicio) {
  {
    auto str = "12.3";
    auto res = from_string<double>(str);
    EXPECT_EQ(res, 12.3);
  }
  {
    auto str1 = "true";
    auto res1 = from_string<bool>(str1);
    EXPECT_EQ(res1, true);
    auto str2 = "TRUE";
    auto res2 = from_string<bool>(str2);
    EXPECT_EQ(res2, true);
    auto str3 = "NOTTRUE";
    auto res3 = from_string<bool>(str3);
    EXPECT_EQ(res3, false);
  }
}

TEST(output, basicio) {
  {
    set<> s = {1, 2, 3);
    stringstream ss;
    ss << s;
    EXPECT_EQ(ss.str(), "1 2 3");
  }
  {
    auto res = prec(1.23456, 1);
    EXPECT_EQ(res, "1.2");
    auto res3 = prec3(1.23456);
    EXPECT_EQ(res3, "1.234");
  }
  {
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,0.1)), false);
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14)), true);
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14,1e-16)), false);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
