#include <gtest/gtest.h>
#include <basicio.hpp>
#include <sstream>

using namespace NRG;

TEST(to_string, basicio) {
  {
    auto res1 = to_string(std::complex(1.0,2.0));
    EXPECT_EQ(res1, "(1,2)");
    auto res2 = to_string(std::complex(1.5,2.5));
    EXPECT_EQ(res2, "(1.5,2.5)");
  }
}

TEST(inserters1, basicio) {
  {
    std::stringstream ss;
    ss << std::make_pair(1,2);
    EXPECT_EQ(ss.str(), "1 2");
  }
  {
    std::stringstream ss;
    ss << std::vector<int>({1,2,3,4});
    EXPECT_EQ(ss.str(), "1 2 3 4 "); // trailing whitespace
  }
  {
    std::stringstream ss;
    ublas::vector<int> a(2);
    a(0) = 1;
    a(1) = 2;
    ss << a;
    EXPECT_EQ(ss.str(), "1 2 "); // trailing ws
  }
  {
    std::stringstream ss;
    ublas::matrix<int> m(2,2);
    m(0,0) = 1;
    m(0,1) = 2;
    m(1,0) = 3;
    m(1,1) = 4;
    ss << m;
    EXPECT_EQ(ss.str(), "1 2 \n3 4 \n"); // trailing ws
  }
  {
    std::stringstream ss;
    std::set s = {1, 2, 3};
    ss << s;
    EXPECT_EQ(ss.str(), "1 2 3 "); // trailing ws
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
    auto res = prec(1.23456, 1);
    EXPECT_EQ(res, "1.2");
    auto res3 = prec3(1.23456);
    EXPECT_EQ(res3, "1.235"); // rounding!
  }
  {
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,0.1)), false);
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14)), true);
    EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14),1e-16), false);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
