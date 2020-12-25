#include <gtest/gtest.h>
#include <sstream>

#include <traits.hpp>
#include <basicio.hpp>

#include "compare.hpp"

using namespace NRG;

TEST(basicio, to_string) {
  auto res1 = to_string(std::complex(1.0,2.0));
  EXPECT_EQ(res1, "(1,2)");
  auto res2 = to_string(std::complex(1.5,2.5));
  EXPECT_EQ(res2, "(1.5,2.5)");
}

TEST(basicio, inserters1) {
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
    std::set s = {1, 2, 3};
    ss << s;
    EXPECT_EQ(ss.str(), "1 2 3 "); // trailing ws
  } 
}

TEST(basicio, from_string) {
  auto str = "12.3";
  auto res = from_string<double>(str);
  EXPECT_EQ(res, 12.3);
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

TEST(basicio, output) {
  auto res = prec(1.23456, 1);
  EXPECT_EQ(res, "1.2");
  auto res3 = prec3(1.23456);
  EXPECT_EQ(res3, "1.235"); // rounding!
  EXPECT_EQ(negligible_imag_part(std::complex(1.0,0.1)), false);
  EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14)), true);
  EXPECT_EQ(negligible_imag_part(std::complex(1.0,1e-14),1e-16), false);
}

TEST(basicio, count_words_in_string) {
  EXPECT_EQ(count_words_in_string("one two three four a"s), 5);
  EXPECT_EQ(count_words_in_string(""s), 0);
}

TEST(basicio, get_dims) {
  auto file = safe_open_for_reading("txt/matrix.txt");
  auto const [dim1, dim2] = get_dims(file);
  EXPECT_EQ(dim1, 3);
  EXPECT_EQ(dim2, 4);
  auto file_err = safe_open_for_reading("txt/matrix_err.txt");
  EXPECT_THROW(get_dims(file_err), std::runtime_error);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
