#include <gtest/gtest.h>
#include <basicio.hpp>
#include <sstream>

#include "compare.hpp"

using namespace NRG;

TEST(basicio, to_string) {
  {
    auto res1 = to_string(std::complex(1.0,2.0));
    EXPECT_EQ(res1, "(1,2)");
    auto res2 = to_string(std::complex(1.5,2.5));
    EXPECT_EQ(res2, "(1.5,2.5)");
  }
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

TEST(basicio, from_string) {
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

TEST(basicio, output) {
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

TEST(basicio, read_matrix){
    std::vector<std::vector<int>> const ref_matrix = {{31,41,53,46},{12,5,1,41},{5,2,4,7}};
    auto matrix = read_matrix("txt/matrix.txt");
    compare(ref_matrix, matrix);
}

TEST(basicio, _eigen_read_matrix){
    Eigen::Matrix<double, 3,4> ref_matrix;
    ref_matrix << 31,41,53,46,
                  12,5,1,41,
                  5,2,4,7;
    auto matrix = _eigen_read_matrix("txt/matrix.txt");
    compare(ref_matrix, matrix);
}

TEST(basicio, save_matrix){
    auto matrix = read_matrix("txt/matrix.txt");
    save_matrix("txt/matrix_temp.txt", matrix);
    auto matrix_temp = read_matrix("txt/matrix_temp.txt");
    compare(matrix, matrix_temp);
    std::remove("txt/matrix_temp.txt");
}

/*
TEST(basicio, _eigen_save_matrix){
    auto matrix = _eigen_read_matrix("txt/matrix.txt");
    save_matrix("txt/matrix_temp.txt", matrix); // generic function!
    auto matrix_temp = _eigen_read_matrix("txt/matrix_temp.txt");
    compare(matrix, matrix_temp);
    std::remove("txt/matrix_temp.txt");
}
*/

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
