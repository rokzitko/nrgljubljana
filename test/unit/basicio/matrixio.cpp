#include <gtest/gtest.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <traits.hpp>
#include <io.hpp>

#include "compare.hpp"

using namespace NRG;

TEST(io, read_matrix){
  Eigen::Matrix<double, 3, 4> ref_matrix;
  ref_matrix << 31,41,53,46,
                12,5,1,41,
                5,2,4,7;
  auto matrix = read_matrix("txt/matrix.txt");
  EXPECT_TRUE(matrix.isApprox(ref_matrix));
}

TEST(io, save_matrix){
  auto matrix = read_matrix("txt/matrix.txt");
  save_matrix("txt/matrix_temp.txt", matrix);
  auto matrix_temp = read_matrix("txt/matrix_temp.txt");
  EXPECT_TRUE(matrix.isApprox(matrix_temp));
  std::remove("txt/matrix_temp.txt");
}

TEST(io, read_matrix_skips_blank_and_comment_lines) {
  Eigen::Matrix<double, 3, 3> ref_matrix;
  ref_matrix << 1,2,3,
                4,5,6,
                7,8,9;
  auto matrix = read_matrix("txt/matrix_comments.txt");
  EXPECT_TRUE(matrix.isApprox(ref_matrix));
}

TEST(io, read_matrix_bin_throws_on_truncated_data) {
  const auto filename = "txt/matrix_truncated.bin";
  {
    std::ofstream file(filename, std::ios::binary | std::ios::out);
    const uint32_t dim = 2;
    const double element = 1.0;
    file.write(reinterpret_cast<const char *>(&dim), sizeof(dim));
    file.write(reinterpret_cast<const char *>(&dim), sizeof(dim));
    file.write(reinterpret_cast<const char *>(&element), sizeof(element));
  }

  EXPECT_THROW(read_matrix(filename, true), std::runtime_error);
  std::remove(filename);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
