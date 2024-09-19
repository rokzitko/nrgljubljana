#include <gtest/gtest.h>
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
  compare(ref_matrix, matrix);
}

TEST(io, save_matrix){
  auto matrix = read_matrix("txt/matrix.txt");
  save_matrix("txt/matrix_temp.txt", matrix);
  auto matrix_temp = read_matrix("txt/matrix_temp.txt");
  compare(matrix, matrix_temp);
  std::remove("txt/matrix_temp.txt");
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
