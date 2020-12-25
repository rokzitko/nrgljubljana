#include <gtest/gtest.h>
#include <sstream>

#include <traits.hpp>
#include <io.hpp>

#include "compare.hpp"

using namespace NRG;

#ifdef INCL_UBLAS
TEST(io, inserters_ublas) {
  std::stringstream ss;
  ublas::matrix<double> m(2,2);
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  ss << m;
  EXPECT_EQ(ss.str(), "1 2 \n3 4 \n"); // trailing ws
}

TEST(io, read_matrix_ublas) {
  std::vector<std::vector<int>> const ref_matrix = {{31,41,53,46},{12,5,1,41},{5,2,4,7}};
  auto matrix = read_matrix_ublas("txt/matrix.txt");
  compare(ref_matrix, matrix);
}
#endif

#ifdef INCL_EIGEN
TEST(io, read_matrix_Eigen){
  Eigen::Matrix<double, 3, 4> ref_matrix;
  ref_matrix << 31,41,53,46,
                12,5,1,41,
                5,2,4,7;
  auto matrix = read_matrix_Eigen("txt/matrix.txt");
  compare(ref_matrix, matrix);
}
#endif
 
TEST(io, save_matrix){
  auto matrix = read_matrix("txt/matrix.txt");
  save_matrix("txt/matrix_temp.txt", matrix);
  auto matrix_temp = read_matrix("txt/matrix_temp.txt");
  compare(matrix, matrix_temp);
  std::remove("txt/matrix_temp.txt");
}

TEST(io, save_matrix_ublas){
  auto matrix = read_matrix_ublas("txt/matrix.txt");
  save_matrix("txt/matrix_temp.txt", matrix);
  auto matrix_temp = read_matrix_ublas("txt/matrix_temp.txt");
  compare(matrix, matrix_temp);
  std::remove("txt/matrix_temp.txt");
}

TEST(io, save_matrix_Eigen){
  auto matrix = read_matrix_Eigen("txt/matrix.txt");
  save_matrix("txt/matrix_temp.txt", matrix);
  auto matrix_temp = read_matrix_Eigen("txt/matrix_temp.txt");
  compare(matrix, matrix_temp);
  std::remove("txt/matrix_temp.txt");
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
