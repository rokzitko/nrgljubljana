#include <gtest/gtest.h>
#include <traits.hpp>
#include <io.hpp>
#include <sstream>

#include "compare.hpp"

using namespace NRG;

TEST(basicio, read_matrix){
    std::vector<std::vector<int>> const ref_matrix = {{31,41,53,46},{12,5,1,41},{5,2,4,7}};
    auto matrix = read_matrix("txt/matrix.txt");
    compare(ref_matrix, matrix);
}

#ifdef INCL_EIGEN
TEST(basicio, read_matrix_Eigen){
    Eigen::Matrix<double, 3,4> ref_matrix;
    ref_matrix << 31,41,53,46,
                  12,5,1,41,
                  5,2,4,7;
    auto matrix = read_matrix_Eigen("txt/matrix.txt");
    compare(ref_matrix, matrix);
}
#endif

#ifdef INCL_UBLAS
TEST(basicio, read_matrix_ublas){
  std::vector<std::vector<int>> const ref_matrix = {{31,41,53,46},{12,5,1,41},{5,2,4,7}};
  auto matrix = read_matrix_ublas("txt/matrix.txt");
  compare(ref_matrix, matrix);
}
#endif

TEST(basicio, save_matrix){
    auto matrix = read_matrix("txt/matrix.txt");
    save_matrix("txt/matrix_temp.txt", matrix);
    auto matrix_temp = read_matrix("txt/matrix_temp.txt");
    compare(matrix, matrix_temp);
    std::remove("txt/matrix_temp.txt");
}

TEST(basicio, save_matrix_ublas){
    auto matrix = read_matrix_ublas("txt/matrix.txt");
    save_matrix("txt/matrix_temp.txt", matrix);
    auto matrix_temp = read_matrix_ublas("txt/matrix_temp.txt");
    compare(matrix, matrix_temp);
    std::remove("txt/matrix_temp.txt");
}

TEST(basicio, save_matrix_Eigen){
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
