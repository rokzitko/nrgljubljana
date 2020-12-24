#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "compare.hpp"

TEST(eigen, basic_vector){
    Eigen::Vector3i a(5,7,9);
    std::vector b = {5,7,9};
    compare(a,b);

    Eigen::Vector3i c;
    c << 5,7,9;
    compare(c,b);
}

TEST(eigen, matrix_ops) {
  Eigen::MatrixX<double> m(2,2);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);

  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;

  EXPECT_DOUBLE_EQ(m(0,0), 1);
  EXPECT_DOUBLE_EQ(m(0,1), 2);
  EXPECT_DOUBLE_EQ(m(1,0), 3);
  EXPECT_DOUBLE_EQ(m(1,1), 4);
}
