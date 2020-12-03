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
