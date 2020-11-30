#include <gtest/gtest.h>
#include <Eigen/Dense>

template <typename T1, typename T2, int N>
void compare_vectors(Eigen::Matrix<T2,N,1> a, std::vector<T1> b){
    ASSERT_EQ(a.size(), b.size());
    for(int i = 0; i < a.size(); i++){
        EXPECT_EQ(a(i), b[i]);
    }
}

template <typename T1, typename T2, int N, int M, int K, int L>
void compare_matrices(Eigen::Matrix<T1,N,M> a, Eigen::Matrix<T2,K,L> b){
    ASSERT_EQ(a.rows(), b.rows());
    ASSERT_EQ(a.cols(), b.cols());
    for(int i = 0; i < a.rows(); i++)
        for(int j = 0; j < a.cols(); j++)
            EXPECT_EQ(a(i,j), b(i,j));
}

TEST(eigen, basic_vector){
    Eigen::Vector3i a(5,7,9);
    std::vector b = {5,7,9};
    compare_vectors(a,b);

    Eigen::Vector3i c;
    c << 5,7,9;
    compare_vectors(c,b);

    Eigen::Vector3i d;
    d << 5;
    d << 7;
    d << 9;
    compare_vectors(d,b);
    // std::cout << c << std::endl;
    // std::cout << d << std::endl;
    // compare_matrices(c,d);
}
