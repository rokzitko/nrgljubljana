#include <gtest/gtest.h>
#include "compare.hpp"

// Use the matrix backend settings from traits.hpp
#include "traits.hpp"
#include "numerics.hpp"
#include "io.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics, reim) {
  const auto [r, i] = reim(std::complex(1.0,2.0));
  EXPECT_EQ(r, 1.0);
  EXPECT_EQ(i, 2.0);
}

TEST(numerics, sum2){
	std::vector v = {std::pair{4,12}, std::pair{51,123}, std::pair{412,441}};
	EXPECT_EQ(sum2(v), 576);

	std::vector<std::pair<int,int>> empty_vec = {};
	EXPECT_EQ(sum2(empty_vec), 0);
}

TEST(numerics, conj_me){
  EXPECT_EQ(conj_me(std::complex(5.0,8.0)), std::complex(5.0,-8.0));
}

TEST(numerics, bucket){
  std::vector v = {std::pair{4.0,12.0}, std::pair{51.0,123.0}, std::pair{412.0,441.0}};
  bucket sum(v);
  EXPECT_EQ(sum, 576);
  sum += 10;
  EXPECT_EQ(sum, 586);
}

TEST(numerics, is_even_is_odd){
  const int odd = 3;
  const int even = 4;
  EXPECT_TRUE(is_even(even));
  EXPECT_FALSE(is_odd(even));
  EXPECT_TRUE(is_odd(odd));
  EXPECT_FALSE(is_even(odd));
}

TEST(numerics, my_fcmp){
  const double a = 1.0;
  const double b = 1.0001;
  EXPECT_EQ(my_fcmp(a, b, 0.1, 0.00015), 0);
  EXPECT_EQ(my_fcmp(a, b, 0.1, 0.00001), -1);
}

class numericsMatrixOperationsTest_complex : public ::testing::Test {
protected:
  int N;
  std::complex<double> factor;
  ublas::matrix<std::complex<double>> m1_ublas;
  ublas::matrix<std::complex<double>> m2_ublas;
  ublas::matrix<std::complex<double>> m3_ublas;
  ublas::matrix<std::complex<double>> m_ublas_result;
  Eigen::MatrixXcd m1_eigen;
  Eigen::MatrixXcd m2_eigen;
  Eigen::MatrixXcd m3_eigen;

  numericsMatrixOperationsTest_complex(): N(3), factor(2.5 + 1.5i), m1_ublas(N, N), m2_ublas(N, N), m3_ublas(N, N) {}


  void SetUp() override{
    for(int j = 0; j < N*N; j++)
      m1_ublas(j) = j + 3.0 + 2i * (double)j;

    for(int j = 0; j < N*N; j++)
      m2_ublas(j) = j + 1.0 + 3i * (double)j;

    for(int j = 0; j < N*N; j++)
      m3_ublas(j) = 2.0 * j + 1 + 5i * (double)j;

    m_ublas_result = ublas::zero_matrix<std::complex<double>>(N,N);
    m1_eigen = ublas_to_eigen(m1_ublas);
    m2_eigen = ublas_to_eigen(m2_ublas);
    m3_eigen = ublas_to_eigen(m3_ublas);
  }
};

class numericsMatrixOperationsTest_real : public numericsMatrixOperationsTest_complex {
protected:
  double factor_real;
  ublas::matrix<double> m1_ublas_real;
  ublas::matrix<double> m2_ublas_real;
  ublas::matrix<double> m3_ublas_real;
  ublas::matrix<double> m_ublas_result_real;
  Eigen::MatrixXd m1_eigen_real;
  Eigen::MatrixXd m2_eigen_real;
  Eigen::MatrixXd m3_eigen_real;

  void SetUp_real(){
    factor_real = factor.real();
    m1_ublas_real = ublas::real(m1_ublas);
    m2_ublas_real = ublas::real(m2_ublas);
    m3_ublas_real = ublas::real(m3_ublas);
    m_ublas_result_real = ublas::real(m_ublas_result);

    m1_eigen_real = m1_eigen.real();
    m2_eigen_real = m2_eigen.real();
    m3_eigen_real = m3_eigen.real();
  }
};



TEST_F(numericsMatrixOperationsTest_real, product){
  SetUp_real();
  product<double>(m_ublas_result_real, factor_real, m1_ublas_real, m2_ublas_real);
  auto m_eigen_result_real = product<double>(factor_real, m1_eigen_real, m2_eigen_real);
  compare(m_ublas_result_real, m_eigen_result_real);
}

TEST_F(numericsMatrixOperationsTest_real, transform){
  SetUp_real();
  transform<double>(m_ublas_result_real, factor_real, m1_ublas_real, m2_ublas_real, m3_ublas_real);
  auto m_eigen_result_real = transform<double>(factor_real, m1_eigen_real, m2_eigen_real, m3_eigen_real);
  compare(m_ublas_result_real, m_eigen_result_real);
}

TEST_F(numericsMatrixOperationsTest_real, rotate){
  SetUp_real();
  rotate<double>(m_ublas_result_real, factor_real, m1_ublas_real, m2_ublas_real);
  auto m_eigen_result_real = rotate<double>(factor_real, m1_eigen_real, m2_eigen_real);
  compare(m_ublas_result_real, m_eigen_result_real);
}

TEST_F(numericsMatrixOperationsTest_complex, product){
  product<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  auto m_eigen_result = product<std::complex<double>>(factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, transform){
  transform<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  auto m_eigen_result = transform<std::complex<double>>(factor, m1_eigen, m2_eigen, m3_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, rotate){
  rotate<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  auto m_eigen_result = rotate<std::complex<double>>(factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

#if defined(INCL_UBLAS) && defined(INCL_EIGEN)
TEST(misc, ublas_to_eigen){
  auto ublas_matrix = read_matrix_ublas("txt/matrix.txt");
  auto eigen_matrix = ublas_to_eigen(ublas_matrix);
  compare(eigen_matrix, ublas_matrix);

  ublas::vector<int> ublas_vector(5);
  for(int i = 0; i < 5; i++) ublas_vector(i) = i;
  auto eigen_vector = ublas_to_eigen(ublas_vector);
  compare(eigen_vector, ublas_vector);
}

TEST(misc, eigen_to_ublas){
  Eigen::Vector4i eigen_vector(3,5,8,14);
  auto ublas_vector = eigen_to_ublas_vector(eigen_vector);
  compare(eigen_vector, ublas_vector);
  {
    Eigen::Matrix3i eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }

  {
    Eigen::Matrix<int, 3, 3, Eigen::RowMajor> eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }
}
#endif
