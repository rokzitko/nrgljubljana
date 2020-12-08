#include <gtest/gtest.h>
#include <numerics.hpp>
#include "compare.hpp"

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
TEST(numerics, Zero_matrix){
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  auto zero_m1 = Zero_matrix<double>(dim1,dim2);
  auto zero_m2 = Zero_matrix<double>(dim3);

  ASSERT_EQ(zero_m1.size1(), dim1);
  ASSERT_EQ(zero_m1.size2(), dim2);
  ASSERT_EQ(zero_m2.size1(), dim3);
  ASSERT_EQ(zero_m2.size2(), dim3);

  for(size_t i = 0; i < dim1; i++){
    for(size_t j = 0; j < dim2; j++){
      EXPECT_EQ(zero_m1(i,j), 0);
    }
  }

  for(size_t i = 0; i < dim3; i++){
    for(size_t j = 0; j < dim3; j++){
      EXPECT_EQ(zero_m2(i,j), 0);
    }
  }

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

TEST(numerics, ublas_trace_exp_real){
  const int N = 3;
  ublas::vector<double> v(N);
  for(int i = 0; i < N ; i++)
    v(i) = i + 1;
  
  ublas::matrix<double> m(N, N);
  for(int i = 0; i < N*N; i++)
    m(i) = i + 1;
  
  double expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v(i)) * m(i, i);
  }

  EXPECT_DOUBLE_EQ(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, ublas_trace_exp_complex){
  const int N = 3;
  ublas::vector<std::complex<double>> v(N);
  for(int j = 0; j < N ; j++)
    v(j) = j + 1 + 3i * j;
  
  ublas::matrix<std::complex<double>> m(N, N);
  for(int j = 0; j < N*N; j++)
    m(j) = j + 3 + 2i * j;
  
  std::complex<double> expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v(i)) * m(i, i);
  }

  compare(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, Eigen_trace_exp_real){
  const int N = 3;
  Eigen::VectorXd v(N);
  for(int i = 0; i < N ; i++)
    v(i) = i + 1;
  
  Eigen::MatrixXd m(N, N);
  for(int i = 0; i < N*N; i++)
    m(i) = i + 1;
  
  double expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v[i]) * m(i, i);
  }

  EXPECT_DOUBLE_EQ(expected, trace_exp(v, m , 2.5));
}

TEST(numerics, Eigen_trace_exp_complex){
  const int N = 3;
  Eigen::VectorXcd v(N);
  for(int j = 0; j < N ; j++)
    v(j) = j + 1 + 3i * j;
  
  Eigen::MatrixXcd m(N, N);
  for(int j = 0; j < N*N; j++)
    m(j) = j + 3 + 2i * j;
  
  std::complex<double> expected = 0;
  for(int i = 0; i < N; i++){
    expected += exp(- 2.5 * v(i)) * m(i, i);
  }

  compare(expected, trace_exp(v, m , 2.5));
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
      m1_ublas(j) = j + 3 + 2i * j;
    
    for(int j = 0; j < N*N; j++)
      m2_ublas(j) = j + 1 + 3i * j;

    for(int j = 0; j < N*N; j++)
      m3_ublas(j) = 2 * j + 1 + 5i * j;

    m_ublas_result = ublas::zero_matrix<std::complex<double>>(N,N);
    m1_eigen = ublas_to_eigen(m1_ublas);
    m2_eigen = ublas_to_eigen(m2_ublas);
    m3_eigen = ublas_to_eigen(m3_ublas);
  }
};

// class numericsMatrixOperationsTest_real : public numericsMatrixOperationsTest_complex {
// protected:
//   double factor_real = factor.real();
//   ublas::matrix<double> m1_ublas_real = ublas::real(m1_ublas);
//   ublas::matrix<double> m2_ublas_real = ublas::real(m2_ublas);
//   ublas::matrix<double> m3_ublas_real = ublas::real(m3_ublas);
//   ublas::matrix<double> m_ublas_result_real = ublas::real(m_ublas_result);

//   Eigen::MatrixXd m1_eigen_real = m1_eigen.real();
//   Eigen::MatrixXd m2_eigen_real = m2_eigen.real();
//   Eigen::MatrixXd m3_eigen_real = m3_eigen.real();
// };

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

TEST(numerics, submatrix_ublas){
  const int N = 3; 
  ublas::matrix<double> a(N, N);
  for(int i = 0; i < N*N; i++) a(i) = i;
  
  ublas::matrix<double> dg(N - 1, N - 1);
  dg(0, 0) = 0;
  dg(0, 1) = 1;
  dg(1, 0) = 3;
  dg(1, 1) = 4;
  ublas::matrix<double> dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics, submatrix_eigen){
  const int N = 3; 
  Eigen::MatrixXd a(N,N);
  for(int i = 0; i < N*N; i++) a(i) = i;

  Eigen::MatrixXd dg(N - 1, N - 1);
  dg << 0, 3, 1, 4;

  Eigen::MatrixXd dg_result = submatrix(a, {0, 2}, {0, 2});
  compare(dg, dg_result);
}

TEST(numerics, sum_of_exp){
  const int N = 3;
  double factor = 1.8;
  ublas::matrix<double> a(N, N);
  for(int i = 0; i < N*N; i++) a(i) = i;
  double result_ublas = sum_of_exp(a.data(), factor);
  double sum = 0;
  for(int i = 0; i < N*N; i++) sum += exp(-factor*a(i));
  compare(result_ublas, sum);

  Eigen::MatrixXd b(N, N);
  for(int i = 0; i < N*N; i++) b(i) = i;
  double result_eigen = sum_of_exp(b, factor);
  compare(result_eigen, sum);

}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
