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

TEST(numerics, generate_matrix) {
  const auto m1 = generate_matrix<double>(2,3);
  EXPECT_EQ(size1(m1), 2);
  EXPECT_EQ(size2(m1), 3);
}

TEST(numerics, zero_matrix) {
  const auto m1 = NRG::zero_matrix<double>(2,3);
  EXPECT_EQ(size1(m1), 2);
  EXPECT_EQ(size2(m1), 3);
  EXPECT_DOUBLE_EQ(m1(1,2), 0.0);
}

TEST(numerics, id_matrix) {
  const auto m1 = id_matrix<double>(2);
  EXPECT_EQ(size1(m1), 2);
  EXPECT_EQ(size2(m1), 2);
  EXPECT_DOUBLE_EQ(m1(0,0), 1.0);
  EXPECT_DOUBLE_EQ(m1(1,1), 1.0);
  EXPECT_DOUBLE_EQ(m1(0,1), 0.0);
  EXPECT_DOUBLE_EQ(m1(1,0), 0.0);
}

TEST(numerics, trans) {
  auto m = generate_matrix<double>(2,2);
  m(0,0) = m(1,1) = 0.0;
  m(0,1) = 2.0;
  m(1,0) = 3.0;
  const auto mt = trans(m);
  EXPECT_DOUBLE_EQ(mt(0,1), 3.0);
  EXPECT_DOUBLE_EQ(mt(1,0), 2.0);
}

TEST(numerics, herm) {
  auto m = generate_matrix<std::complex<double>>(2,2);
  m(0,0) = m(1,1) = 0.0;
  m(0,1) = 2.0+1.0i;
  m(1,0) = 3.0+4.0i;
  const auto mt = herm(m);
  EXPECT_EQ(mt(0,1), std::complex(3.0,-4.0));
  EXPECT_EQ(mt(1,0), std::complex(2.0,-1.0));
}

TEST(numerics, save_r) {
  auto m1 = generate_matrix<double>(10,20);
  for (size_t i = 0; i < 10*20; i++) m1(i) = (double)i;
  {
    std::ofstream f("test", std::ios::binary | std::ios::out);
    EXPECT_TRUE(f);
    boost::archive::binary_oarchive oa(f);
    NRG::save(oa, m1);
    EXPECT_FALSE(f.bad());
  }
  {
    std::ifstream f("test", std::ios::binary | std::ios::in);
    EXPECT_TRUE(f);
    boost::archive::binary_iarchive ia(f);
    const auto m2 = NRG::load<double>(ia);
    compare(m1, m2);
    EXPECT_FALSE(f.bad());
  }
}

TEST(numerics, save_c) {
  auto m1 = generate_matrix<std::complex<double>>(10,20);
  for (size_t i = 0; i < 10*20; i++) m1(i) = (double)i + 1.0i;
  {
    std::ofstream f("test", std::ios::binary | std::ios::out);
    EXPECT_TRUE(f);
    boost::archive::binary_oarchive oa(f);
    NRG::save(oa, m1);
    EXPECT_FALSE(f.bad());
  }
  {
    std::ifstream f("test", std::ios::binary | std::ios::in);
    EXPECT_TRUE(f);
    boost::archive::binary_iarchive ia(f);
    const auto m2 = NRG::load<std::complex<double>>(ia);
    compare(m1, m2);
    EXPECT_FALSE(f.bad());
  }
}

TEST(numerics, product) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto r = NRG::zero_matrix<double>(2,2);
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 2;
  product<double>(r, 1.0, a, b);
  compare(r, ref);
}
 
TEST(numerics, transform) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  auto c = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  c(0,0) = c(0,1) = c(1,0) = c(1,1) = 1;
  auto r = NRG::zero_matrix<double>(2,2);
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 4;
  transform<double>(r, 1.0, a, b, c);
  compare(r, ref);
}
 
TEST(numerics, rotate) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto r = NRG::zero_matrix<double>(2,2);
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 4;
  rotate<double>(r, 1.0, a, b);
  compare(r, ref);
}

TEST(numerics, matrix_prod) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 2;
  const auto r = matrix_prod<double>(a, b);
  compare(r, ref);
}

TEST(numerics, matrix_adj_prod) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 2;
  const auto r = matrix_adj_prod<double>(a, b);
  compare(r, ref);
}

#if defined(INCL_UBLAS) && defined(INCL_EIGEN)
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
  Eigen::MatrixXcd m_eigen_result;

  numericsMatrixOperationsTest_complex(): N(3), factor(2.5 + 1.5i), m1_ublas(N, N), m2_ublas(N, N), m3_ublas(N, N) {}

  void SetUp() override {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        m1_ublas(i,j) = 1.0 + 2.0 * (double)i + 3.0i * (double)j;
        m2_ublas(i,j) = 4.0 + 5.0 * (double)i + 6.0i * (double)j;
        m3_ublas(i,j) = 7.0 + 8.0 * (double)i + 9.0i * (double)j;
      }
    }
    m_ublas_result = ublas::zero_matrix<std::complex<double>>(N,N);
    m1_eigen = ublas_to_eigen(m1_ublas);
    m2_eigen = ublas_to_eigen(m2_ublas);
    m3_eigen = ublas_to_eigen(m3_ublas);
    m_eigen_result = Eigen::MatrixXcd::Zero(N,N);
  }
};

TEST_F(numericsMatrixOperationsTest_complex, product){
  product<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, transform){
  transform<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, rotate){
  rotate<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, product_twice) {
  product<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  product<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, transform_twice) {
  transform<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  transform<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, rotate_twice) {
  rotate<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<std::complex<double>>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  rotate<std::complex<double>>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, matrix_prod) {
  const auto m_ublas_result = matrix_prod<std::complex<double>>(m1_ublas, m2_ublas);
  const auto m_eigen_result = matrix_prod<std::complex<double>>(m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_complex, matrix_adj_prod) {
  const auto m_ublas_result = matrix_adj_prod<std::complex<double>>(m1_ublas, m2_ublas);
  const auto m_eigen_result = matrix_adj_prod<std::complex<double>>(m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

class numericsMatrixOperationsTest_real : public ::testing::Test {
protected:
  int N;
  double factor;
  ublas::matrix<double> m1_ublas;
  ublas::matrix<double> m2_ublas;
  ublas::matrix<double> m3_ublas;
  ublas::matrix<double> m_ublas_result;
  Eigen::MatrixXd m1_eigen;
  Eigen::MatrixXd m2_eigen;
  Eigen::MatrixXd m3_eigen;
  Eigen::MatrixXd m_eigen_result;

  numericsMatrixOperationsTest_real(): N(3), factor(2.5), m1_ublas(N, N), m2_ublas(N, N), m3_ublas(N, N) {}

  void SetUp() override {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        m1_ublas(i,j) = 1.0 + 2.0 * (double)i + 3.0 * (double)j;
        m2_ublas(i,j) = 4.0 + 5.0 * (double)i + 6.0 * (double)j;
        m3_ublas(i,j) = 7.0 + 8.0 * (double)i + 9.0 * (double)j;
      }
    }
    m_ublas_result = ublas::zero_matrix<double>(N,N);
    m1_eigen = ublas_to_eigen(m1_ublas);
    m2_eigen = ublas_to_eigen(m2_ublas);
    m3_eigen = ublas_to_eigen(m3_ublas);
    m_eigen_result = Eigen::MatrixXd::Zero(N,N);
  }
};

TEST_F(numericsMatrixOperationsTest_real, product){
  product<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, transform){
  transform<double>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<double>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, rotate){
  rotate<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, product_twice) {
  product<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  product<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  product<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, transform_twice) {
  transform<double>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<double>(m_ublas_result, factor, m1_ublas, m2_ublas, m3_ublas);
  transform<double>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  transform<double>(m_eigen_result, factor, m1_eigen, m2_eigen, m3_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, rotate_twice) {
  rotate<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<double>(m_ublas_result, factor, m1_ublas, m2_ublas);
  rotate<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  rotate<double>(m_eigen_result, factor, m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, matrix_prod) {
  const auto m_ublas_result = matrix_prod<double>(m1_ublas, m2_ublas);
  const auto m_eigen_result = matrix_prod<double>(m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST_F(numericsMatrixOperationsTest_real, matrix_adj_prod) {
  const auto m_ublas_result = matrix_adj_prod<double>(m1_ublas, m2_ublas);
  const auto m_eigen_result = matrix_adj_prod<double>(m1_eigen, m2_eigen);
  compare(m_ublas_result, m_eigen_result);
}

TEST(misc, ublas_to_eigen){
  const auto ublas_matrix = read_matrix_ublas("txt/matrix.txt");
  const auto eigen_matrix = ublas_to_eigen(ublas_matrix);
  compare(eigen_matrix, ublas_matrix);

  ublas::vector<int> ublas_vector(5);
  for(int i = 0; i < 5; i++) ublas_vector(i) = i;
  const auto eigen_vector = ublas_to_eigen(ublas_vector);
  compare(eigen_vector, ublas_vector);
}

TEST(misc, eigen_to_ublas){
  Eigen::Vector4i eigen_vector(3,5,8,14);
  const auto ublas_vector = eigen_to_ublas_vector(eigen_vector);
  compare(eigen_vector, ublas_vector);
  {
    Eigen::Matrix3i eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    const auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }

  {
    Eigen::Matrix<int, 3, 3, Eigen::RowMajor> eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    const auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }
}
#endif
