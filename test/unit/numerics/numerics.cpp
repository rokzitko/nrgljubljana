#include <gtest/gtest.h>

// Use the matrix backend settings from traits.hpp
#include <traits.hpp>
#include <numerics.hpp>
#include <io.hpp>

#include "compare.hpp"

using namespace NRG;
using namespace std::complex_literals;

TEST(numerics, reim) {
  const auto [r, i] = reim(std::complex(1.0, 2.0));
  EXPECT_EQ(r, 1.0);
  EXPECT_EQ(i, 2.0);
}

/*
TEST(numerics, reim_non_const) {
  auto z = std::complex(1.0, 2.0);
  auto [r, i] = reim(z);
  EXPECT_EQ(r, 1.0);
  EXPECT_EQ(i, 2.0);
  r = 3.0;
  i = 4.0;
  EXPECT_EQ(r, 3.0);
  EXPECT_EQ(i, 4.0);
  EXPECT_EQ(z.real(), 3.0);
  EXPECT_EQ(z.imag(), 4.0);
}
 */

TEST(numerics, sum2){
	std::vector v = {std::pair{4,12}, std::pair{51,123}, std::pair{412,441}};
	EXPECT_EQ(sum2(v), 576);

	std::vector<std::pair<int,int>> empty_vec = {};
	EXPECT_EQ(sum2(empty_vec), 0);
}

TEST(numerics, conj_me){
  EXPECT_EQ(conj_me(std::complex(5.0,8.0)), std::complex(5.0,-8.0));
  EXPECT_EQ(conj_me(5.0), 5.0);
}

TEST(numerica, empty_matrix) {
  const auto m = empty_matrix<double>();
  EXPECT_EQ(size1(m), 0);
  EXPECT_EQ(size2(m), 0);
}

TEST(numerics, finite_size) {
  const auto m = empty_matrix<double>();
  EXPECT_FALSE(finite_size(m));
  const auto f = generate_matrix<double>(1,1);
  EXPECT_TRUE(finite_size(f));
}

TEST(numerics, has_lesseq_rows) {
  const auto a = generate_matrix<double>(2,4);
  const auto b = generate_matrix<double>(3,4);
  EXPECT_TRUE(has_lesseq_rows(a,b));
  EXPECT_FALSE(has_lesseq_rows(b,a));
  const auto c = generate_matrix<double>(2,3);
  EXPECT_FALSE(has_lesseq_rows(a,c));
}

TEST(numerics, is_square) {
  const auto sq = generate_matrix<double>(2,2);
  EXPECT_TRUE(is_square(sq));
  const auto re = generate_matrix<double>(2,3);
  EXPECT_FALSE(is_square(re));
}

TEST(numerics, is_real) {
  EXPECT_TRUE(is_real(2.0));
  EXPECT_TRUE(is_real(std::complex(2.0,0.0)));
  EXPECT_FALSE(is_real(std::complex(2.0,1.0)));
}

TEST(numerics, check_is_matrix_upper) {
  auto m = generate_matrix<double>(2,2);
  m(0,0) = m(0,1) = m(1,0) = m(1,1) = 1.0;
  EXPECT_FALSE(is_matrix_upper(m));
  m(1,0) = 0.0;
  EXPECT_TRUE(is_matrix_upper(m));
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

TEST(numerics, num_equal) {
  const double a = 1.0;
  const double b = a + 1e-4;
  
  EXPECT_FALSE(num_equal(a,b));
  EXPECT_TRUE(num_equal(a,b,1e-2));
}

TEST(numerics, chop) {
  EXPECT_DOUBLE_EQ(chop(1e-10), 0.0);
  EXPECT_DOUBLE_EQ(chop(1.0), 1.0);
}

TEST(numerics, Power) {
  EXPECT_EQ(Power(2, 0), 1);
  EXPECT_EQ(Power(2, 1), 2);
  EXPECT_EQ(Power(2, 2), 4);
  EXPECT_EQ(Power(2, 3), 8);
  EXPECT_DOUBLE_EQ(Power(2, 0.5), std::sqrt(2.0));
}

TEST(numerics, are_conjugate_re) {
  const double x = 1.0;
  EXPECT_TRUE(are_conjugate(x,x));
}

TEST(numerics, are_conjugate_cmpl) {
  const auto z = std::complex<double>(1.0, 2.0);
  const auto zbar = std::conj(z);
  EXPECT_TRUE(are_conjugate(z,zbar));
  EXPECT_FALSE(are_conjugate(z,z));
}

TEST(numerics, is_unitary_real) {
  using T = double;
  
  auto i = generate_matrix<T>(2,2);
  i(0,0) = i(1,1) = 1.0;
  i(0,1) = i(1,0) = 0.0;
  EXPECT_TRUE(is_unitary<T>(i));

  auto z = generate_matrix<T>(2,2);
  z(0,0) = z(0,1) = z(1,0) = z(1,1) = 0.0;
  EXPECT_FALSE(is_unitary<T>(z));

  auto m = generate_matrix<T>(2,2);
  m(0,0) = 1.0;
  m(0,1) = 1.0/sqrt(2.0);
  m(1,0) = 0.0;
  m(1,1) = 1.0/sqrt(2.0);
  EXPECT_FALSE(is_unitary<T>(m));

  auto t = generate_matrix<T>(2,2);
  t(0,0) = 1.0;
  t(0,1) = 0.0;
  t(1,0) = 1.0/sqrt(2.0);
  t(1,1) = 1.0/sqrt(2.0);
  EXPECT_FALSE(is_unitary<T>(t));
}

TEST(numerics, psgn) {
  EXPECT_DOUBLE_EQ(psgn(0), 1.0);
  EXPECT_DOUBLE_EQ(psgn(1), -1.0);
}

TEST(numerics, frobenius_norm_1) {
  auto m1 = generate_matrix<double>(1,1);
  m1(0,0) = 1;
  EXPECT_DOUBLE_EQ(frobenius_norm(m1), 1.0);
}

TEST(numerics, frobenius_norm_2) {
  auto m2 = generate_matrix<double>(2,2);
  m2(0,0) = m2(0,1) = m2(1,0) = m2(1,1) = 1.0;
  EXPECT_DOUBLE_EQ(frobenius_norm(m2), 4.0);
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

TEST(numerics, dump_matrix) {
  const auto m2 = id_matrix<double>(2);
  std::ostringstream str;
  dump_matrix(m2, str);
  // 
}

TEST(numerics, dump_diagonal_matrix) {
  const auto m2 = id_matrix<double>(2);
  std::ostringstream str;
  dump_diagonal_matrix(m2, size1(m2), str);
  EXPECT_EQ(str.str(), "1 1 \n");
}

TEST(numerics, csqrt) {
  const auto z = std::complex<double>(4.0, 0.0);
  const auto r = csqrt(z);
  EXPECT_DOUBLE_EQ(r.real(), 2.0);
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

TEST(numerics, trace_real) {
  auto m = generate_matrix<std::complex<double>>(2,2);
  m(0,0) = 1.0;
  m(1,1) = 2.0;
  EXPECT_DOUBLE_EQ(trace_real(m), 3.0);
}

TEST(numerics, trace_contract) {
  auto a = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1.0;
  EXPECT_DOUBLE_EQ(trace_contract(a, a, size1(a)), 4.0);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;
  EXPECT_DOUBLE_EQ(trace_contract(a, a, size1(a)), 29.0);
}

TEST(numerics, trace_exp) {
  std::vector<double> v(2);
  v[0] = 1.0;
  v[1] = 2.0;
  auto m = generate_matrix<double>(2,2);
  m(0,0) = m(1,1) = 2.0;
  EXPECT_DOUBLE_EQ(trace_exp(v, m, 2.0), 2.0*(exp(-2.0)+exp(-4.0)));
}

TEST(numerics, sum_of_exp) {
  std::vector<double> v(2);
  v[0] = 1.0;
  v[1] = 2.0;
  EXPECT_DOUBLE_EQ(sum_of_exp(v, 2.0), exp(-2.0) + exp(-4.0));
}

TEST(numerics, trim_matrix) {
  auto a = generate_matrix<double>(2,2);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;

  const auto b22 = trim_matrix(a, 2, 2);
  EXPECT_EQ(size1(b22), 2);
  EXPECT_EQ(size2(b22), 2);
  EXPECT_DOUBLE_EQ(a(0,0), b22(0,0));
  EXPECT_DOUBLE_EQ(a(0,1), b22(0,1));
  EXPECT_DOUBLE_EQ(a(1,0), b22(1,0));
  EXPECT_DOUBLE_EQ(a(1,1), b22(1,1));

  const auto b12 = trim_matrix(a, 1, 2);
  EXPECT_EQ(size1(b12), 1);
  EXPECT_EQ(size2(b12), 2);
  EXPECT_DOUBLE_EQ(a(0,0), b12(0,0));
  EXPECT_DOUBLE_EQ(a(0,1), b12(0,1));

  const auto b21 = trim_matrix(a, 2, 1);
  EXPECT_EQ(size1(b21), 2);
  EXPECT_EQ(size2(b21), 1);
  EXPECT_DOUBLE_EQ(a(0,0), b21(0,0));
  EXPECT_DOUBLE_EQ(a(1,0), b21(1,0));

  const auto b11 = trim_matrix(a, 1, 1);
  EXPECT_EQ(size1(b11), 1);
  EXPECT_EQ(size2(b11), 1);
  EXPECT_DOUBLE_EQ(a(0,0), b11(0,0));
}

TEST(numerics, prod_fit_left) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(1,1);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 2.0;
  b(0,0) = 2.0;
  const auto m = prod_fit_left(a, b);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 1);
  EXPECT_DOUBLE_EQ(m(0,0), 4.0);
  EXPECT_DOUBLE_EQ(m(1,0), 4.0);
}

TEST(numerics, prod_fit_left_22) {
  auto a = generate_matrix<double>(2,2);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;
  const auto m = prod_fit_left(a, a);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0), 7.0);
  EXPECT_DOUBLE_EQ(m(0,1), 10.0);
  EXPECT_DOUBLE_EQ(m(1,0), 15.0);
  EXPECT_DOUBLE_EQ(m(1,1), 22.0);
}

TEST(numerics, prod_fit_right) {
  auto a = generate_matrix<double>(1,1);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = 2.0;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 2.0;
  const auto m = prod_fit_right(a, b);
  EXPECT_EQ(size1(m), 1);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0), 4.0);
  EXPECT_DOUBLE_EQ(m(0,1), 4.0);
}

TEST(numerics, prod_fit_right_22) {
  auto a = generate_matrix<double>(2,2);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;
  const auto m = prod_fit_right(a, a);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0), 7.0);
  EXPECT_DOUBLE_EQ(m(0,1), 10.0);
  EXPECT_DOUBLE_EQ(m(1,0), 15.0);
  EXPECT_DOUBLE_EQ(m(1,1), 22.0);
}

TEST(numerics, prod_fit) {
  auto a = generate_matrix<double>(2,2);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;
  const auto m = prod_fit(a, a);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0), 7.0);
  EXPECT_DOUBLE_EQ(m(0,1), 10.0);
  EXPECT_DOUBLE_EQ(m(1,0), 15.0);
  EXPECT_DOUBLE_EQ(m(1,1), 22.0);
}

TEST(numerics, prod_adj_fit_left) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(1,1);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 2.0;
  b(0,0) = 2.0;
  const auto m = prod_adj_fit_left(a, b);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 1);
  EXPECT_DOUBLE_EQ(m(0,0), 4.0);
  EXPECT_DOUBLE_EQ(m(1,0), 4.0);
}

TEST(numerics, prod_adj_fit_left_re) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = 1.0;
  a(0,1) = 2.0;
  a(1,0) = 3.0;
  a(1,1) = 4.0;
  b(0,0) = 5.0;
  b(0,1) = 6.0;
  b(1,0) = 7.0;
  b(1,1) = 8.0;
  const auto m = prod_adj_fit_left(a, b);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0), 26.0);
  EXPECT_DOUBLE_EQ(m(0,1), 30.0);
  EXPECT_DOUBLE_EQ(m(1,0), 38.0);
  EXPECT_DOUBLE_EQ(m(1,1), 44.0);
}

TEST(numerics, prod_adj_fit_left_cx) {
  using T = std::complex<double>;
  auto a = generate_matrix<T>(2,2);
  auto b = generate_matrix<T>(2,2);
  a(0,0) = T(0.0,1.0);
  a(0,1) = T(0.0,2.0);
  a(1,0) = T(0.0,3.0);
  a(1,1) = T(0.0,4.0);
  b(0,0) = 5.0;
  b(0,1) = 6.0;
  b(1,0) = 7.0;
  b(1,1) = 8.0;
  const auto m = prod_adj_fit_left(a, b);
  EXPECT_EQ(size1(m), 2);
  EXPECT_EQ(size2(m), 2);
  EXPECT_DOUBLE_EQ(m(0,0).imag(), -26.0);
  EXPECT_DOUBLE_EQ(m(0,1).imag(), -30.0);
  EXPECT_DOUBLE_EQ(m(1,0).imag(), -38.0);
  EXPECT_DOUBLE_EQ(m(1,1).imag(), -44.0);
}

TEST(numerics, chit_weight) {
  EXPECT_DOUBLE_EQ(chit_weight(2.0, 1.0, 1.0), exp(-1.0)-exp(-2.0));
  EXPECT_DOUBLE_EQ(chit_weight(2.0, 1.0, 0.5), (exp(-0.5)-exp(-1.0))/0.5);
  EXPECT_DOUBLE_EQ(chit_weight(1.0, 1.0, 1.0), exp(-1.0));
  EXPECT_DOUBLE_EQ(chit_weight(3.0, 1.0, 1.0), (exp(-1.0)-exp(-3.0))/2.0);
  EXPECT_DOUBLE_EQ(chit_weight(3.0, 3.0, 1.0), exp(-3.0));
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
    EXPECT_TRUE(m1.isApprox(m2));
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
    EXPECT_TRUE(m1.isApprox(m2));
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
  EXPECT_TRUE(r.isApprox(ref));
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
  EXPECT_TRUE(r.isApprox(ref));
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
  EXPECT_TRUE(r.isApprox(ref));
}

TEST(numerics, matrix_prod) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 2;
  const auto r = matrix_prod<double>(a, b);
  EXPECT_TRUE(r.isApprox(ref));
}

TEST(numerics, matrix_adj_prod) {
  auto a = generate_matrix<double>(2,2);
  auto b = generate_matrix<double>(2,2);
  a(0,0) = a(0,1) = a(1,0) = a(1,1) = 1;
  b(0,0) = b(0,1) = b(1,0) = b(1,1) = 1;
  auto ref = generate_matrix<double>(2,2);
  ref(0,0) = ref(0,1) = ref(1,0) = ref(1,1) = 2;
  const auto r = matrix_adj_prod<double>(a, b);
  EXPECT_TRUE(r.isApprox(ref));
}
