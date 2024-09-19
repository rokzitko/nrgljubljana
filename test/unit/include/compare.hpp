#pragma once
#include <complex>
#include <string>
#include <vector>
#include <map>
#include <type_traits>
#include <exception>

#include <traits.hpp>

using namespace NRG;

void compare(const double a, const double b) {
  EXPECT_DOUBLE_EQ(a, b);
}

void compare(const std::complex<double> a, const std::complex<double> b) {
  EXPECT_DOUBLE_EQ(a.real(), b.real());
  EXPECT_DOUBLE_EQ(a.imag(), b.imag());
}

void compare(const std::string &a, const std::string &b) {
  EXPECT_EQ(a, b);
}

template<typename U, typename V>
void VECTOR_EQ(const U &A, const V &B)
{
  EXPECT_EQ(A.size(), B.size());
  for(int i = 0; i < A.size(); i++) EXPECT_EQ(A[i], B[i]);
}
    
template<typename U, typename V> // XXX: concept matrix: size1(), size2(), (i,j) accessor
void MATRIX_EQ(const U &A, const V &B)
{
  EXPECT_EQ(size1(A), size1(B));
  EXPECT_EQ(size2(A), size2(B));
  for(int i = 0; i < size1(A); i++)
    for(int j = 0; j < size2(A); j++)
      EXPECT_EQ(A(i,j), B(i,j));
}

template<typename U, typename V> // XXX: concept vector: .size, [] accessor; constrain to double
void VECTOR_DOUBLE_EQ(const U &A, const V &B)
{
  EXPECT_EQ(A.size(), B.size());
  for(int i = 0; i < A.size(); i++) EXPECT_DOUBLE_EQ(A[i], B[i]);
}
    
template<typename U, typename V> // XXX: concept matrix: size1(), size2(), (i,j) accessor; constrain to double
void MATRIX_DOUBLE_EQ(const U &A, const V &B)
{
  EXPECT_EQ(size1(A), size1(B));
  EXPECT_EQ(size2(A), size2(B));
  for(int i = 0; i < size1(A); i++)
    for(int j = 0; j < size2(A); j++)
      EXPECT_DOUBLE_EQ(A(i,j), B(i,j));
}

template<typename T1, typename T2, typename S1, typename S2>
void compare(const std::map<T1,T2> &a, const std::map<S1,S2> &b) {
  ASSERT_EQ(a.size(), b.size());
  for(auto const& [key, value] : a)
    compare(b.at(key), value);
}

template <typename T1, typename T2, int N, int M, int K, int L>
void compare(Eigen::Matrix<T1,N,M> a, Eigen::Matrix<T2,K,L> b) {
  ASSERT_EQ(a.rows(), b.rows());
  ASSERT_EQ(a.cols(), b.cols());
  for(int i = 0; i < a.rows(); i++)
    for(int j = 0; j < a.cols(); j++)
      compare(a(i,j), b(i,j));
}

template <typename T1, typename T2, int N>
void compare(Eigen::Matrix<T2,N,1> a, std::vector<T1> b) {
  ASSERT_EQ(a.size(), b.size());
  for(int i = 0; i < a.size(); i++)
    compare(a(i), b[i]);
}

template <typename T1, typename T2, int N>
void compare(std::vector<T1> b, Eigen::Matrix<T2,N,1> a) {
  compare(a,b);
}
