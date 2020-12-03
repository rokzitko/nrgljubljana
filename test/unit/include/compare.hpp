#pragma once


template<typename T1, typename T2, typename S1, typename S2>
void compare(const std::map<T1,T2> &a, const std::map<S1,S2> &b) {
  ASSERT_EQ(a.size(), b.size());
  for(auto const& [key, value] : a)
    EXPECT_EQ(b.at(key), value);
}

#ifdef _BOOST_UBLAS_MATRIX_

using namespace boost::numeric;

template<typename T1, typename T2>
void compare(const ublas::matrix<T1> &m1, const ublas::matrix<T2> &m2) {
    ASSERT_EQ(m1.size1(), m2.size1());
    ASSERT_EQ(m1.size2(), m2.size2());
    for(int i = 0; i < m1.size1(); i++)
        for(int j = 0; j < m1.size2(); j++)
            EXPECT_DOUBLE_EQ(m1(i, j), m2(i, j));
}

template<typename T1, typename T2>
void compare(const std::vector<std::vector<T1>> m1, const ublas::matrix<T2> &m2){
    ASSERT_EQ(m2.size1(), m1.size());
    ASSERT_EQ(m2.size2(), m1[0].size());
    for(int i = 0; i < m2.size1(); i++)
        for(int j = 0; j < m2.size2(); j++)
            EXPECT_DOUBLE_EQ(m2(i, j), m1[i][j]);
}

template<typename T1, typename T2>
void compare(const ublas::matrix<T2> &m2, const std::vector<std::vector<T1>> m1)
{compare(m1 ,m2);}

#ifdef EIGEN_MATRIX_H

template <typename T1, typename T2, int N, int M>
void compare(const Eigen::Matrix<T1,N,M> a, const ublas::matrix<T2> b){
  ASSERT_EQ(a.rows(), b.size1());
  ASSERT_EQ(a.cols(), b.size2());
  for(int i = 0; i < a.rows(); i++)
      for(int j = 0; j < a.cols(); j++)
          EXPECT_DOUBLE_EQ(a(i,j), b(i,j));
}

template <typename T1, typename T2, int N, int M>
void compare(const ublas::matrix<T2> b, const Eigen::Matrix<T1,N,M> a)
{compare(a, b);}

#endif //EIGEN_MATRIX_H
#endif //_BOOST_UBLAS_MATRIX_

#ifdef _BOOST_UBLAS_VECTOR_

template<typename T1, typename T2>
void compare(const std::vector<std::vector<T1>> m1, const ublas::vector<T2> &m2){
    ASSERT_EQ(m2.size1(), m1.size());
    ASSERT_EQ(m2.size2(), m1[0].size());
    for(int i = 0; i < m2.size1(); i++)
        for(int j = 0; j < m2.size2(); j++)
            EXPECT_DOUBLE_EQ(m2(i, j), m1[i][j]);
}

template<typename T1, typename T2>
void compare(const ublas::vector<T2> &m2, const std::vector<std::vector<T1>> m1)
{compare(m1, m2);}

#ifdef EIGEN_MATRIX_H

template <typename T1, typename T2, int N>
void compare(const Eigen::Matrix<T1,N,1> a, const ublas::vector<T2> b){
  ASSERT_EQ(a.size(), b.size());
  for(int i = 0; i < a.size(); i++)
    EXPECT_DOUBLE_EQ(a(i), b(i));
}

template <typename T1, typename T2, int N>
void compare(const ublas::vector<T2> b, const Eigen::Matrix<T1,N,1> a)
{compare(a, b);}

#endif //EIGEN_MATRIX_H
#endif //_BOOST_UBLAS_VECTOR_

#ifdef EIGEN_MATRIX_H

template <typename T1, typename T2, int N, int M, int K, int L>
void compare(Eigen::Matrix<T1,N,M> a, Eigen::Matrix<T2,K,L> b){
    ASSERT_EQ(a.rows(), b.rows());
    ASSERT_EQ(a.cols(), b.cols());
    for(int i = 0; i < a.rows(); i++)
        for(int j = 0; j < a.cols(); j++)
            EXPECT_DOUBLE_EQ(a(i,j), b(i,j));
}

template <typename T1, typename T2, int N>
void compare(Eigen::Matrix<T2,N,1> a, std::vector<T1> b){
    ASSERT_EQ(a.size(), b.size());
    for(int i = 0; i < a.size(); i++){
        EXPECT_DOUBLE_EQ(a(i), b[i]);
    }
}
template <typename T1, typename T2, int N>
void compare(std::vector<T1> b, Eigen::Matrix<T2,N,1> a)
{compare(a,b);}

#endif //EIGEN_MATRIX_H