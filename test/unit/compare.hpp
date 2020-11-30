
template<typename T1, typename T2, typename S1, typename S2>
void compare(const std::map<T1,T2> &a, const std::map<S1,S2> &b) {
  ASSERT_EQ(a.size(), b.size());
  for(auto const& [key, value] : a)
    EXPECT_EQ(b.at(key), value);
}

#ifdef _BOOST_UBLAS_MATRIX_
template<typename T1, typename T2>
void compare(const ublas::matrix<T1> &m1, const ublas::matrix<T2> &m2) {
    ASSERT_EQ(m1.size1(), m2.size1());
    ASSERT_EQ(m1.size2(), m2.size2());
    for(int i = 0; i < m1.size1(); i++)
        for(int j = 0; j < m1.size2(); j++)
            EXPECT_EQ(m1(i, j), m2(i, j));
}
#endif
#ifdef _BOOST_UBLAS_VECTOR_
template<typename T1, typename T2>
void compare(const std::vector<std::vector<T1>> m1, const ublas::matrix<T2> &m2){
    ASSERT_EQ(m2.size1(), m1.size());
    ASSERT_EQ(m2.size2(), m1[0].size());
    for(int i = 0; i < m2.size1(); i++)
        for(int j = 0; j < m2.size2(); j++)
            EXPECT_EQ(m2(i, j), m1[i][j]);
}
#endif
#ifdef EIGEN_MATRIX_H

#endif