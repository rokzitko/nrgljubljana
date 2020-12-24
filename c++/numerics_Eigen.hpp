//Don't include directly, include numerics.hpp with define USE_EIGEN instead

#ifndef _NUMERICS_EIGEN_HPP_
#define _NUMERICS_EIGEN_HPP_

// Generators
#ifdef USE_EIGEN
template<scalar S>
[[nodiscard]] auto zero_matrix(const size_t size1, const size_t size2) {
  return Eigen::MatrixX<S>::Zero(size1, size2);
}

template<scalar S>
[[nodiscard]] auto id_matrix(const size_t size) { 
  return Eigen::MatrixX<S>::Identity(size, size); 
}
#endif

// Access the low-level data storage in the matrix (used in diag.hpp)
template<scalar S> S * data(Eigen::MatrixX<S> &m) { return m.data(); }

template <scalar T>
void save(boost::archive::binary_oarchive &oa, const Eigen::MatrixX<T> &m) {
  oa << size1(m) << size2(m);
  for (const auto row : m.rowwise())
    oa << row;
}

template <scalar T>
void load(boost::archive::binary_iarchive &ia, Eigen::MatrixX<T> &m) { // XXX
  const auto size1 = read_one<size_t>(ia);
  const auto size2 = read_one<size_t>(ia);
  m = Eigen::MatrixX<T>(size1, size2);
  for (const auto i : range0(size1))
    m.row(i) = read_one<Eigen::MatrixX<T>>(ia);
}

// Read 'size' values of type T into an Eigen vector<T>.
template <scalar T> auto read_Eigen_vector(std::istream &F, const size_t size) {
  Eigen::VectorX<T> vec(size);
  for (auto j = 0; j < size; j++)
    vec[j] = read_one<T>(F);
  if (F.fail()) throw std::runtime_error("read_vector() error. Input file is corrupted.");
  return vec;
}

// Read values of type T into an Eigen vector<T>. 'nr' is either vector dimension or the value of maximum index
template <scalar T> auto read_Eigen_vector(std::istream &F, const bool nr_is_max_index = false) {
  const auto nr = read_one<size_t>(F);
  const auto len = nr_is_max_index ? nr+1 : nr;
  return read_Eigen_vector<T>(F, len);
}

// Read 'size1' x 'size2' Eigen matrix of type T.
template <scalar T> auto read_Eigen_matrix(std::istream &F, const size_t size1, const size_t size2) {
  Eigen::MatrixX<T> m(size1, size2);
  for (auto j1 = 0; j1 < size1; j1++)
    for (auto j2 = 0; j2 < size2; j2++)
      m(j1, j2) = assert_isfinite( read_one<T>(F) );
  if (F.fail()) std::runtime_error("read_matrix() error. Input file is corrupted.");
  return m;
}

#ifdef USE_EIGEN
template <scalar T> auto read_matrix(std::istream &F, const size_t size1, const size_t size2) {
  return read_Eigen_matrix<T>(F, size1, size2);
}
#endif

template<scalar S, typename T, typename t_coef = coef_traits<S>>
Eigen::MatrixX<T> product(const t_coef factor, const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &B) {
  my_assert(my_isfinite(factor));
  return factor * A * B.adjoint();
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product(EM &M, const t_coef factor, const EM &A, const EM &B) {
  my_assert(my_isfinite(factor));
  M += factor * A * B.adjoint();
}

template<scalar S, typename t_coef = coef_traits<S>, typename T>
Eigen::MatrixX<T> transform(const t_coef factor, const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &O, const Eigen::MatrixX<T> &B) {
  my_assert(my_isfinite(factor));
  return factor * A * O * B.adjoint();
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform(EM &M, const t_coef factor, const EM &A, const EM &O, const EM &B) {
  my_assert(my_isfinite(factor));
  M += factor * A * O * B.adjoint();
}

template<scalar S, typename t_coef = coef_traits<S>, typename T, typename T1>
Eigen::MatrixX<T> rotate(const t_coef factor, const Eigen::MatrixX<T1> &U, const Eigen::MatrixX<T> &O) {
  my_assert(my_isfinite(factor));
  return factor * U.adjoint() * O * U;
}

template<scalar S, typename U_type, Eigen_matrix EM, typename t_coef = coef_traits<S>> // XXX: U_type
void rotate(EM &M, const t_coef factor, const U_type &U, const EM &O) {
  my_assert(my_isfinite(factor));
  M += factor * U.adjoint() * O * U;
}

template<scalar S>
Eigen::MatrixX<S> submatrix(const Eigen::MatrixX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first); // XXX: avoid copy?
}

template<scalar S>
Eigen::MatrixX<S> submatrix(Eigen::ArrayX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first); // XXX
}

template<typename T>
auto sum_of_exp(Eigen::MatrixX<T> A, const double factor) // sum exp(-factor*x)
{
  return exp(-factor * A.array()).sum();
}

template<typename T>
Eigen::MatrixX<T> matrix_prod(const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &B) {
  return A * B;
}

template<typename T>
Eigen::MatrixX<T> matrix_adj_prod(const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &B) {
  return A.adjoint() * B;
}

#endif
