//Don't include directly, include numerics.hpp with define USE_EIGEN instead

#ifndef _NUMERICS_EIGEN_HPP_
#define _NUMERICS_EIGEN_HPP_

template <scalar S> 
[[nodiscard]] Eigen::MatrixX<S> generate_Eigen(const size_t size1, const size_t size2) {
  return Eigen::MatrixX<S>(size1, size2);
}

// Generators
#ifdef USE_EIGEN
template <scalar S> 
[[nodiscard]] Eigen::MatrixX<S> generate_matrix(const size_t size1, const size_t size2) {
  return Eigen::MatrixX<S>(size1, size2);
}

template <scalar S>
[[nodiscard]] Eigen::MatrixX<S> zero_matrix(const size_t size1, const size_t size2) {
  return Eigen::MatrixX<S>::Zero(size1, size2);
}

template <scalar S>
[[nodiscard]] Eigen::MatrixX<S> id_matrix(const size_t size) { 
  return Eigen::MatrixX<S>::Identity(size, size);
}
#endif

// XXX: return view instead?
template <scalar S> Eigen::MatrixX<S> herm(const Eigen::MatrixX<S> &m) { return m.adjoint(); }
template <scalar S> Eigen::MatrixX<S> trans(const Eigen::MatrixX<S> &m) { return m.transpose(); }

// Access the low-level data storage in the matrix (used in diag.hpp)
template<scalar S> S * data(Eigen::MatrixX<S> &m) { return m.data(); }

template <scalar T>
void save(boost::archive::binary_oarchive &oa, const Eigen::MatrixX<T> &m) {
  oa << size1(m) << size2(m);
  for (const auto row : m.rowwise())
    oa << row;
}

template <scalar T>
auto load_Eigen(boost::archive::binary_iarchive &ia) {
  const auto size1 = read_one<size_t>(ia);
  const auto size2 = read_one<size_t>(ia);
  auto m = Eigen::MatrixX<T>(size1, size2);
  for (const auto i : range0(size1))
    m.row(i) = read_one<Eigen::MatrixX<T>>(ia);
  return m;
}

#ifdef USE_EIGEN
template <scalar T>
auto load(boost::archive::binary_iarchive &ia) { return load_Eigen<T>(ia); }
#endif

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

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void product(EM &M, const t_coef factor, const EM &A, const EM &B) {
  if (finite_size(A) && finite_size(B)) {
    my_assert(size1(M) == size1(A) && size2(A) == size2(B) && size1(B) == size2(M));   
    my_assert(my_isfinite(factor));
    M += factor * A * B.adjoint();
  }
}

template<scalar S, Eigen_matrix EM, typename t_coef = coef_traits<S>>
void transform(EM &M, const t_coef factor, const EM &A, const EM &O, const EM &B) {
  if (finite_size(A) && finite_size(B)) {
    my_assert(size1(M) == size1(A) && size2(A) == size1(O) && size2(O) == size2(B) && size1(B) == size2(M));
    my_assert(my_isfinite(factor));
    M += factor * A * O * B.adjoint();
  }
}

template<scalar S, typename U_type, Eigen_matrix EM, typename t_coef = coef_traits<S>> // XXX: U_type
void rotate(EM &M, const t_coef factor, const U_type &U, const EM &O) {
  if (finite_size(U)) {
    my_assert(size1(M) == size2(U) && size1(U) == size1(O) && size2(O) == size1(U) && size2(U) == size2(M));
    my_assert(my_isfinite(factor));
    M += factor * U.adjoint() * O * U;
  }
}

template<scalar S>
Eigen::Block<const Eigen::MatrixX<S>> submatrix(const Eigen::MatrixX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first);
}

template<scalar S>
Eigen::Block<Eigen::MatrixX<S>> submatrix(Eigen::MatrixX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first);
}

template<scalar T, Eigen_matrix U, Eigen_matrix V>
Eigen::MatrixX<T> matrix_prod(const U &A, const V &B) {
  my_assert(size2(A) == size1(B));
  return A * B;
}

template<scalar T, Eigen_matrix U, Eigen_matrix V>
Eigen::MatrixX<T> matrix_adj_prod(const U &A, const V &B) {
  my_assert(size1(A) == size1(B));
  return A.adjoint() * B;
}

#endif
