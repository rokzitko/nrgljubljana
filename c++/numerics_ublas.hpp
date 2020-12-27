//Don't include directly, include numerics.hpp with define USE_UBLAS instead

#ifndef _NUMERICS_UBLAS_HPP_
#define _NUMERICS_UBLAS_HPP_

// ublas matrix & vector containers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "basicio.hpp"

template <scalar S>
[[nodiscard]] ublas::matrix<S> generate_ublas(const size_t size1, const size_t size2) {
  return ublas::matrix<S>(size1, size2);
}

// Generators
#ifdef USE_UBLAS
template <scalar S>
[[nodiscard]] ublas::matrix<S> generate_matrix(const size_t size1, const size_t size2) {
  return ublas::matrix<S>(size1, size2);
}

template<scalar S>
[[nodiscard]] ublas::matrix<S> zero_matrix(const size_t size1, const size_t size2) {
  return ublas::matrix<S>(size1, size2, 0);
}

template<scalar S>
[[nodiscard]] ublas::matrix<S> id_matrix(const size_t size) {
  return ublas::identity_matrix<S>(size);
}
#endif

template<scalar S> ublas::matrix<S> herm(const ublas::matrix<S> &m) { return ublas::herm(m); }
template<scalar S> ublas::matrix<S> trans(const ublas::matrix<S> &m) { return ublas::trans(m); }

// Access the low-level data storage in the matrix (used in diag.hpp)
template<scalar S> S * data(ublas::matrix<S> &m) { return bindings::traits::matrix_storage(m); }

template <scalar T>
void save(boost::archive::binary_oarchive &oa, const ublas::matrix<T> &m) {
  oa << size1(m) << size2(m);
  for (const auto i : range0(size1(m)))
    oa << ublas::vector<T>(ublas::matrix_row<const ublas::matrix<T>>(m, i));
}

template <scalar T>
auto load_ublas(boost::archive::binary_iarchive &ia) {
  const auto size1 = read_one<size_t>(ia);
  const auto size2 = read_one<size_t>(ia);
  auto m = ublas::matrix<T>(size1, size2);
  for (const auto i : range0(size1))
    ublas::matrix_row<ublas::matrix<T>>(m, i) = read_one<ublas::vector<T>>(ia);
  return m;
}

#ifdef USE_UBLAS
template <scalar T>
auto load(boost::archive::binary_iarchive &ia) { return load_ublas<T>(ia); }
#endif

// Read 'size' values of type T into a ublas vector<T>.
template <scalar T> auto read_ublas_vector(std::istream &F, const size_t size) {
  ublas::vector<T> vec(size);
  for (auto j = 0; j < size; j++)
    vec[j] = read_one<T>(F);
  if (F.fail()) throw std::runtime_error("read_vector() error. Input file is corrupted.");
  return vec;
}

// Read values of type T into a ublas vector<T>. 'nr' is either vector dimension or the value of maximum index
template <scalar T> auto read_ublas_vector(std::istream &F, const bool nr_is_max_index = false) {
  const auto nr = read_one<size_t>(F);
  const auto len = nr_is_max_index ? nr+1 : nr;
  return read_ublas_vector<T>(F, len);
}

// Read 'size1' x 'size2' ublas matrix of type T.
template <scalar T> auto read_ublas_matrix(std::istream &F, const size_t size1, const size_t size2) {
  ublas::matrix<T> m(size1, size2);
  for (auto j1 = 0; j1 < size1; j1++)
    for (auto j2 = 0; j2 < size2; j2++)
      m(j1, j2) = assert_isfinite( read_one<T>(F) );
  if (F.fail()) std::runtime_error("read_matrix() error. Input file is corrupted.");
  return m;
}

#ifdef USE_UBLAS
template <scalar T> ublas::matrix<T> read_matrix(std::istream &F, const size_t size1, const size_t size2) {
  return read_ublas_matrix<T>(F, size1, size2);
}
#endif
   
// M += factor * A * B^\dag
template<scalar S, ublas_matrix UM, typename t_coef = coef_traits<S>>
void product(UM &M, const t_coef factor, const UM &A, const UM &B) {
  if (finite_size(A) && finite_size(B)) { // if this contributes at all...
    my_assert(size1(M) == size1(A) && size2(A) == size2(B) && size1(B) == size2(M));
    my_assert(my_isfinite(factor));
    atlas::gemm(CblasNoTrans, CblasConjTrans, factor, A, B, t_coef(1.0), M);
  }
}

// M += factor * A * O * B^\dag
template<scalar S, ublas_matrix UM, typename t_coef = coef_traits<S>>
void transform(UM &M, const t_coef factor, const UM &A, const UM &O, const UM &B) {
  if (finite_size(A) && finite_size(B)) {
    my_assert(size1(M) == size1(A) && size2(A) == size1(O) && size2(O) == size2(B) && size1(B) == size2(M));
    my_assert(my_isfinite(factor));
    ublas::matrix<S> T(size1(O), size1(B));
    atlas::gemm(CblasNoTrans, CblasConjTrans, t_coef(1.0), O, B, t_coef(0.0), T); // T = O * B^\dag
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, A, T, t_coef(1.0), M); // M += factor * A * T
  }
}

// M += factor * U^\dag * O * U
template<scalar S, ublas_matrix UM, typename U_type,  typename t_coef = coef_traits<S>>
void rotate(UM &M, const t_coef factor, const U_type &U, const UM &O) { /// XXX: U_type
  if (finite_size(U)) {
    my_assert(size1(M) == size2(U) && size1(U) == size1(O) && size2(O) == size1(U) && size2(U) == size2(M));
    my_assert(my_isfinite(factor));
    ublas::matrix<S> T(size2(U), size2(O));
    atlas::gemm(CblasConjTrans, CblasNoTrans, t_coef(1.0), U, O, t_coef(0.0), T); // T = U^\dag * O
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, T, U, t_coef(1.0), M); // M += factor * T * U
  }
}

inline auto to_ublas_range(const std::pair<size_t,size_t> &p) { return ublas::range(p.first, p.second); }

template<scalar S>
auto submatrix_const(const ublas::matrix<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return ublas::matrix_range<const ublas::matrix<S>>(M, to_ublas_range(r1), to_ublas_range(r2));
}

template<scalar S>
auto submatrix(ublas::matrix<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return ublas::matrix_range<ublas::matrix<S>>(M, to_ublas_range(r1), to_ublas_range(r2));
}

template<scalar S>
void resize(ublasMatrix<S> &m, const size_t new_size1, const size_t new_size2) {
  my_assert(new_size1 <= size1(m) && new_size2 <= size2(m));
  m.resize(new_size1, new_size2);
}

template<scalar T, ublas_matrix U, ublas_matrix V> // U and/or V may be matrix views
auto matrix_prod(const U &A, const V &B) {
  my_assert(size2(A) == size1(B));
  auto M = ublas::matrix<T>(size1(A), size2(B));
  atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
  return M;
}

template<scalar T, ublas_matrix U, ublas_matrix V> // U and/or V may be matrix views
auto matrix_adj_prod(const U &A, const V &B) {
  my_assert(size1(A) == size1(B));
  auto M = ublas::matrix<T>(size2(A), size2(B));
  atlas::gemm(CblasConjTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
  return M;
}

#endif
