// numerics.h - Miscelaneous numerical routines
// Copyright (C) 2005-2020 Rok Zitko

// This header should be included in all other headers where vector/matrix
// objects are manipulated.

#ifndef _numerics_hpp_
#define _numerics_hpp_

#include <complex>
#include <iomanip>
#include <vector>
#include <fstream>
#include <range/v3/all.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/math/special_functions/sign.hpp>

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

#include <Eigen/Dense>

// Serialization support (used for storing to files and for MPI)
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include "portabil.hpp"
#include "misc.hpp"
#include "traits.hpp"

namespace NRG {

using namespace boost::numeric;
using namespace boost::numeric::ublas; // keep this!
namespace atlas = boost::numeric::bindings::atlas;

template <typename T>
  using complex_array_ref_t = T(&)[2];

template<typename T>
  complex_array_ref_t<T> reim(std::complex<T>& z) {
    return reinterpret_cast<T(&)[2]>(z);
  }

template <typename T>
  using complex_array_const_ref_t = const T(&)[2];

template<typename T>
  complex_array_const_ref_t<T> reim(const std::complex<T>& z) {
    return reinterpret_cast<const T(&)[2]>(z);
  }

template<typename U, typename V>
  V sum2(const std::vector<std::pair<U,V>> &v) { // sum second elements of a vector of pairs
    return ranges::accumulate(v, V{}, {}, [](const auto el) { return el.second; }); // XXX
  }
 
[[nodiscard]] inline std::complex<double> conj_me(const std::complex<double> &z) { return conj(z); } // conjugation
[[nodiscard]] inline double conj_me(const double x) { return x; }    // no op

template<scalar S, typename Matrix = Matrix_traits<S>>
[[nodiscard]] auto zero_matrix(const size_t size1, const size_t size2) {
  return Matrix(size1, size2, 0);
}

template<scalar S>
[[nodiscard]] auto zero_matrix(const size_t size) { 
  return NRG::zero_matrix<S>(size, size); 
}

template<scalar S, typename Matrix = Matrix_traits<S>>
auto empty_matrix() { return Matrix(); }

// Access the low-level data storage in the matrix (used in diag.hpp)
template<scalar S> S * data(ublas::matrix<S> &m) { return bindings::traits::matrix_storage(m); }

// Accumulator abstraction: automatically initialized to 0, result checked for finiteness.
template <scalar T> class generic_bucket {
private:
  T value{};
public:
  generic_bucket() = default;
  // Can be constructured from a STL vector of pairs, by summing the second elements.
  template <scalar T1> 
  explicit generic_bucket(std::vector<std::pair<T1, T>> v) {
    for (const auto &i : v) value += i.second;
  }
  inline T operator+=(T x) { return value += x; }
  inline operator T() const { return value; }
};
using bucket = generic_bucket<double>;

template <typename T>
inline constexpr bool is_odd(const T n) { return n & 1; } // must return bool
template <typename T>
inline constexpr bool is_even(const T n) { return !is_odd(n); } // must return bool

inline CONSTFNC int my_fcmp(const double x, const double y, const double small_epsilon, const double rel_epsilon) {
  if (x == 0.0 && y == 0.0) return 0.0; // evidently equal
  if (std::abs(x) < small_epsilon && std::abs(y) < small_epsilon) return 0; // If both x and y are small, we ASSUME them to be equivalent
  if (std::abs(x-y) < rel_epsilon * (std::abs(x)+std::abs(y))) return 0;
  return boost::math::sign(x-y);
}

inline CONSTFNC int my_fcmp(const double x, const double y, const double epsilon) { return my_fcmp(x, y, epsilon, epsilon); }

// Test if two numbers are equal to within numerical errors. (Use this for comparing values that are expected to be
// of order 1.)
inline CONSTFNC auto num_equal(const double a, const double b, const double check_precision = 1.e-12) {
  return my_fcmp(a, b, check_precision) == 0;
}

inline CONSTFNC auto num_equal(const std::complex<double> &a, const std::complex<double> &b, const double check_precision = 1.e-12) {
  return (my_fcmp(a.real(), b.real(), check_precision) == 0) && (my_fcmp(a.imag(), b.imag(), check_precision) == 0);
}

inline CONSTFNC auto are_conjugate(const double a, const double b) { return num_equal(a, b); }

inline CONSTFNC auto are_conjugate(const std::complex<double> &a, const std::complex<double> &b) { return num_equal(a.real(), b.real()) && num_equal(a.imag(), -b.imag()); }

template<matrix M> auto frobenius_norm(const M &m) { // Frobenius norm (without taking the final square root!)
  double sum{};
  for (auto i = 0; i < size1(m); i++)
    for (auto j = 0; j < size2(m); j++) sum += pow(abs(m(i, j)),2);
  return sum;
}

template<matrix M> bool is_square(const M &m) { return size1(m) == size2(m); }

// Check if m is upper triangular. In the lower triangle, all elements must be 0.
template<matrix M> void check_is_matrix_upper(const M &m) {
  my_assert(is_square(m));
  for (auto i = 1; i < size1(m); i++)
    for (auto j = 0; j < i; j++) // j < i
      my_assert(num_equal(m(i, j), 0));
}

// (-1)^n
CONSTFNC inline auto psgn(const int n) { return n % 2 == 0 ? 1.0 : -1.0; }

// Dump a matrix with full numerical precision. The columns are aligned for easier inspection. Expect large output!
template<matrix M> inline void dump_matrix(const M &m, std::ostream &F = std::cout, 
                                           const int header_width = 7, const int column_width = 23) {
  boost::io::ios_base_all_saver ofs(F);
  F << std::setprecision(std::numeric_limits<double>::max_digits10);
  F << fmt::format("Matrix: {}x{}\n", size1(m), size2(m));
  for (auto r1 = 0; r1 < size1(m); r1++) {
    F << std::setw(header_width) << r1 << ":";
    for (auto r2 = 0; r2 < size2(m); r2++) F << std::setw(column_width) << m(r1, r2) << " ";
    F << std::endl;
  }
}

template<matrix M> inline void dump_diagonal_matrix(const M &m, const size_t max_nr, std::ostream &F = std::cout) {
  for (const auto r : range0(std::min(size1(m), max_nr))) F << m(r,r) << ' ';
  F << std::endl;
}

template <scalar T>
void save(boost::archive::binary_oarchive &oa, const ublas::matrix<T> &m) {
  oa << size1(m) << size2(m);
  for (const auto i : range0(size1(m)))
    oa << ublas::vector<T>(ublas::matrix_row<const ublas::matrix<T>>(m, i));
}

template <scalar T>
void load(boost::archive::binary_iarchive &ia, ublas::matrix<T> &m) { // XXX
  const auto size1 = read_one<size_t>(ia);
  const auto size2 = read_one<size_t>(ia);
  m = ublas::matrix<T>(size1, size2);
  for (const auto i : range0(size1))
    ublas::matrix_row<ublas::matrix<T>>(m, i) = read_one<ublas::vector<T>>(ia);
}

// Chop numerical noise
template <scalar T> CONSTFNC inline T chop(const T x, const double xlimit = 1.e-8) { return std::abs(x) < xlimit ? 0.0 : x; }

// Powers, such as (-1)^n, appear in the coupling coefficients.
CONSTFNC inline double Power(const double i, const double nn) { return std::pow(i, nn); }

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

template <scalar T> auto read_matrix(std::istream &F, const size_t size1, const size_t size2) {
  return read_ublas_matrix<T>(F, size1, size2);
}

// Check if the value x is real [for complex number calculations].
CONSTFNC inline auto is_real(const double x) { return x; }
CONSTFNC inline auto is_real(const std::complex<double> z, const double check_real_tolerance = 1e-8) {
  return abs(z.imag()) <= check_real_tolerance;
}

// Check if x is real and return the real part, i.e. x.real().
CONSTFNC inline auto check_real(double x) { return x; }
CONSTFNC inline auto check_real(std::complex<double> z) {
  if (!is_real(z)) std::cout << "Warning: expected real number, but got " << z << std::endl;
  return z.real();
}

template <matrix M> CONSTFNC auto trace_real(const M &m) {
  my_assert(is_square(m));
  return ranges::accumulate(range0(size1(m)), 0.0, {}, [&m](const auto i){ return check_real(m(i, i)); });
}

inline auto csqrt(const std::complex<double> z) { return std::sqrt(z); }

template<matrix R> // 2D matrix or matrix view
auto finite_size(const R &m) { return size1(m) && size2(m); }
   
// M += factor * A * B^\dag
template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void product(Matrix &M, const t_coef factor, const Matrix &A, const Matrix &B) {
  if (finite_size(A) && finite_size(B)) { // if this contributes at all...
    my_assert(size1(M) == size1(A) && size2(A) == size2(B) && size1(B) == size2(M));
    my_assert(my_isfinite(factor));
    atlas::gemm(CblasNoTrans, CblasConjTrans, factor, A, B, t_coef(1.0), M);
  }
}

template<scalar S, typename T, typename t_coef = coef_traits<S>>
Eigen::MatrixX<T> product(const t_coef factor, const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &B) {
  my_assert(my_isfinite(factor));
  return factor * A * B.adjoint();
}

// M += factor * A * O * B^\dag
template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void transform(Matrix &M, const t_coef factor, const Matrix &A, const Matrix &O, const Matrix &B) {
  if (finite_size(A) && finite_size(B)) {
    my_assert(size1(M) == size1(A) && size2(A) == size1(O) && size2(O) == size2(B) && size1(B) == size2(M));
    my_assert(my_isfinite(factor));
    Matrix T(size1(O), size1(B));
    atlas::gemm(CblasNoTrans, CblasConjTrans, t_coef(1.0), O, B, t_coef(0.0), T); // T = O * B^\dag
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, A, T, t_coef(1.0), M); // M += factor * A * T
  }
}

template<scalar S, typename t_coef = coef_traits<S>, typename T>
Eigen::MatrixX<T> transform(const t_coef factor, const Eigen::MatrixX<T> &A, const Eigen::MatrixX<T> &O, const Eigen::MatrixX<T> &B) {
  my_assert(my_isfinite(factor));
  return factor * A * O * B.adjoint();
}

// M += factor * U^\dag * O * U
template<scalar S, typename U_type, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void rotate(Matrix &M, const t_coef factor, const U_type &U, const Matrix &O) {
  if (finite_size(U)) {
    my_assert(size1(M) == size2(U) && size1(U) == size1(O) && size2(O) == size1(U) && size2(U) == size2(M));
    my_assert(my_isfinite(factor));
    Matrix T(size2(U), size2(O));
    atlas::gemm(CblasConjTrans, CblasNoTrans, t_coef(1.0), U, O, t_coef(0.0), T); // T = U^\dag * O
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, T, U, t_coef(1.0), M); // M += factor * T * U
  }
}

template<scalar S, typename t_coef = coef_traits<S>, typename T, typename T1>
Eigen::MatrixX<T> rotate(const t_coef factor, const Eigen::MatrixX<T1> &U, const Eigen::MatrixX<T> &O) {
  my_assert(my_isfinite(factor));
  return factor * U.adjoint() * O * U;
}

inline auto to_ublas_range(const std::pair<size_t,size_t> &p) { return ublas::range(p.first, p.second); }

template<scalar S>
auto submatrix(const ublas::matrix<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return ublas::matrix_range<const ublas::matrix<S>>(M, to_ublas_range(r1), to_ublas_range(r2));
}

template<scalar S>
auto submatrix(ublas::matrix<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return ublas::matrix_range<ublas::matrix<S>>(M, to_ublas_range(r1), to_ublas_range(r2));
}

template<scalar S>
auto submatrix(const Eigen::MatrixX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first);
}

template<scalar S>
auto submatrix(Eigen::ArrayX<S> &M, const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2)
{
  return M.block(r1.first, r2.first, r1.second - r1.first, r2.second - r2.first);
}

// 'v' is any 1D range we can iterate over
template<typename R, matrix M>
auto trace_exp(R && v, const M &m, const double factor) { // Tr[exp(-factor*v) m]
  my_assert(v.size() == size1(m) && v.size() == size2(m));
  return ranges::accumulate(range0(v.size()), typename M::value_type{}, {}, [&v, &m, factor](const auto i){ return exp(-factor * v[i]) * m(i, i); });
}

// 'values' is any 1D range we can iterate over
template<typename R>
auto sum_of_exp(R && values, const double factor) // sum exp(-factor*x)
{
  return ranges::accumulate(values, 0.0, {}, [factor](const auto &x){ return exp(-factor*x); });
}      

template<typename T>
auto sum_of_exp(Eigen::MatrixX<T> A, const double factor) // sum exp(-factor*x)
{
  return exp(-factor * A.array()).sum();
}

template<matrix M>
auto trace_contract(const M &A, const M &B, const size_t range) // Tr[AB]
{
  typename M::value_type sum{};
  for (const auto i : range0(range))
       for (const auto j : range0(range))
      sum += A(i, j) * B(j, i);
  return sum;
}

template<matrix M>
auto trim_matrix(M &mat, const size_t new_size1, const size_t new_size2) {
  const auto old_size1 = size1(mat);
  const auto old_size2 = size2(mat);
  if (old_size1 == 0 || old_size2 == 0) return;
  my_assert(new_size1 <= old_size1 && new_size2 <= old_size2);
  if (new_size1 == old_size1 && new_size2 == old_size2) return; // Trimming not necessary!!
  const auto sub = submatrix(mat, {0, new_size1}, {0, new_size2});
  M mat2 = sub;
  mat.swap(mat2);
} 

template<matrix M>
  bool has_lesseq_rows(const M &A, const M &B) {
    return size1(A) <= size1(B) && size2(A) == size2(B);
  }

template<typename T, typename U, typename V> // U and/or V may be matrix views
auto matrix_prod(const U &A, const V &B) {
  auto M = ublas::matrix<T>(size1(A), size2(B));
  atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
  return M;
}

template<typename T, typename U, typename V> // U and/or V may be matrix views
auto matrix_adj_prod(const U &A, const V &B) {
  auto M = ublas::matrix<T>(size2(A), size2(B));
  atlas::gemm(CblasConjTrans, CblasNoTrans, 1.0, A, B, 0.0, M);
  return M;
}

// M = A*B, size of A is adapted to the size of B
template<matrix M>
auto prod_fit_left(const M &A, const M &B) {
  using T = typename M::value_type;
  my_assert(size1(B) <= size2(A));
  if (size1(A) == 0 || size2(B) == 0) return empty_matrix<T>();
  const auto Asub = submatrix(A, {0, size1(A)}, {0, size1(B)});
  return matrix_prod<T>(Asub, B);
}

// M = A*B, size of B is adapted to the size of A
template<matrix M>
auto prod_fit_right(const M &A, const M &B) {
  using T = typename M::value_type;
  my_assert(size1(B) >= size2(A));
  if (size1(A) == 0 || size2(B) == 0) return empty_matrix<T>();
  const auto Bsub = submatrix(B, {0, size2(A)}, {0, size2(B)});
  return matrix_prod<T>(A, Bsub);
}

template<matrix M>
auto prod_fit(const M &A, const M &B) {
  return size1(B) <= size2(A) ? prod_fit_left(A, B) : prod_fit_right(A, B);
}

// M = A^\dag*B, size of A is adapted to the size of B
template<matrix M>
inline auto prod_adj_fit_left(const M &A, const M &B) {
  using T = typename M::value_type;
  my_assert(size1(B) <= size1(A));
  if (size2(A) == 0 || size2(B) == 0) return empty_matrix<T>();
  const auto Asub = submatrix(A, {0, size1(B)}, {0, size2(A)});
  return matrix_adj_prod<T>(Asub, B);
}

inline constexpr double WEIGHT_TOL = 1e-8; // where to switch to l'Hospital rule form

// weight=(exp(-beta Em)-exp(-beta En))/(beta En-beta Em). NOTE: arguments En, Em are order omega_N, while beta is
// order 1/omega_N, thus the combinations betaEn and betaEm are order 1. Also En>0, Em>0, since these are excitation
// energies !
inline auto chit_weight(const double En, const double Em, const double beta) {
  const auto betaEn = beta * En;
  const auto betaEm = beta * Em;
  const auto x      = betaEn - betaEm;
  if (abs(x) > WEIGHT_TOL) {
    // If one of {betaEm,betaEn} is small, one of exp() will have a value around 1, the other around 0, thus the
    // overall result will be approximately +-1/x.
    return (exp(-betaEm) - exp(-betaEn)) / x;
  } else {
    // Special case for Em~En. In this case, we are integrating a constant over tau\in{0,\beta}, and dividing this by
    // beta we get 1. What remains is the Boltzmann weight exp(-betaEm).
    return exp(-betaEm);
  }
}

} // namespace

#endif
