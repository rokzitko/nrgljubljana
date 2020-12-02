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
   
template<typename S, typename Matrix = Matrix_traits<S>>
[[nodiscard]] auto Zero_matrix(const size_t size1, const size_t size2) {
  return Matrix(size1, size2, 0);
}

template<typename S>
[[nodiscard]] auto Zero_matrix(const size_t size) { 
  return Zero_matrix<S>(size, size); 
}  
   
// Accumulator abstraction: automatically initialized to 0, result checked for finiteness.
template <typename T> class generic_bucket {
private:
  T value{};
public:
  generic_bucket() = default;
  // Can be constructured from a STL vector of pairs, by summing the second elements.
  template <typename T1> 
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

template<typename M> auto frobenius_norm(const ublas::matrix<M> &m) { // Frobenius norm (without taking the final square root!)
  double sum{};
  for (auto i = 0; i < m.size1(); i++)
    for (auto j = 0; j < m.size2(); j++) sum += pow(abs(m(i, j)),2);
  return sum;
}

// Check if matrix m is upper triangular. In the lower triangle, all elements must be 0. NOTE: we store the upper
// triangular part of the symmetric Hamiltonian matrix. In FORTRAN convention, this is the lower part !!
template<typename M> void check_is_matrix_upper(const ublas::matrix<M> &m) {
  my_assert(m.size1() == m.size2() && m.size1() >= 1);
  for (auto i = 1; i < m.size1(); i++)
    for (auto j = 0; j < i; j++) // j < i
      my_assert(m(i, j) == 0.);
}

// (-1)^n
CONSTFNC inline auto psgn(const int n) { return n % 2 == 0 ? 1.0 : -1.0; }

// Dump a matrix with full numerical precision. The columns are aligned for easier inspection. Expect large output!
template<typename M> inline void dump_matrix(const ublas::matrix<M> &m, std::ostream &F = std::cout, 
                                             const int header_width = 7, const int column_width = 23) {
  boost::io::ios_base_all_saver ofs(F);
  F << std::setprecision(std::numeric_limits<double>::max_digits10);
  F << fmt::format("Matrix: {}x{}\n", m.size1(), m.size2());
  for (auto r1 = 0; r1 < m.size1(); r1++) {
    F << std::setw(header_width) << r1 << ":";
    for (auto r2 = 0; r2 < m.size2(); r2++) F << std::setw(column_width) << m(r1, r2) << " ";
    F << std::endl;
  }
}

template<typename M> inline void dump_diagonal_matrix(const ublas::matrix<M> &m, const size_t max_nr, std::ostream &F = std::cout) {
  for (const auto r : range0(std::min(m.size1(), max_nr))) F << m(r,r) << ' ';
  F << std::endl;
}

template <typename T>
void save(boost::archive::binary_oarchive &oa, const ublas::matrix<T> &m) {
  oa << m.size1() << m.size2();
  for (const auto i : range0(m.size1()))
    oa << ublas::vector<T>(ublas::matrix_row<const ublas::matrix<T>>(m, i));
}

template <typename T>
void load(boost::archive::binary_iarchive &ia, ublas::matrix<T> &m) { // XXX
  const auto size1 = read_one<size_t>(ia);
  const auto size2 = read_one<size_t>(ia);
  m = ublas::matrix<T>(size1, size2);
  for (const auto i : range0(size1))
    ublas::matrix_row<ublas::matrix<T>>(m, i) = read_one<ublas::vector<T>>(ia);
}

// Chop numerical noise
template <typename T> CONSTFNC inline T chop(const T x, const double xlimit = 1.e-8) { return std::abs(x) < xlimit ? 0.0 : x; }

template<typename T>
void assert_issquare(const ublas::matrix<T> &m) { my_assert(m.size1() == m.size2()); }

// Powers, such as (-1)^n, appear in the coupling coefficients.
CONSTFNC inline double Power(const double i, const double nn) { return std::pow(i, nn); }

// Read 'size' values of type T into a ublas vector<T>.
template <typename T> ublas::vector<T> read_vector(std::istream &F, const size_t size) {
  ublas::vector<T> vec(size);
  for (auto j = 0; j < size; j++)
    vec[j] = read_one<T>(F);
  if (F.fail()) throw std::runtime_error("read_vector() error. Input file is corrupted.");
  return vec;
}

// Read values of type T into a ublas vector<T>. 'nr' is either vector dimension or the value of maximum index
template <typename T> ublas::vector<T> read_vector(std::istream &F, const bool nr_is_max_index = false) {
  const auto nr = read_one<size_t>(F);
  const auto len = nr_is_max_index ? nr+1 : nr;
  return read_vector<T>(F, len);
}

// Read 'size1' x 'size2' ublas matrix of type T.
template <typename T> ublas::matrix<T> read_matrix(std::istream &F, const size_t size1, const size_t size2) {
  ublas::matrix<T> m(size1, size2);
  for (auto j1 = 0; j1 < size1; j1++)
    for (auto j2 = 0; j2 < size2; j2++)
      m(j1, j2) = assert_isfinite( read_one<T>(F) );
  if (F.fail()) std::runtime_error("read_matrix() error. Input file is corrupted.");
  return m;
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

template <typename M> CONSTFNC auto trace_real(const ublas::matrix<M> &m) {
  assert_issquare(m);
  return ranges::accumulate(range0(m.size2()), 0.0, {}, [&m](const auto i){ return check_real(m(i, i)); });
}

inline auto csqrt(const std::complex<double> z) { return std::sqrt(z); }

template<typename T>
auto finite_size(const T &M) {
  return M.size1() && M.size2();
}
   
// M += factor * A * B^\dag
template<typename S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void product(Matrix &M, const t_coef factor, const Matrix &A, const Matrix &B) {
  if (finite_size(A) && finite_size(B)) { // if this contributes at all...
    my_assert(M.size1() == A.size1() && A.size2() == B.size2() && B.size1() == M.size2());
    my_assert(my_isfinite(factor));
    atlas::gemm(CblasNoTrans, CblasConjTrans, factor, A, B, t_coef(1.0), M);
  }
}
   
// M += factor * A * O * B^\dag
template<typename S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void rotate(Matrix &M, const t_coef factor, const Matrix &A, const Matrix &O, const Matrix &B) {
  if (finite_size(A) && finite_size(B)) {
    my_assert(M.size1() == A.size1() && A.size2() == O.size1() && O.size2() == B.size2() && B.size1() == M.size2());
    my_assert(my_isfinite(factor));
    Matrix T(O.size1(), B.size1());
    atlas::gemm(CblasNoTrans, CblasConjTrans, t_coef(1.0), O, B, t_coef(0.0), T); // T = O * B^\dag
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, A, T, t_coef(1.0), M); // M += factor * A * T
  }
}

// void rotateU(Matrix &M, const t_coef factor, const Matrix &U, const Matrix &O) {
// M += factor * U^\dag * O * U
template<typename S, typename U_type, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void rotateU(Matrix &M, const t_coef factor, const U_type &U, const Matrix &O) {
  if (finite_size(U)) {
    my_assert(M.size1() == U.size2() && U.size1() == O.size1() && O.size2() == U.size1() && U.size2() == M.size2());
    my_assert(my_isfinite(factor));
    Matrix T(U.size2(), O.size2());
    atlas::gemm(CblasConjTrans, CblasNoTrans, t_coef(1.0), U, O, t_coef(0.0), T); // T = U^\dag * O
    atlas::gemm(CblasNoTrans, CblasNoTrans, factor, T, U, t_coef(1.0), M); // M += factor * T * U
  }
}

} // namespace

#endif
