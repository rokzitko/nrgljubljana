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
#include <stdexcept>
#include <range/v3/all.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/math/special_functions/sign.hpp>

#include "traits.hpp" // defines INCL_UBLAS and/or INCL_EIGEN

#if !(defined(USE_UBLAS) || defined(USE_EIGEN))
#error "Pick one matrix backend"
#endif

#if defined(INCL_UBLAS) || defined(USE_UBLAS)
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
#endif

#if defined(INCL_EIGEN) || defined(USE_EIGEN)
#include <Eigen/Dense>
#endif

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

#ifdef INCL_UBLAS
using namespace boost::numeric;
using namespace boost::numeric::ublas; // keep this!
namespace atlas = boost::numeric::bindings::atlas;
#endif

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
 
// XXX: conj() is constexpr in C++20
[[nodiscard]] inline std::complex<double> conj_me(const std::complex<double> &z) { return conj(z); } // conjugation
[[nodiscard]] inline double conj_me(const double x) { return x; }    // no op

template<scalar S, typename Matrix = Matrix_traits<S>>
auto empty_matrix() { return Matrix(); }

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
  inline constexpr T operator+=(T x) { return value += x; }
  [[nodiscard]] inline constexpr operator T() const { return value; }
};
using bucket = generic_bucket<double>;

template <typename T>
inline constexpr bool is_odd(const T n) { return n & 1; } // must return bool
template <typename T>
inline constexpr bool is_even(const T n) { return !is_odd(n); } // must return bool

inline int my_fcmp(const double x, const double y, const double small_epsilon, const double rel_epsilon) {
  if (x == 0.0 && y == 0.0) return 0.0; // evidently equal
  if (std::abs(x) < small_epsilon && std::abs(y) < small_epsilon) return 0; // If both x and y are small, we ASSUME them to be equivalent
  if (std::abs(x-y) < rel_epsilon * (std::abs(x)+std::abs(y))) return 0;
  return boost::math::sign(x-y);
}

inline int my_fcmp(const double x, const double y, const double epsilon) { return my_fcmp(x, y, epsilon, epsilon); }

// Test if two numbers are equal to within numerical errors. (Use this for comparing values that are expected to be
// of order 1.)
inline auto num_equal(const double a, const double b, const double check_precision = 1.e-12) {
  return my_fcmp(a, b, check_precision) == 0;
}

inline auto num_equal(const std::complex<double> &a, const std::complex<double> &b, const double check_precision = 1.e-12) {
  return (my_fcmp(a.real(), b.real(), check_precision) == 0) && (my_fcmp(a.imag(), b.imag(), check_precision) == 0);
}

inline auto are_conjugate(const double a, const double b) { return num_equal(a, b); }

inline auto are_conjugate(const std::complex<double> &a, const std::complex<double> &b) { return num_equal(a.real(), b.real()) && num_equal(a.imag(), -b.imag()); }

template<matrix M> auto frobenius_norm(const M &m) { // Frobenius norm (without taking the final square root!)
  double sum{};
  for (auto i = 0; i < size1(m); i++)
    for (auto j = 0; j < size2(m); j++) sum += pow(abs(m(i, j)),2);
  return sum;
}

template<matrix M> bool is_square(const M &m) { return size1(m) == size2(m); }

// Is m upper triangular? In the lower triangle, all elements must be 0.
template<matrix M> auto is_matrix_upper(const M &m) {
  if (!is_square(m)) return false;
  for (auto i = 1; i < size1(m); i++)
    for (auto j = 0; j < i; j++) // j < i
      if (!num_equal(m(i, j), 0.0)) return false;
  return true;
}

// (-1)^n
inline constexpr auto psgn(const int n) { return n % 2 == 0 ? 1.0 : -1.0; }

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

// Chop numerical noise
template <scalar T> inline constexpr T chop(const T x, const double xlimit = 1.e-8) { return std::abs(x) < xlimit ? 0.0 : x; }

// Powers, such as (-1)^n, appear in the coupling coefficients.
inline constexpr double Power(const double i, const double nn) { return std::pow(i, nn); }

// Check if the value x is real [for complex number calculations].
constexpr inline auto is_real(const double x) { return true; }
constexpr inline auto is_real(const std::complex<double> z, const double check_real_tolerance = 1e-8) {
  return abs(z.imag()) <= check_real_tolerance;
}

// Check if x is real and return the real part, i.e. x.real().
constexpr inline auto real_part_with_check(double x) { return x; }
constexpr inline auto real_part_with_check(std::complex<double> z) {
  if (!is_real(z)) std::cout << "Warning: expected real number, but got " << z << std::endl;
  return z.real();
}

template <matrix M> auto trace_real(const M &m) {
  my_assert(is_square(m));
  return ranges::accumulate(range0(size1(m)), 0.0, {}, [&m](const auto i){ return real_part_with_check(m(i, i)); });
}

inline auto csqrt(const std::complex<double> z) { return std::sqrt(z); } // sqrt() not constexpr for complex (C++17)

template<matrix R> // 2D matrix or matrix view
auto finite_size(const R &m) { return size1(m) && size2(m); }

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

template<matrix M>
auto trace_contract(const M &A, const M &B, const size_t range) // Tr[AB]
{
  typename M::value_type sum{};
  for (const auto i : range0(range))
       for (const auto j : range0(range))
      sum += A(i, j) * B(j, i);
  return sum;
}

#if defined(INCL_UBLAS) && defined(INCL_EIGEN)

template <typename T>
Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> ublas_to_eigen(ublas::matrix<T> m){
  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor>> m_eigen(m.data().begin(),m.size1(),m.size2());
  return m_eigen;
}

template <typename T>
Eigen::Matrix<T,Eigen::Dynamic,1> ublas_to_eigen(ublas::vector<T> m){
  Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> m_eigen(m.data().begin(),m.size(),1);
  return m_eigen;
}

template<typename T, int N, int M>
auto eigen_to_ublas_matrix(Eigen::Matrix<T,N,M> m){
  ublas::matrix<T> m_ublas(m.rows(),m.cols());
  auto m1 = !m.IsRowMajor ? m.transpose() : m;
  std::copy(m1.data(), m1.data() + m1.size(),m_ublas.data().begin());
  return m_ublas;
}

template<typename T, int N, int M>
auto eigen_to_ublas_vector(Eigen::Matrix<T,N,M> m){
  ublas::vector<T> m_ublas(m.size());
  std::copy(m.data(), m.data() + m.size(),m_ublas.data().begin());
  return m_ublas;
}

#endif

#ifdef INCL_UBLAS
#include "numerics_ublas.hpp"
#endif

#ifdef INCL_EIGEN
#include "numerics_Eigen.hpp"
#endif

template<scalar S>
[[nodiscard]] auto zero_matrix(const size_t size) {
  return NRG::zero_matrix<S>(size, size);
}

template<matrix M>
[[nodiscard]] auto trim_matrix(const M &mat, const size_t new_size1, const size_t new_size2) {
  const auto old_size1 = size1(mat);
  const auto old_size2 = size2(mat);
  if (old_size1 == 0 || old_size2 == 0) return mat;                 // trimming not necessary
  my_assert(new_size1 <= old_size1 && new_size2 <= old_size2);
  if (new_size1 == old_size1 && new_size2 == old_size2) return mat; // trimming not necessary
  const auto sub = submatrix_const(mat, {0, new_size1}, {0, new_size2});
  M new_mat = sub;
  return new_mat;
} 

template<matrix M, matrix N>
  bool has_lesseq_rows(const M &A, const N &B) {
    return size1(A) <= size1(B) && size2(A) == size2(B);
}

// M = A*B, size of A is adapted to the size of B
template<matrix M>
auto prod_fit_left(const M &A, const M &B) {
  using T = typename M::value_type;
  my_assert(size1(B) <= size2(A));
  if (size1(A) == 0 || size2(B) == 0) return empty_matrix<T>();
  const auto Asub = submatrix_const(A, {0, size1(A)}, {0, size1(B)});
  return matrix_prod<T>(Asub, B);
}

// M = A*B, size of B is adapted to the size of A
template<matrix M>
auto prod_fit_right(const M &A, const M &B) {
  using T = typename M::value_type;
  my_assert(size1(B) >= size2(A));
  if (size1(A) == 0 || size2(B) == 0) return empty_matrix<T>();
  const auto Bsub = submatrix_const(B, {0, size2(A)}, {0, size2(B)});
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
  const auto Asub = submatrix_const(A, {0, size1(B)}, {0, size2(A)});
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

// Note: with dsyevr I have experienced orthogonality between eigenvectors below 1e-12. We thus use a more conservative
// epsilon for orthogonality tests of 1e-10.
// Addendum (2021): with dysevr, orthogonality can even go below 1e-10. Is it even safe to go beyond this point??
template <scalar S, typename Matrix = Matrix_traits<S>>
bool is_unitary(const Matrix &vec,
                const double NORMALIZATION_EPSILON = 1e-12,
                const double ORTHOGONALITY_EPSILON = 1e-10) {
  const auto M = nrvec(vec);
  const auto d = dim(vec);
  // Check normalization
  for (const auto r : range0(M)) {
    S sumabs{};
    for (const auto j : range0(d)) sumabs += conj_me(vec(r, j)) * vec(r, j);
    if (!num_equal(abs(sumabs), 1.0, NORMALIZATION_EPSILON)) {
      std::cout << "is_unitary() r=" << r << " : sumabs=" << sumabs << std::endl;
      return false;
    }
  }
  // Check orthogonality
  for (const auto r1 : range0(M)) {
    for (const auto r2 : boost::irange(r1 + 1, M)) {
      S skpdt{};
      for (const auto j : range0(d)) skpdt += conj_me(vec(r1, j)) * vec(r2, j);
      if (!num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON)) {
        std::cout << "is_unitary() r1=" << r1 << " r2=" << r2 << " : skpdt=" << skpdt << std::endl;
        return false;
      }
    }
  }
  return true;
}

template <scalar S, typename Matrix = Matrix_traits<S>>
bool is_unitary_blocks(const std::vector<Matrix> &U,
                       const double NORMALIZATION_EPSILON = 1e-12,
                       const double ORTHOGONALITY_EPSILON = 1e-10) {
  my_assert(U.size() > 0);
  const auto M = nrvec(U[0]);
  // Check normalization
  for (const auto r : range0(M)) {
    S sumabs{};
    for (const auto i : range0(U.size())) {
      const auto d = dim(U[i]);
      my_assert(M == nrvec(U[i]));
      for (const auto j : range0(d)) sumabs += conj_me(U[i](r, j)) * U[i](r, j);
    }
    if (!num_equal(abs(sumabs), 1.0, NORMALIZATION_EPSILON)) return false;
  }
  // Check orthogonality
  for (const auto r1 : range0(M)) {
    for (const auto r2 : boost::irange(r1 + 1, M)) {
      S skpdt{};
      for (const auto i : range0(U.size())) {
        const auto d = dim(U[i]);
        for (const auto j : range0(d)) skpdt += conj_me(U[i](r1, j)) * U[i](r2, j);
      }
      if (!num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON)) return false;
    }
  }
  return true;
}

template <scalar S, typename RVector = RVector_traits<S>, typename Matrix = Matrix_traits<S>>
void check_diag(const RVector &val, const Matrix &vec) {
  my_assert(val.size() == nrvec(vec));
  for (const auto v: val) assert_isfinite(v);
  my_assert(is_unitary<S>(vec));
}

} // namespace

#endif
