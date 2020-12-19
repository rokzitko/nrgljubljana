#ifndef _traits_hpp_
#define _traits_hpp_

#define INCL_UBLAS
#define INCL_EIGEN
#define USE_UBLAS

//#include <concepts> // C++20
#include <complex>
#include <type_traits> // is_same_v, is_floating_point_v

#ifdef INCL_UBLAS
#include <boost/numeric/ublas/matrix.hpp>
#endif

#ifdef INCL_EIGEN
#include <Eigen/Dense>
#endif

namespace NRG {

#ifdef INCL_UBLAS
using namespace boost::numeric;
#endif

template <typename T> concept floating_point = std::is_floating_point_v<T>;
template <typename T> struct is_complex : std::false_type {};
template <floating_point T> struct is_complex<std::complex<T>> : std::true_type {};
template <typename T> concept scalar = floating_point<T> || is_complex<T>::value;

#ifdef INCL_UBLAS
template <typename S> const auto generate_ublas = [](const size_t dim1, const size_t dim2) { return ublas::matrix<S>(dim1, dim2); };
template <typename S> auto size1(const ublas::matrix<S> &m) { return m.size1(); }
template <typename S> auto size2(const ublas::matrix<S> &m) { return m.size2(); }
template <typename S> auto size1(const ublas::matrix_range<const ublas::matrix<S>> &m) { return m.size1(); }
template <typename S> auto size2(const ublas::matrix_range<const ublas::matrix<S>> &m) { return m.size2(); }
#endif

#ifdef INCL_EIGEN
template <typename S> const auto generate_Eigen = [](const size_t dim1, const size_t dim2) { return Eigen::MatrixX<S>(dim1, dim2); };
template <typename S> auto size1(const Eigen::MatrixX<S> &m) { return m.rows(); }
template <typename S> auto size2(const Eigen::MatrixX<S> &m) { return m.cols(); }
#endif

template <typename T>
  concept matrix = requires(T a, T b, size_t i, size_t j) {
     { size1(a) }; // -> std::convertible_to<std::size_t>;
     { size2(a) }; // -> std::convertible_to<std::size_t>;
     { a(i,j) };
     { a = b };
     { a.swap(b) };
     typename T::value_type;
  };

// XXX: real, imag, conj for complex matrix?

template <typename T> concept real_matrix = matrix<T> && floating_point<typename T::value_type>;
template <typename T> concept complex_matrix = matrix<T> && is_complex<typename T::value_type>::value;

template <typename T> struct is_ublas_object : std::false_type {};
template <scalar S> struct is_ublas_object<ublas::matrix<S>> : std::true_type {};
template <scalar S> struct is_ublas_object<ublas::matrix_range<const ublas::matrix<S>>> : std::true_type {};

template <typename T> struct is_Eigen_object : std::false_type {};
template <scalar S> struct is_Eigen_object<Eigen::MatrixX<S>> : std::true_type {};

template <typename T> concept ublas_matrix = matrix<T> && is_ublas_object<T>::value;
template <typename T> concept Eigen_matrix = matrix<T> && is_Eigen_object<T>::value;
template <typename T> concept real_ublas_matrix = real_matrix<T> && is_ublas_object<T>::value;
template <typename T> concept real_Eigen_matrix = real_matrix<T> && is_Eigen_object<T>::value;
template <typename T> concept complex_ublas_matrix = complex_matrix<T> && is_ublas_object<T>::value;
template <typename T> concept complex_Eigen_matrix = complex_matrix<T> && is_Eigen_object<T>::value;

template <typename T>
  concept vector = requires(T a, size_t i) {
     { a.size() };
     { a[i] };
     { a.data() };
     { a.begin() };
     { a.end() };
     { a.resize(i) };
     typename T::value_type;
  };

// We encapsulate the differences between real-value and complex-value versions of the code in class traits.

template <scalar S> struct traits {};

template <> struct traits<double> {
  using t_matel = double;  // type for the matrix elements
  using t_coef = double;   // type for the Wilson chain coefficients & various prefactors
  using t_expv = double;   // type for expectation values of operators
  using t_eigen = double;  // type for the eigenvalues (always real)
  using t_temp = t_eigen;  // type for temperatures
  using t_weight = std::complex<double>;  // spectral weight accumulators (always complex)
  using evec = std::vector<double>;     // vector of eigenvalues type (always real) // YYY
  using RVector = std::vector<double>;    // vector of eigenvalues type (always real)
#ifdef USE_UBLAS
  using Matrix = ublas::matrix<t_matel>;  // matrix type
#endif
#ifdef USE_EIGEN
  using Matrix = Eigen::MatrixX<t_matel>; // matrix type
#endif
};

template <> struct traits<std::complex<double>> {
  using t_matel = std::complex<double>;
  using t_coef = std::complex<double>;
  using t_expv = std::complex<double>;     // we allow the calculation of expectation values of non-Hermitian operators!
  using t_eigen = double;
  using t_temp = t_eigen;
  using t_weight = std::complex<double>;
  using evec = std::vector<double>;
  using RVector = std::vector<double>;
#ifdef USE_UBLAS
  using Matrix = ublas::matrix<t_matel>;
#endif
#ifdef USE_EIGEN
  using Matrix = Eigen::MatrixX<t_matel>; // matrix type
#endif
};

template <scalar S> using matel_traits   = typename traits<S>::t_matel;
template <scalar S> using coef_traits    = typename traits<S>::t_coef;
template <scalar S> using expv_traits    = typename traits<S>::t_expv;
template <scalar S> using eigen_traits   = typename traits<S>::t_eigen;
template <scalar S> using weight_traits  = typename traits<S>::t_weight;
template <scalar S> using evec_traits    = typename traits<S>::evec;
template <scalar S> using RVector_traits = typename traits<S>::RVector;
template <scalar S> using Matrix_traits  = typename traits<S>::Matrix;

#ifdef USE_UBLAS
template <scalar S> const auto generate_matrix = generate_ublas<S>;
#endif
#ifdef USE_EIGEN
template <scalar S> const auto generate_matrix = generate_Eigen<S>;
#endif

template <matrix M> auto nrvec(const M &m) { return size1(m); }
template <matrix M> auto dim(const M &m) { return size2(m); }

} // namespace

#endif
