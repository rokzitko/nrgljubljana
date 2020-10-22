#ifndef _traits_h_
#define _traits_h_

// We encapsulate the differences between real-value and complex-value versions of the code in class traits.

template <typename S> struct traits {};

template <> struct traits<double> {
  using t_matel = double;  // type for the matrix elements
  using t_coef = double;   // type for the Wilson chain coefficients & various prefactors
  using t_expv = double;   // type for expectation values of operators
  using t_eigen = double;  // type for the eigenvalues (always real)
  using t_temp = t_eigen;  // type for temperatures
  using t_weight = std::complex<double>;  // spectral weight accumulators (always complex)
  using Matrix = ublas::matrix<t_matel>;  // matrix type
};

template <> struct traits<cmpl> {
  using t_matel = std::complex<double>;
  using t_coef = std::complex<double>;
  using t_expv = std::complex<double>;     // we allow the calculation of expectation values of non-Hermitian operators!
  using t_eigen = double;
  using t_temp = t_eigen;
  using t_weight = std::complex<double>;
  using Matrix = ublas::matrix<t_matel>;
};

inline cmpl conj_me(const cmpl &z) { return conj(z); } // conjugation
inline double conj_me(const double x) { return x; }    // no op

#endif