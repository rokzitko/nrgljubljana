// diag.h - Diagonalisation code
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _diag_hpp_
#define _diag_hpp_

#include <type_traits> // is_same_v
#include <complex>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip> // std::setprecision
#include <stdexcept>

#include "traits.hpp"
#include "params.hpp"
#include "eigen.hpp"
#include "time_mem.hpp"
#include "numerics.hpp" // is_matrix_upper
#include "debug.hpp" // nrglogdp

#define LAPACK_COMPLEX_STRUCTURE
#include "lapack.h"

namespace NRG {

// Handle complex-type conversions (used in copy_vec)
[[nodiscard]] inline auto to_matel(const double x) { return x; }
[[nodiscard]] inline auto to_matel(const lapack_complex_double &z) { return std::complex<double>(z.real, z.imag); }

template<vector SV, vector DV> requires std::is_convertible_v<typename SV::value_type, typename DV::value_type>
void copy_val(const SV &source, DV &dest, const size_t M) {
  using S = typename SV::value_type;
  my_assert(source.size() >= M);
  if (std::adjacent_find(source.begin(), source.begin() + M, std::greater<S>()) != source.begin() + M)
    std::cout << "WARNING: Values are not in ascending order. Bug in LAPACK dsyev* routines." << std::endl;
  dest.resize(M);
  std::copy_n(source.begin(), M, dest.begin());
}

template<typename U, matrix MM> // U may be _lapack_complex_double
void copy_vec(U* eigenvectors, MM & diagvectors, const size_t dim, const size_t M)
{
  diagvectors.resize(M, dim);
  for (const auto r : range0(M))
    for (const auto j : range0(dim)) diagvectors(r, j) = to_matel(eigenvectors[dim * r + j]);
}

template<scalar S, vector V, typename U>
auto copy_results(const V &eigenvalues, U* eigenvectors, const char jobz, const size_t dim, const size_t M)
{
  RawEigen<S> d(M, dim);
  copy_val(eigenvalues, d.val, M);
  if (jobz == 'V') copy_vec(eigenvectors, d.vec, dim, M);
  my_assert(d.val.size() == d.vec.size1());
  return d;
}

// Perform diagonalisation: wrappers for LAPACK. jobz: 'N' for values only, 'V' for values and vectors
template<real_matrix RM>
auto diagonalise_dsyev(RM &m, const char jobz = 'V') { // XXXXX: auto
  const auto dim = m.size1();
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array
  double WORK0 = 0;         // on exit: optimal WORK size
  // Step 1: determine optimal LWORK
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK    = int(WORK0);
  std::vector<double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevd(RM &m, const char jobz = 'V')
{
  const auto dim = m.size1();
  auto ham       = data(m);
  std::vector<double> eigenvalues(dim);
  char UPLO  = 'L';
  int NN     = dim;
  int LDA    = dim;
  int INFO   = 0;
  int LWORK  = -1;
  int LIWORK = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  int IWORK0 = 0;   // on exit: optimal IWORK size
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK, &IWORK0, &LIWORK, &INFO);
  my_assert(INFO == 0);
  LWORK      = int(WORK0);
  LIWORK     = IWORK0;
  std::vector<double> WORK(LWORK);
  std::vector<int> IWORK(LIWORK);
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) {
    // dsyevd sometimes fails to converge (INFO>0). In such cases we do not trigger
    // an error but return 0, to permit error recovery.
    if (INFO > 0)
      return RawEigen<double>();
    else
      throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  }
  return copy_results<double>(eigenvalues, ham, jobz, dim, dim);
}

template<real_matrix RM>
auto diagonalise_dsyevr(RM &m, const double ratio = 1.0, const char jobz = 'V') {
  const auto dim = m.size1();
  // M is the number of the eigenvalues that we will attempt to
  // calculate using dsyevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = static_cast<size_t>(ceil(ratio * M)); // round up
    M     = std::clamp<size_t>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = data(m);
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';     // lower triangle of a is stored
  int NN        = dim;     // the order of the matrix
  int LDA       = dim;     // the leading dimension of the array a
  int INFO      = 0;       // 0 on successful exit
  double VL     = 0;       // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  int IL        = 1; // index range
  int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  int MM{}; // total number of eigenvalues found
  int LDZ = dim;
  std::vector<int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices
  //  indicating the nonzero elements in Z.  The i-th eigenvector is
  //  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<double> Z(LDZ * M); // eigenvectors
  int LWORK0  = -1;
  int LIWORK0 = -1;
  double WORK0 = 0; // on exit: optimal WORK size
  int IWORK0 = 0;   // on exist: optimal IWORK size
  // Step 1: determine optimal LWORK and LIWORK
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK  = int(WORK0);
  int LIWORK = IWORK0;
  std::vector<double> WORK(LWORK);
  std::vector<int> IWORK(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(),
                &LWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  if (MM != int(M)) {
    std::cout << "dsyevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<double>(eigenvalues, Z.data(), jobz, dim, M);
}

template<complex_matrix CM>
auto diagonalise_zheev(CM &m, const char jobz = 'V') {
  const auto dim = m.size1();
  auto ham       = reinterpret_cast<lapack_complex_double*>(data(m));
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  int RWORKdim = std::max(1ul, 3 * dim - 2);
  std::vector<double> RWORK(RWORKdim);
  // Step 1: determine optimal LWORK
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), &WORK0, &LWORK0, RWORK.data(), &INFO);
  my_assert(INFO == 0);
  int LWORK  = int(WORK0.real);
  std::vector<lapack_complex_double> WORK(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("dsyev failed. INFO={}", INFO));
  return copy_results<std::complex<double>>(eigenvalues, ham, jobz, dim, dim);
}

template<complex_matrix CM>
auto diagonalise_zheevr(CM &m, const double ratio = 1.0, const char jobz = 'V') {
  const auto dim = m.size1();
  // M is the number of the eigenvalues that we will attempt to
  // calculate using zheevr.
  auto M = dim;
  char RANGE = 'A'; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = static_cast<size_t>(ceil(ratio * M)); // round up
    M     = std::clamp<size_t>(M, 1, dim);        // at least 1, at most dim
    RANGE = 'I';
  }
  auto ham = reinterpret_cast<lapack_complex_double*>(data(m));
  std::vector<double> eigenvalues(dim); // eigenvalues on exit
  char UPLO     = 'L';      // lower triangle of a is stored
  int NN        = dim;      // the order of the matrix
  int LDA       = dim;      // the leading dimension of the array a
  int INFO      = 0;        // 0 on successful exit
  double VL     = 0;        // value range; not referenced if RANGE != 'V'
  double VU     = 0;
  int IL        = 1; // index range
  int IU        = M;
  double ABSTOL = 0;
  // If ABSTOL=0, EPS*|T| where |T| is the 1-norm of the tridiagonal
  // matrix obtained by reducing m to tridiagonal form.
  int MM = 0; // total number of eigenvalues found
  int LDZ = dim;
  std::vector<int> ISUPPZ(2 * M);
  //  The support of the eigenvectors in Z, i.e., the indices indicating the nonzero elements in Z.  The i-th
  //  eigenvector is nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<lapack_complex_double> Z(LDZ * M); // eigenvectors
  int LWORK0 = -1;                 // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0;
  int LRWORK0 = -1;  // query
  double RWORK0 = 0; // on exit: optimal RWORK size
  int LIWORK0 = -1;  // query
  int IWORK0 = 0;    // on exit: optimal IWORK size
  // Step 1: determine optimal LWORK, LRWORK, and LIWORK
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), &WORK0, &LWORK0,
                &RWORK0, &LRWORK0, &IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK   = int(WORK0.real);
  std::vector<lapack_complex_double> WORK(LWORK);
  int LRWORK  = int(RWORK0);
  std::vector<double> RWORK(LRWORK);
  int LIWORK  = IWORK0;
  std::vector<int> IWORK(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, eigenvalues.data(), Z.data(), &LDZ, ISUPPZ.data(), WORK.data(), &LWORK,
                RWORK.data(), &LRWORK, IWORK.data(), &LIWORK, &INFO);
  if (INFO != 0) throw std::runtime_error(fmt::format("zheevr failed. INFO={}", INFO));
  if (MM != int(M)) {
    std::cout << "zheevr computed " << MM << "/" << M << std::endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  return copy_results<std::complex<double>>(eigenvalues, Z.data(), jobz, dim, M);
}

// Wrapper for the diagonalization of the Hamiltonian matrix. The number of eigenpairs returned does NOT need to be
// equal to the dimension of the matrix h. Matrix m is destroyed in the process, thus no const attribute!
template<matrix M> auto diagonalise(M &m, const DiagParams &DP, const int myrank) {
  using S = typename M::value_type;
  mpilog("diagonalise " << m.size1() << "x" << m.size2() << " " << DP.diag << " " << DP.diagratio);
  Timing timer;
  check_is_matrix_upper(m);
  RawEigen<S> d;
  if constexpr (std::is_same_v<S, double>) {
    if (DP.diag == "dsyev"s || DP.diag == "default"s) d = diagonalise_dsyev(m);
    if (DP.diag == "dsyevd"s) {
      d = diagonalise_dsyevd(m);
      if (d.getnrcomputed() == 0) {
        std::cout << "dsyevd failed, falling back to dsyev" << std::endl;
        d = diagonalise_dsyev(m);
      }
    }
    if (DP.diag == "dsyevr"s) d = diagonalise_dsyevr(m, DP.diagratio);
  }
  if constexpr (std::is_same_v<S, std::complex<double>>) {
    if (DP.diag == "zheev"s || DP.diag == "default"s) d = diagonalise_zheev(m);
    if (DP.diag == "zheevr"s) d = diagonalise_zheevr(m, DP.diagratio);
  }
  const auto nr_computed = d.getnrcomputed();
  my_assert(nr_computed > 0); // zero computed eigenvalues signals serious failure
  my_assert(has_lesseq_rows(d.vec, m)); // sanity check
  if (DP.logletter('e'))
    d.dump_eigenvalues();
  const std::string rank_string = myrank >= 0 ? " [rank=" + std::to_string(myrank) + "]" : "";
  nrglogdp('A', "LAPACK, dim=" << m.size1() << " M=" << nr_computed << rank_string);
  nrglogdp('t', "Elapsed: " << std::setprecision(3) << timer.total_in_seconds() << rank_string);
  d.check_diag();
  return d;
}

} // namespace

#endif
