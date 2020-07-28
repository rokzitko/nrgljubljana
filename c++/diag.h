// tridiag.h - Diagonalisation code
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _diag_h_
#define _diag_h_

#define LAPACK_COMPLEX_STRUCTURE
#include "lapack.h"

bool logletter(char);

// Used by diagonalise_dsyev and diagonalise_dsyevr.
void copy_values(t_eigen *eigenvalues, EVEC &diagvalue, int M) {
  // make sure eigenvalues are in ascending order!!
  t_eigen *ptr = adjacent_find(eigenvalues, eigenvalues + M, [](t_eigen x, t_eigen y) { return x > y; });
  if (ptr != eigenvalues + M) {
    ptrdiff_t index = ptr - eigenvalues;
    cout << "WARNING: Values are not in ascending order "
         << "[index=" << index << ", M=" << M << "]: " << *ptr << " " << *(ptr + 1) << endl;
    cout << "This likely indicates a bug in LAPACK dsyev* routines." << endl;
  }
  diagvalue.resize(M);
  copy(eigenvalues, eigenvalues + M, begin(diagvalue));
}

// Perform diagonalisation: wrappers for LAPACK
// m: matrix to be diagonalised
// diag: eigenvalues and eigenvectors
// jobz: 'N' for values only, 'V' for values and vectors

#ifdef NRG_REAL
size_t diagonalise_dsyev(Matrix &m, Eigen &d, char jobz = 'V') {
  const size_t dim = m.size1();
  t_matel *ham = bindings::traits::matrix_storage(m);
  t_eigen eigenvalues[dim]; // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array
  double WORK0[1];
  // Step 1: determine optimal LWORK
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK0, &LWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK    = int(WORK0[0]);
  auto *WORK = new double[LWORK];
  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK, &LWORK, &INFO);
  if (INFO != 0) my_error("dsyev failed. INFO=%i", INFO);
  d = Eigen(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = ham[dim * r + j];
    d.perform_checks();
  }
  delete[] WORK;
  return dim;
}
#endif

#ifdef NRG_REAL
size_t diagonalise_dsyevd(Matrix &m, Eigen &d, char jobz = 'V')
{
  const size_t dim = m.size1();
  t_matel *ham = bindings::traits::matrix_storage(m);
  t_eigen eigenvalues[dim];
  char UPLO  = 'L';
  int NN     = dim;
  int LDA    = dim;
  int INFO   = 0;
  int LWORK  = -1;
  int LIWORK = -1;
  double WORK0[1];
  int IWORK0[1];
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK0, &LWORK,
                IWORK0, &LIWORK, &INFO);
  my_assert(INFO == 0);
  LWORK = int(WORK0[0]);
  LIWORK = IWORK0[0];
  auto *WORK = new double[LWORK];
  auto *IWORK = new int[LIWORK];
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK, &LWORK,
                IWORK, &LIWORK, &INFO);
  if (INFO != 0) my_error("dsyevd failed. INFO=%i", INFO);
  d = Eigen(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = ham[dim * r + j];
    d.perform_checks();
  }
  delete[] WORK;
  delete[] IWORK;
  return dim;
}
#endif  

#ifdef NRG_REAL
size_t diagonalise_dsyevr(Matrix &m, Eigen &d, char jobz = 'V',
                          double ratio = 1.0) // reduction ratio for dsyevr
{
  const size_t dim = m.size1();
  // M is the number of the eigenvalues that we will attempt to
  // calculate using dsyevr.
  size_t M = dim;
  char RANGE; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = static_cast<size_t>(ceil(ratio * M)); // round up
    M     = CLIP(M, 1ul, dim);                    // at least 1, at most dim
    RANGE = 'I';
  } else
    RANGE = 'A';
  t_matel *ham = bindings::traits::matrix_storage(m);
  t_eigen eigenvalues[dim]; // eigenvalues on exit
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
  int MM; // total number of eigenvalues found
  int LDZ = dim;
  int ISUPPZ[2 * M];
  //  The support of the eigenvectors in Z, i.e., the indices
  //  indicating the nonzero elements in Z.  The i-th eigenvector is
  //  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<t_matel> Z(LDZ * M); // eigenvectors
  int LWORK0  = -1;
  int LIWORK0 = -1;
  double WORK0[1];
  int IWORK0[1];
  // Step 1: determine optimal LWORK and LIWORK
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK0, &LWORK0,
                IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK    = int(WORK0[0]);
  int LIWORK   = IWORK0[0];
  auto *WORK = new double[LWORK];
  int *IWORK   = new int[LIWORK];
  // Step 2: perform the diagonalisation
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK, &LWORK,
                IWORK, &LIWORK, &INFO);
  if (INFO != 0) my_error("dsyevr failed. INFO=%i", INFO);
  if (MM != int(M)) {
    cout << "dsyevr computed " << MM << "/" << M << endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  d = Eigen(M, dim);
  copy_values(eigenvalues, d.value, M);
  if (jobz == 'V') {
    for (size_t r = 0; r < M; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = Z[dim * r + j];
    d.perform_checks();
  }
  delete[] WORK;
  delete[] IWORK;
  return M;
}
#endif

#ifdef NRG_COMPLEX
size_t diagonalise_zheev(Matrix &m, Eigen &d, char jobz = 'V') {
  const size_t dim = m.size1();
  lapack_complex_double *ham = (lapack_complex_double*)bindings::traits::matrix_storage(m);
  t_eigen eigenvalues[dim]; // eigenvalues on exit
  char UPLO  = 'L';         // lower triangle of a is stored
  int NN     = dim;         // the order of the matrix
  int LDA    = dim;         // the leading dimension of the array a
  int INFO   = 0;           // 0 on successful exit
  int LWORK0 = -1;          // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0[1];
  int RWORKdim = max(1ul, 3 * dim - 2);
  double RWORK[RWORKdim];
  // Step 1: determine optimal LWORK
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK0, &LWORK0, RWORK, &INFO);
  my_assert(INFO == 0);
  int LWORK  = int(WORK0[0].real);
  lapack_complex_double *WORK = new lapack_complex_double[LWORK];
  // Step 2: perform the diagonalisation
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK, &LWORK, RWORK, &INFO);
  if (INFO != 0) my_error("zheev failed. INFO=%i", INFO);
  d = Eigen(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) {
        lapack_complex_double v = ham[dim * r + j];
        d.vektor(r, j) = cmpl(v.real, v.imag);
      }
    d.perform_checks();
  }
  delete[] WORK;
  return dim;
}
#endif
  
#ifdef NRG_COMPLEX
size_t diagonalise_zheevr(Matrix &m, Eigen &d, char jobz = 'V', double ratio = 1.0) {
  const size_t dim = m.size1();
  // M is the number of the eigenvalues that we will attempt to
  // calculate using zheevr.
  size_t M = dim;
  char RANGE; // 'A'=all, 'V'=interval, 'I'=part
  if (ratio != 1.0) {
    M     = static_cast<size_t>(ceil(ratio * M)); // round up
    M     = CLIP(M, 1ul, dim);                    // at least 1, at most dim
    RANGE = 'I';
  } else
    RANGE = 'A';
  lapack_complex_double *ham = (lapack_complex_double*)bindings::traits::matrix_storage(m);
  t_eigen eigenvalues[dim]; // eigenvalues on exit
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
  int MM; // total number of eigenvalues found
  int LDZ = dim;
  int ISUPPZ[2 * M];
  //  The support of the eigenvectors in Z, i.e., the indices
  //  indicating the nonzero elements in Z.  The i-th eigenvector is
  //  nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ(2*i).
  std::vector<lapack_complex_double> Z(LDZ * M); // eigenvectors
  int LWORK0 = -1;                 // length of the WORK array (-1 == query!)
  lapack_complex_double WORK0[1];
  int LRWORK0 = -1; // query
  double RWORK0[1];
  int LIWORK0 = -1; // query
  int IWORK0[1];
  // Step 1: determine optimal LWORK, LRWORK, and LIWORK
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK0, &LWORK0,
                RWORK0, &LRWORK0, IWORK0, &LIWORK0, &INFO);
  my_assert(INFO == 0);
  int LWORK   = int(WORK0[0].real);
  lapack_complex_double *WORK  = new lapack_complex_double[LWORK];
  int LRWORK  = int(RWORK0[0]);
  double *RWORK = new double[LRWORK];
  int LIWORK  = IWORK0[0];
  int *IWORK  = new int[LIWORK];
  // Step 2: perform the diagonalisation
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK, &LWORK,
                RWORK, &LRWORK, IWORK, &LIWORK, &INFO);
  if (INFO != 0) my_error("zheevr failed. INFO=%i", INFO);
  if (MM != int(M)) {
    cout << "zheevr computed " << MM << "/" << M << endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  d = Eigen(M, dim);
  copy_values(eigenvalues, d.value, M);
  if (jobz == 'V') {
    for (size_t r = 0; r < M; r++)
      for (size_t j = 0; j < dim; j++) {
        lapack_complex_double v = Z[dim * r + j];
        d.vektor(r, j) = cmpl(v.real, v.imag);
      }
    d.perform_checks();
  }
  delete[] WORK;
  delete[] RWORK;
  delete[] IWORK;
  return M;
}
#endif

const double NORMALIZATION_EPSILON = 1e-12;
const double ORTHOGONALITY_EPSILON = 1e-12;

// Wrapper for the diagonalization of the Hamiltonian matrix in the
// NRG iteration. The number of pairs returned does NOT need to be
// equal to the dimension of the matrix h. m is destroyed in the
// process, thus no const attribute!
Eigen diagonalise(Matrix &m) {
  time_mem::Timing t;
  check_is_matrix_upper(m);
  const size_t mdim = m.size1();
  Eigen d;
  size_t M = 0; // number of computed eigenvalues
  // dr is the preferred LAPACK diagonalization routine
  dr_value dr = sP.diagroutine;
#ifdef NRG_REAL
  if (dr == diagdsyev) M = diagonalise_dsyev(m, d, 'V');
  if (dr == diagdsyevd) M = diagonalise_dsyevd(m, d, 'V');
  if (dr == diagdsyevr) M = diagonalise_dsyevr(m, d, 'V', sP.diagratio);
#endif
#ifdef NRG_COMPLEX
  if (dr == diagzheev) M = diagonalise_zheev(m, d, 'V');
  if (dr == diagzheevr) M = diagonalise_zheevr(m, d, 'V', sP.diagratio);
#endif
  my_assert(M > 0);
  my_assert(M == d.value.size());
  my_assert(M == d.matrix0.size1());
  my_assert(d.matrix0.size1() <= m.size1());
  my_assert(d.matrix0.size2() == m.size2());
  if (logletter('e')) { // Dump the first P::logenumber eigenvalues
    cout << "eig= ";
    for_each(begin(d.value), begin(d.value) + min(P::logenumber.value(), M), [](t_eigen x) { cout << x << ' '; });
    cout << endl;
  }
  if (P::checkdiag) { // default is false (since Nov 2019)
    const size_t dim = d.getrmax();
    my_assert(d.matrix0.size2() == dim);
    for (size_t r = 0; r < M; r++) {
      assert_isfinite(d.value(r));
      double sumabs = 0.0;
      for (size_t j = 0; j < dim; j++) {
        assert_isfinite(d.vektor(r, j));
        sumabs += sqr(abs(d.vektor(r, j)));
      }
      my_assert(num_equal(sumabs, 1.0, NORMALIZATION_EPSILON));
    }
    // Check orthogonality
    for (size_t r1 = 0; r1 < M; r1++)
      for (size_t r2 = r1 + 1; r2 < M; r2++) {
        t_matel skpdt = 0.0;
        for (size_t j = 0; j < dim; j++) skpdt += CONJ_ME(d.vektor(r1, j)) * d.vektor(r2, j);
        my_assert(num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON));
      }
  }
  nrglog('A', "LAPACK, dim=" << mdim << " diag=" << dr_to_string(dr)
                        << " M=" << M << " [" << myrank() << "]");
  nrglog('%', "max=" << d.value[M - 1] << " min=" << d.value[0]);
  nrglog('t', "Elapsed: " << setprecision(3) << t.total_in_seconds() << " [" << myrank() << "]");
  return d;
}

#endif // _diag_h_
