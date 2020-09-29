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
Eigen diagonalise_dsyev(Matrix &m, char jobz = 'V') {
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
  auto WORK = std::make_unique<double[]>(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK.get(), &LWORK, &INFO);
  if (INFO != 0) my_error("dsyev failed. INFO=%i", INFO);
  Eigen d(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = ham[dim * r + j];
    d.perform_checks();
  }
  return d;
}
#endif

#ifdef NRG_REAL
Eigen diagonalise_dsyevd(Matrix &m, char jobz = 'V')
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
  auto WORK = std::make_unique<double[]>(LWORK);
  auto IWORK = std::make_unique<int[]>(LIWORK);
  LAPACK_dsyevd(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK.get(), &LWORK,
                IWORK.get(), &LIWORK, &INFO);
  if (INFO != 0) {
    // dsyevd sometimes fails to converge (INFO>0). In such cases we do not trigger
    // an error but return 0, to permit error recovery.
    if (INFO > 0) return Eigen{};
    my_error("dsyevd failed. INFO=%i", INFO);
  }
  Eigen d(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = ham[dim * r + j];
    d.perform_checks();
  }
  return d;
}
#endif  

#ifdef NRG_REAL
Eigen diagonalise_dsyevr(Matrix &m, double ratio = 1.0,  char jobz = 'V')
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
  auto WORK = std::make_unique<double[]>(LWORK);
  auto IWORK  = std::make_unique<int[]>(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_dsyevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK.get(), &LWORK,
                IWORK.get(), &LIWORK, &INFO);
  if (INFO != 0) my_error("dsyevr failed. INFO=%i", INFO);
  if (MM != int(M)) {
    cout << "dsyevr computed " << MM << "/" << M << endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  Eigen d(M, dim);
  copy_values(eigenvalues, d.value, M);
  if (jobz == 'V') {
    for (size_t r = 0; r < M; r++)
      for (size_t j = 0; j < dim; j++) d.vektor(r, j) = Z[dim * r + j];
    d.perform_checks();
  }
  return d;
}
#endif

#ifdef NRG_COMPLEX
Eigen diagonalise_zheev(Matrix &m, char jobz = 'V') {
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
  auto WORK = std::make_unique<lapack_complex_double[]>(LWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheev(&jobz, &UPLO, &NN, ham, &LDA, (double *)eigenvalues, WORK.get(), &LWORK, RWORK, &INFO);
  if (INFO != 0) my_error("zheev failed. INFO=%i", INFO);
  Eigen d(dim, dim);
  copy_values(eigenvalues, d.value, dim);
  if (jobz == 'V') {
    for (size_t r = 0; r < dim; r++)
      for (size_t j = 0; j < dim; j++) {
        lapack_complex_double v = ham[dim * r + j];
        d.vektor(r, j) = cmpl(v.real, v.imag);
      }
    d.perform_checks();
  }
  return d;
}
#endif
  
#ifdef NRG_COMPLEX
Eigen diagonalise_zheevr(Matrix &m, double ratio = 1.0, char jobz = 'V') {
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
  auto WORK  = std::make_unique<lapack_complex_double[]>(LWORK);
  int LRWORK  = int(RWORK0[0]);
  auto RWORK = std::make_unique<double[]>(LRWORK);
  int LIWORK  = IWORK0[0];
  auto IWORK  = std::make_unique<int[]>(LIWORK);
  // Step 2: perform the diagonalisation
  LAPACK_zheevr(&jobz, &RANGE, &UPLO, &NN, ham, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &MM, (double *)eigenvalues, &Z[0], &LDZ, ISUPPZ, WORK.get(), &LWORK,
                RWORK.get(), &LRWORK, IWORK.get(), &LIWORK, &INFO);
  if (INFO != 0) my_error("zheevr failed. INFO=%i", INFO);
  if (MM != int(M)) {
    cout << "zheevr computed " << MM << "/" << M << endl;
    M = MM;
    my_assert(M > 0); // at least one
  }
  Eigen d(M, dim);
  copy_values(eigenvalues, d.value, M);
  if (jobz == 'V') {
    for (size_t r = 0; r < M; r++)
      for (size_t j = 0; j < dim; j++) {
        lapack_complex_double v = Z[dim * r + j];
        d.vektor(r, j) = cmpl(v.real, v.imag);
      }
    d.perform_checks();
  }
  return d;
}
#endif

void checkdiag(const Eigen &d, 
               const double NORMALIZATION_EPSILON = 1e-12,
               const double ORTHOGONALITY_EPSILON = 1e-12)
{
  const auto M = d.value.size(); // number of eigenpairs
  const auto dim = d.getrmax(); // dimension of the eigenvector
  my_assert(d.matrix0.size2() == dim);
  // Check normalization
  for (auto r = 0; r < M; r++) {
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

void dump_eigenvalues(const Eigen &d, size_t max_nr = std::numeric_limits<size_t>::max())
{
  cout << "eig= ";
  for_each_n(cbegin(d.value), min(d.value.size(), max_nr), 
           [](const t_eigen x) { cout << x << ' '; });
  cout << endl;
}

// Wrapper for the diagonalization of the Hamiltonian matrix. The number of eigenpairs returned does NOT need to be
// equal to the dimension of the matrix h. m is destroyed in the process, thus no const attribute!
Eigen diagonalise(Matrix &m) {
  time_mem::Timing t;
  check_is_matrix_upper(m);
  Eigen d;
#ifdef NRG_REAL
  if (sP.diag == "dsyev"s) 
    d = diagonalise_dsyev(m);
  if (sP.diag == "dsyevd"s) {
    d = diagonalise_dsyevd(m);
    if (d.value.size() == 0) {
      std::cout << "dsyevd failed, falling back to dsyev" << std::endl;
      d = diagonalise_dsyev(m);
    }
  }
  if (sP.diag == "dsyevr"s) 
    d = diagonalise_dsyevr(m, sP.diagratio);
#endif
#ifdef NRG_COMPLEX
  if (sP.diag == "zheev"s)  
    d = diagonalise_zheev(m);
  if (sP.diag == "zheevr"s) 
    d = diagonalise_zheevr(m, sP.diagratio);
#endif
  my_assert(d.value.size() > 0);
  my_assert(d.value.size() == d.matrix0.size1());
  my_assert(d.matrix0.size1() <= m.size1());
  my_assert(d.matrix0.size2() == m.size2());
  if (logletter('e'))
    dump_eigenvalues(d);
  if (P.checkdiag)
    checkdiag(d);
  nrglog('A', "LAPACK, dim=" << m.size1() << " M=" << d.value.size() << " [" << myrank() << "]");
  nrglog('t', "Elapsed: " << setprecision(3) << t.total_in_seconds() << " [" << myrank() << "]");
  return d;
}

#endif // _diag_h_
