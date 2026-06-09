// Minimal BLAS/LAPACK ABI declarations used by this project.

#ifndef _linalg_hpp_
#define _linalg_hpp_

#include <complex>

#if defined(LAPACK_ILP64)
using lapack_int = long;
#else
using lapack_int = int;
#endif

using lapack_logical = lapack_int;
using lapack_complex_double = std::complex<double>;

#define NRG_APPEND_UNDERSCORE(name) name##_

#if defined(NRG_BLAS_SYMBOL_SUFFIX_NONE)
#define NRG_BLAS_SYMBOL(name) name
#else
#define NRG_BLAS_SYMBOL(name) NRG_APPEND_UNDERSCORE(name)
#endif

#if defined(NRG_LAPACK_SYMBOL_SUFFIX_NONE)
#define NRG_LAPACK_SYMBOL(name) name
#else
#define NRG_LAPACK_SYMBOL(name) NRG_APPEND_UNDERSCORE(name)
#endif

#if defined(__APPLE__)
#define LAPACK_dgemm NRG_BLAS_SYMBOL(dgemm)
#define LAPACK_zgemm NRG_BLAS_SYMBOL(zgemm)
#define LAPACK_dsyev NRG_LAPACK_SYMBOL(dsyev)
#define LAPACK_dsyevd NRG_LAPACK_SYMBOL(dsyevd)
#define LAPACK_dsyevr NRG_LAPACK_SYMBOL(dsyevr)
#define LAPACK_zheev NRG_LAPACK_SYMBOL(zheev)
#define LAPACK_zheevd NRG_LAPACK_SYMBOL(zheevd)
#define LAPACK_zheevr NRG_LAPACK_SYMBOL(zheevr)
#else
#define NRG_LINALG_STRINGIFY_IMPL(name) #name
#define NRG_LINALG_STRINGIFY(name) NRG_LINALG_STRINGIFY_IMPL(name)
#define NRG_BLAS_ASM_SYMBOL(name) NRG_LINALG_STRINGIFY(NRG_BLAS_SYMBOL(name))
#define NRG_LAPACK_ASM_SYMBOL(name) NRG_LINALG_STRINGIFY(NRG_LAPACK_SYMBOL(name))
#define LAPACK_dgemm NRG_LAPACK_dgemm
#define LAPACK_zgemm NRG_LAPACK_zgemm
#define LAPACK_dsyev NRG_LAPACK_dsyev
#define LAPACK_dsyevd NRG_LAPACK_dsyevd
#define LAPACK_dsyevr NRG_LAPACK_dsyevr
#define LAPACK_zheev NRG_LAPACK_zheev
#define LAPACK_zheevd NRG_LAPACK_zheevd
#define LAPACK_zheevr NRG_LAPACK_zheevr
#endif

extern "C" {

int LAPACK_dgemm(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 const double *alpha, const double *a, const lapack_int *lda, const double *b, const lapack_int *ldb,
                 const double *beta, double *c, const lapack_int *ldc)
#if !defined(__APPLE__)
  __asm__(NRG_BLAS_ASM_SYMBOL(dgemm))
#endif
  ;

int LAPACK_zgemm(const char *transa, const char *transb, const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 const double *alpha, const double *a, const lapack_int *lda, const double *b, const lapack_int *ldb,
                 const double *beta, double *c, const lapack_int *ldc)
#if !defined(__APPLE__)
  __asm__(NRG_BLAS_ASM_SYMBOL(zgemm))
#endif
  ;

void LAPACK_dsyev(const char *jobz, const char *uplo, const lapack_int *n, double *a, const lapack_int *lda, double *w,
                  double *work, const lapack_int *lwork, lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(dsyev))
#endif
  ;

void LAPACK_dsyevd(const char *jobz, const char *uplo, const lapack_int *n, double *a, const lapack_int *lda, double *w,
                   double *work, const lapack_int *lwork, lapack_int *iwork, const lapack_int *liwork, lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(dsyevd))
#endif
  ;

void LAPACK_dsyevr(const char *jobz, const char *range, const char *uplo, const lapack_int *n, double *a,
                   const lapack_int *lda, const double *vl, const double *vu, const lapack_int *il,
                   const lapack_int *iu, const double *abstol, lapack_int *m, double *w, double *z,
                   const lapack_int *ldz, lapack_int *isuppz, double *work, const lapack_int *lwork,
                   lapack_int *iwork, const lapack_int *liwork, lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(dsyevr))
#endif
  ;

void LAPACK_zheev(const char *jobz, const char *uplo, const lapack_int *n, lapack_complex_double *a,
                  const lapack_int *lda, double *w, lapack_complex_double *work, const lapack_int *lwork,
                  double *rwork, lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(zheev))
#endif
  ;

void LAPACK_zheevd(const char *jobz, const char *uplo, const lapack_int *n, lapack_complex_double *a,
                   const lapack_int *lda, double *w, lapack_complex_double *work, const lapack_int *lwork,
                   double *rwork, const lapack_int *lrwork, lapack_int *iwork, const lapack_int *liwork,
                   lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(zheevd))
#endif
  ;

void LAPACK_zheevr(const char *jobz, const char *range, const char *uplo, const lapack_int *n, lapack_complex_double *a,
                   const lapack_int *lda, const double *vl, const double *vu, const lapack_int *il,
                   const lapack_int *iu, const double *abstol, lapack_int *m, double *w, lapack_complex_double *z,
                   const lapack_int *ldz, lapack_int *isuppz, lapack_complex_double *work, const lapack_int *lwork,
                   double *rwork, const lapack_int *lrwork, lapack_int *iwork, const lapack_int *liwork,
                   lapack_int *info)
#if !defined(__APPLE__)
  __asm__(NRG_LAPACK_ASM_SYMBOL(zheevr))
#endif
  ;

}

#endif
