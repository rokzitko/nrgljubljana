#include <omp.h>

#ifdef MKL
#include <mkl_service.h>
#include <mkl_types.h>
#endif

// Report OpenMP parallelization settings
inline void report_openMP(std::ostream &s = std::cout) {
  s << "[OpenMP] Max. number of threads: " << omp_get_max_threads() << std::endl;
  s << "[OpenMP] Number of processors: " << omp_get_num_procs() << std::endl;
  s << "[OpenMP] Dynamic thread adjustment: " << omp_get_dynamic() << std::endl;
  s << "[OpenMP] Nested parallelism: " << omp_get_nested() << std::endl << std::endl;
#ifdef MKL
  MKLVersion version;
  mkl_get_version(&version);
  s << "Using Intel MKL library " << version.MajorVersion << "." << version.MinorVersion << "." << version.UpdateVersion << std::endl;
  s << "Processor optimization: " << version.Processor << std::endl;
  const int max_threads = mkl_get_max_threads();
  // Portability hack
# ifdef MKL_DOMAIN_BLAS
#  define NRG_MKL_DOMAIN_BLAS MKL_DOMAIN_BLAS
# else
#  define NRG_MKL_DOMAIN_BLAS MKL_BLAS
# endif
  const int blas_max_threads = mkl_domain_get_max_threads(NRG_MKL_DOMAIN_BLAS);
  const int dynamic          = mkl_get_dynamic();
  s << "max_threads=" << max_threads << " blas_max_threads=" << blas_max_threads << " dynamic=" << dynamic << std::endl << std::endl;
#endif
  s << std::endl;
}
