#ifndef _openmp_hpp_
#define _openmp_hpp_

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#ifndef NRG_ENABLE_APP_OPENMP
#define NRG_ENABLE_APP_OPENMP 0
#endif

#ifndef NRG_BLAS_VENDOR
#define NRG_BLAS_VENDOR "unknown"
#endif

#ifndef NRG_HAVE_MKL
#define NRG_HAVE_MKL 0
#endif

#ifndef NRG_HAVE_OPENBLAS
#define NRG_HAVE_OPENBLAS 0
#endif

#ifndef NRG_HAVE_MKL_RT
#define NRG_HAVE_MKL_RT 0
#endif

#ifndef NRG_MKL_THREADING_LAYER
#define NRG_MKL_THREADING_LAYER "AUTO"
#endif

#ifndef NRG_MKL_REQUIRES_OPENMP
#define NRG_MKL_REQUIRES_OPENMP 0
#endif

#if NRG_ENABLE_APP_OPENMP || NRG_MKL_REQUIRES_OPENMP
#include <omp.h>
#endif

namespace NRG {

namespace detail {

inline auto lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

inline auto env_value(const char *name) -> std::string {
  if (const char *value = std::getenv(name)) return value;
  return "<unset>";
}

inline auto normalize_mkl_threading_layer(std::string layer) {
  layer = lower(std::move(layer));
  if (layer == "auto") return std::string{"AUTO"};
  if (layer == "compiler") return std::string{"COMPILER"};
  if (layer == "gnu" || layer == "gnu_thread" || layer == "gnu openmp") return std::string{"GNU"};
  if (layer == "intel" || layer == "intel_thread" || layer == "intel openmp") return std::string{"INTEL"};
  if (layer == "llvm" || layer == "llvm_thread" || layer == "llvm openmp") return std::string{"LLVM"};
  if (layer == "tbb") return std::string{"TBB"};
  if (layer == "sequential" || layer == "seq") return std::string{"SEQUENTIAL"};
  return layer;
}

inline auto configured_mkl_threading_layer() { return normalize_mkl_threading_layer(NRG_MKL_THREADING_LAYER); }

inline auto environment_mkl_threading_layer() {
  if (const char *value = std::getenv("MKL_THREADING_LAYER")) return normalize_mkl_threading_layer(value);
  return std::string{};
}

inline auto active_mkl_threading_layer() {
  const auto env_layer = environment_mkl_threading_layer();
  if (!env_layer.empty()) return env_layer;
  const auto configured_layer = configured_mkl_threading_layer();
  return configured_layer == "AUTO" ? std::string{} : configured_layer;
}

inline auto mkl_layer_openmp_family(const std::string &layer) -> std::string {
  if (layer == "GNU") return "GNU";
  if (layer == "INTEL") return "Intel";
  if (layer == "LLVM") return "LLVM";
  return "";
}

inline auto mkl_layer_description(const std::string &layer) -> std::string {
  if (layer == "GNU") return "GNU OpenMP";
  if (layer == "INTEL") return "Intel OpenMP";
  if (layer == "LLVM") return "LLVM OpenMP";
  if (layer == "TBB") return "Intel TBB";
  if (layer == "SEQUENTIAL") return "sequential";
  if (layer == "AUTO") return "automatic";
  return layer.empty() ? "<unset>" : layer;
}

inline void set_mkl_threading_layer_from_build(std::ostream *s, const bool verbose) {
#if NRG_HAVE_MKL
  const auto configured_layer = configured_mkl_threading_layer();
  if (configured_layer == "AUTO" || configured_layer == "COMPILER" || std::getenv("MKL_THREADING_LAYER") != nullptr) return;
#if defined(_WIN32)
  _putenv_s("MKL_THREADING_LAYER", configured_layer.c_str());
#else
  setenv("MKL_THREADING_LAYER", configured_layer.c_str(), 0);
#endif
  if (verbose && s != nullptr) *s << "[MKL] MKL_THREADING_LAYER was unset; using build-time layer " << configured_layer << std::endl;
#else
  (void)s;
  (void)verbose;
#endif
}

inline auto symbol(const char *name) -> void * { return dlsym(RTLD_DEFAULT, name); }

template <typename Function>
inline auto symbol_as(const char *name) -> Function { return reinterpret_cast<Function>(symbol(name)); }

inline auto symbol_library(void *sym) -> std::string {
  Dl_info info{};
  if (sym != nullptr && dladdr(sym, &info) != 0 && info.dli_fname != nullptr) return info.dli_fname;
  return "<unknown>";
}

inline auto openmp_runtime_family(const std::string &path) -> std::string {
  const auto lower_path = lower(path);
  if (lower_path.find("libiomp5") != std::string::npos || lower_path.find("libguide") != std::string::npos) return "Intel";
  if (lower_path.find("libgomp") != std::string::npos) return "GNU";
  if (lower_path.find("libomp.so") != std::string::npos || lower_path.find("libomp.dylib") != std::string::npos) return "LLVM";
  return "";
}

inline auto loaded_openmp_runtimes() {
  std::ifstream maps("/proc/self/maps");
  std::set<std::string> paths;
  std::string line;
  while (std::getline(maps, line)) {
    const auto first_slash = line.find('/');
    if (first_slash == std::string::npos) continue;
    auto path = line.substr(first_slash);
    if (!openmp_runtime_family(path).empty()) paths.insert(std::move(path));
  }
  return std::vector<std::string>(paths.begin(), paths.end());
}

inline auto online_cpu_count() -> long {
  const auto cpus = sysconf(_SC_NPROCESSORS_ONLN);
  return cpus > 0 ? cpus : 0;
}

inline auto mkl_threading_layer_name(const int layer) -> const char * {
  switch (layer) {
  case 0: return "Intel OpenMP";
  case 1: return "sequential";
  case 2: return "PGI OpenMP";
  case 3: return "GNU OpenMP";
  case 4: return "TBB";
  default: return "unknown";
  }
}

inline auto openblas_parallel_name(const int model) -> const char * {
  switch (model) {
  case 0: return "sequential";
  case 1: return "pthreads";
  case 2: return "OpenMP";
  default: return "unknown";
  }
}

inline void reference_linked_openmp_runtime() {
#if NRG_MKL_REQUIRES_OPENMP
  (void)omp_get_num_procs();
#endif
}

inline auto validate_mkl_threading_runtime(std::ostream *s, const bool verbose) -> bool {
#if NRG_HAVE_MKL
  reference_linked_openmp_runtime();
  const auto layer = active_mkl_threading_layer();
  const auto expected_family = mkl_layer_openmp_family(layer);
  if (expected_family.empty()) return true;

#if NRG_HAVE_MKL_RT
  const auto configured_layer = configured_mkl_threading_layer();
  if (configured_layer == "AUTO" || configured_layer == "COMPILER" || !NRG_MKL_REQUIRES_OPENMP) {
    if (verbose && s != nullptr) {
      *s << "ERROR: MKL_THREADING_LAYER=" << layer << " selects " << mkl_layer_description(layer)
         << ", but this mkl_rt build did not link that OpenMP runtime." << std::endl;
      *s << "ERROR: Reconfigure with -DNRGLJUBLJANA_MKL_THREADING_LAYER=" << layer
         << " so CMake links the matching OpenMP runtime through OpenMP::OpenMP_CXX." << std::endl;
    }
    return false;
  }
  if (layer != configured_layer) {
    if (verbose && s != nullptr) {
      *s << "ERROR: MKL_THREADING_LAYER=" << layer << " conflicts with build-time NRGLJUBLJANA_MKL_THREADING_LAYER=" << configured_layer << "." << std::endl;
      *s << "ERROR: Use MKL_THREADING_LAYER=" << configured_layer << " at runtime, or reconfigure with -DNRGLJUBLJANA_MKL_THREADING_LAYER=" << layer << "." << std::endl;
    }
    return false;
  }
#endif

  void *num_procs_symbol = symbol("omp_get_num_procs");
  void *max_threads_symbol = symbol("omp_get_max_threads");
  if (num_procs_symbol == nullptr || max_threads_symbol == nullptr) {
    if (verbose && s != nullptr) {
      *s << "ERROR: MKL_THREADING_LAYER=" << layer << " requires " << mkl_layer_description(layer)
         << " runtime symbols, but omp_get_num_procs/omp_get_max_threads are not visible in this process." << std::endl;
      *s << "ERROR: Reconfigure with -DNRGLJUBLJANA_MKL_THREADING_LAYER=" << layer
         << " so the executable links the matching OpenMP runtime, or choose an MKL_THREADING_LAYER compatible with the linked runtime." << std::endl;
    }
    return false;
  }
#else
  (void)s;
  (void)verbose;
#endif
  return true;
}

inline auto mkl_thread_queries_are_safe(std::ostream *s, const bool verbose) -> bool {
#if NRG_HAVE_MKL
  if (!validate_mkl_threading_runtime(s, verbose)) return false;
#if NRG_HAVE_MKL_RT
  const auto layer = active_mkl_threading_layer();
  if (layer.empty()) {
    if (verbose && s != nullptr)
      *s << "WARNING: Intel MKL is linked through mkl_rt without an explicit MKL_THREADING_LAYER. Skipping MKL thread-count queries because they can initialize an unvalidated threading backend." << std::endl;
    return false;
  }
  if (layer == "TBB") {
    if (verbose && s != nullptr)
      *s << "WARNING: MKL_THREADING_LAYER=TBB is not validated by NRG Ljubljana startup checks. Skipping MKL thread-count queries." << std::endl;
    return false;
  }
#endif
#else
  (void)s;
  (void)verbose;
#endif
  return true;
}

inline auto blas_thread_count() -> int {
  if (symbol("MKL_Get_Max_Threads") != nullptr || symbol("mkl_get_max_threads") != nullptr) {
    if (!mkl_thread_queries_are_safe(nullptr, false)) return 0;
    if (const auto mkl_get_max_threads = symbol_as<int (*)()>("MKL_Get_Max_Threads")) return mkl_get_max_threads();
    if (const auto mkl_get_max_threads = symbol_as<int (*)()>("mkl_get_max_threads")) return mkl_get_max_threads();
  }
  if (const auto openblas_get_num_threads = symbol_as<int (*)()>("openblas_get_num_threads")) return openblas_get_num_threads();
  return 0;
}

inline void report_environment(std::ostream &s) {
  constexpr const char *vars[] = {
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OMP_MAX_ACTIVE_LEVELS",
    "OMP_DYNAMIC",
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "MKL_NUM_THREADS",
    "MKL_DOMAIN_NUM_THREADS",
    "MKL_DYNAMIC",
    "MKL_THREADING_LAYER",
    "OPENBLAS_NUM_THREADS",
    "GOTO_NUM_THREADS",
    "BLIS_NUM_THREADS"
  };
  s << "[Parallel] Relevant environment:" << std::endl;
  for (const auto *var : vars) s << "[Parallel]   " << var << '=' << env_value(var) << std::endl;
}

inline void report_loaded_openmp_runtimes(std::ostream &s) {
  const auto runtimes = loaded_openmp_runtimes();
  if (runtimes.empty()) {
    s << "[OpenMP runtime] No OpenMP runtime library currently mapped" << std::endl;
    return;
  }
  s << "[OpenMP runtime] Loaded runtime libraries:" << std::endl;
  std::set<std::string> families;
  for (const auto &runtime : runtimes) {
    const auto family = openmp_runtime_family(runtime);
    families.insert(family);
    s << "[OpenMP runtime]   " << family << ": " << runtime << std::endl;
  }
  if (families.size() > 1) s << "WARNING: multiple OpenMP runtime families are loaded; this is unsafe and can crash threaded BLAS/LAPACK calls." << std::endl;
}

inline void report_openmp_symbols(std::ostream &s) {
  void *max_threads_symbol = symbol("omp_get_max_threads");
  if (max_threads_symbol == nullptr) return;
  const auto omp_get_max_threads = reinterpret_cast<int (*)()>(max_threads_symbol);
  s << "[OpenMP runtime] omp_get_max_threads provider: " << symbol_library(max_threads_symbol) << std::endl;
  s << "[OpenMP runtime] omp_get_max_threads: " << omp_get_max_threads() << std::endl;
  if (const auto omp_get_dynamic = symbol_as<int (*)()>("omp_get_dynamic"))
    s << "[OpenMP runtime] omp_get_dynamic: " << omp_get_dynamic() << std::endl;
  if (const auto omp_get_num_procs = symbol_as<int (*)()>("omp_get_num_procs"))
    s << "[OpenMP runtime] omp_get_num_procs: " << omp_get_num_procs() << std::endl;
}

inline void report_mkl(std::ostream &s) {
  void *version_symbol = symbol("MKL_Get_Version_String");
  void *max_threads_symbol = symbol("MKL_Get_Max_Threads");
  if (max_threads_symbol == nullptr) max_threads_symbol = symbol("mkl_get_max_threads");
  if (version_symbol == nullptr && max_threads_symbol == nullptr) {
#if NRG_HAVE_MKL
    s << "[MKL] MKL was detected at build time, but MKL service symbols are not visible for runtime reporting." << std::endl;
#endif
    return;
  }

  s << "[MKL] Runtime library: " << symbol_library(version_symbol != nullptr ? version_symbol : max_threads_symbol) << std::endl;
  if (version_symbol != nullptr) {
    const auto mkl_get_version_string = reinterpret_cast<void (*)(char *, int)>(version_symbol);
    std::array<char, 512> version{};
    mkl_get_version_string(version.data(), static_cast<int>(version.size() - 1));
    version.back() = '\0';
    s << "[MKL] Version: " << version.data() << std::endl;
  }
  if (!mkl_thread_queries_are_safe(&s, true)) {
    s << "[MKL] Threading diagnostics requiring MKL initialization were skipped." << std::endl;
    return;
  }
  if (max_threads_symbol != nullptr) {
    const auto mkl_get_max_threads = reinterpret_cast<int (*)()>(max_threads_symbol);
    s << "[MKL] Max threads: " << mkl_get_max_threads() << std::endl;
  }
  if (const auto mkl_domain_get_max_threads = symbol_as<int (*)(int)>("MKL_Domain_Get_Max_Threads")) {
    constexpr int mkl_domain_blas = 1;
    s << "[MKL] BLAS domain max threads: " << mkl_domain_get_max_threads(mkl_domain_blas) << std::endl;
  }
  if (const auto mkl_get_dynamic = symbol_as<int (*)()>("MKL_Get_Dynamic"))
    s << "[MKL] Dynamic adjustment: " << mkl_get_dynamic() << std::endl;
  if (const auto mkl_get_threading_layer = symbol_as<int (*)()>("mkl_get_threading_layer")) {
    const int layer = mkl_get_threading_layer();
    s << "[MKL] Threading layer: " << mkl_threading_layer_name(layer) << " (" << layer << ")" << std::endl;
    if (layer == 1) s << "WARNING: MKL reports the sequential threading layer; BLAS/LAPACK calls will not use numerical multithreading." << std::endl;
  }
}

inline void report_openblas(std::ostream &s) {
  void *config_symbol = symbol("openblas_get_config");
  void *threads_symbol = symbol("openblas_get_num_threads");
  if (config_symbol == nullptr && threads_symbol == nullptr) {
#if NRG_HAVE_OPENBLAS
    s << "[OpenBLAS] OpenBLAS was detected at build time, but OpenBLAS reporting symbols are not visible." << std::endl;
#endif
    return;
  }

  s << "[OpenBLAS] Runtime library: " << symbol_library(config_symbol != nullptr ? config_symbol : threads_symbol) << std::endl;
  if (config_symbol != nullptr) {
    const auto openblas_get_config = reinterpret_cast<const char *(*)()>(config_symbol);
    s << "[OpenBLAS] Config: " << openblas_get_config() << std::endl;
  }
  if (const auto openblas_get_corename = symbol_as<const char *(*)()>("openblas_get_corename"))
    s << "[OpenBLAS] Core: " << openblas_get_corename() << std::endl;
  if (threads_symbol != nullptr) {
    const auto openblas_get_num_threads = reinterpret_cast<int (*)()>(threads_symbol);
    s << "[OpenBLAS] Threads: " << openblas_get_num_threads() << std::endl;
  }
  if (const auto openblas_get_parallel = symbol_as<int (*)()>("openblas_get_parallel")) {
    const int model = openblas_get_parallel();
    s << "[OpenBLAS] Threading model: " << openblas_parallel_name(model) << " (" << model << ")" << std::endl;
    if (model == 0) s << "WARNING: OpenBLAS reports the sequential threading model; BLAS/LAPACK calls will not use numerical multithreading." << std::endl;
  }
}

inline void report_application_openmp(std::ostream &s) {
#if NRG_ENABLE_APP_OPENMP
  s << "[Application OpenMP] Enabled at build time" << std::endl;
  s << "[Application OpenMP] _OPENMP: " << _OPENMP << std::endl;
  s << "[Application OpenMP] Max threads: " << omp_get_max_threads() << std::endl;
  s << "[Application OpenMP] Number of processors: " << omp_get_num_procs() << std::endl;
  s << "[Application OpenMP] Dynamic thread adjustment: " << omp_get_dynamic() << std::endl;
# if _OPENMP >= 200805
  s << "[Application OpenMP] Max active levels: " << omp_get_max_active_levels() << std::endl;
# endif
#else
  s << "[Application OpenMP] Disabled at build time; BLAS/LAPACK owns numerical threading." << std::endl;
#endif
}

inline void report_parallel_warnings(std::ostream &s, const int mpi_size) {
  const int blas_threads = blas_thread_count();
  const long cpus = online_cpu_count();
  if (blas_threads == 1) s << "WARNING: BLAS/LAPACK currently reports one thread; numerical kernels will run serially unless this is changed with the BLAS threading environment." << std::endl;
  if (mpi_size > 1 && blas_threads > 0) {
    const long total_threads = static_cast<long>(mpi_size) * blas_threads;
    s << "[Parallel] MPI ranks x BLAS threads: " << mpi_size << " x " << blas_threads << " = " << total_threads << std::endl;
    if (cpus > 0 && total_threads > cpus)
      s << "WARNING: MPI ranks times BLAS threads exceeds online CPUs (" << cpus << "); expect oversubscription unless CPU binding/cgroups intentionally limit this job." << std::endl;
  }
#if NRG_ENABLE_APP_OPENMP
  if (blas_threads > 1 && omp_get_max_threads() > 1)
    s << "WARNING: application OpenMP and BLAS/LAPACK threading are both enabled. Use this only intentionally; nested parallelism can oversubscribe and may expose OpenMP runtime incompatibilities." << std::endl;
#endif
}

} // namespace detail

inline auto prepare_parallel_runtime(std::ostream &s = std::cerr, const bool verbose = true) -> bool {
  detail::set_mkl_threading_layer_from_build(&s, verbose);
  return detail::validate_mkl_threading_runtime(&s, verbose);
}

// Report parallelization settings and the numerical library threading backend.
inline void report_openMP(std::ostream &s = std::cout, const int mpi_size = 1) {
  const bool runtime_ok = prepare_parallel_runtime(s, true);
  s << "[Parallel] Build-time BLAS/LAPACK vendor: " << NRG_BLAS_VENDOR << std::endl;
  s << "[Parallel] Online CPUs: " << detail::online_cpu_count() << std::endl;
  detail::report_environment(s);
  detail::report_application_openmp(s);
  if (!runtime_ok) {
    s << "[Parallel] Startup runtime validation failed; skipping BLAS/LAPACK runtime queries that could trigger a loader error." << std::endl << std::endl;
    return;
  }
  detail::report_mkl(s);
  detail::report_openblas(s);
  detail::report_loaded_openmp_runtimes(s);
  detail::report_openmp_symbols(s);
  detail::report_parallel_warnings(s, mpi_size);
  s << std::endl;
}

inline void warn_application_openmp_request(const std::string &diag_mode, const int diagth, std::ostream &s = std::cout) {
#if NRG_ENABLE_APP_OPENMP
  if (diag_mode == "OpenMP" && diagth > 1 && detail::blas_thread_count() > 1)
    s << "WARNING: diag_mode=OpenMP with diagth=" << diagth << " runs simultaneous diagonalisations while BLAS/LAPACK is threaded; this is nested parallelism and should be used only intentionally." << std::endl;
#else
  if (diag_mode == "OpenMP")
    s << "WARNING: diag_mode=OpenMP was requested, but application OpenMP is disabled at build time; diagonalisation scheduling is serial and BLAS/LAPACK remains threaded." << std::endl;
  if (diagth > 1)
    s << "WARNING: diagth=" << diagth << " was requested, but application OpenMP is disabled at build time; diagth is ignored." << std::endl;
#endif
}

} // namespace

#endif
