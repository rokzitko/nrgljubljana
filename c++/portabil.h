// portabil.h - Portability code. Low-level stuff belongs here.
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _portabil_h_
#define _portabil_h_

using namespace std;

// Square raises x to the power of two
#ifndef HAVE_SQR
#define sqr(x) ((x) * (x))
#endif

#define BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED
#include <boost/stacktrace.hpp>

inline void print_trace() {
  std::cout << "Backtrace:" << std::endl << boost::stacktrace::stacktrace() << std::endl;
}

#define assert_isfinite(x) finite_test_fnc(x, __FILE__, __LINE__)

inline bool my_isfinite(const double x) { return std::isfinite(x); }

using cmpl = complex<double>;

inline bool my_isfinite(cmpl z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }

inline int isfinite(cmpl z) { return (std::isfinite(z.real()) && std::isfinite(z.imag()) ? 1 : 0); }

template <typename T> 
inline  T finite_test_fnc(T x, const char *file, int line) {
  if (!my_isfinite(x)) {
    std::cout << "#### EXITING DUE TO FAILED ASSERTION." << std::endl;
    std::cout << "File " << file << ", line " << line << "." << std::endl;
    std::cout << "Finiteness assertion" << std::endl;
    print_trace();
    exit(1);
  }
  return x;
}

#define my_assert(x)                                                                                                                                 \
  do {                                                                                                                                               \
    if (!(x)) {                                                                                                                                      \
      std::cout << "#### EXITING DUE TO FAILED ASSERTION." << std::endl;                                                                                       \
      std::cout << "File " << __FILE__ << ", line " << __LINE__ << "." << std::endl;                                                                           \
      std::cout << #x << std::endl;                                                                                                                            \
      print_trace();                                                                                                                                 \
      exit(1);                                                                                                                                       \
    }                                                                                                                                                \
  } while (0)

// Assert a==b
#define my_assert_equal(a, b)                                                                                                                        \
  do {                                                                                                                                               \
    if (!(a == b)) {                                                                                                                                 \
      std::cout << "#### EXITING DUE TO FAILED ASSERTION." << std::endl;                                                                                       \
      std::cout << "File " << __FILE__ << ", line " << __LINE__ << "." << std::endl;                                                                           \
      std::cout << #a << "=" << a << std::endl;                                                                                                                \
      std::cout << #b << "=" << b << std::endl;                                                                                                                \
      std::cout << "Exiting." << std::endl;                                                                                                                    \
      print_trace();                                                                                                                                 \
      exit(1);                                                                                                                                       \
    }                                                                                                                                                \
  } while (0)

#define my_assert_not_reached()                                                                                                                      \
  do {                                                                                                                                               \
    std::cout << "#### EXITING DUE TO FAILED ASSERTION." << std::endl;                                                                                         \
    std::cout << "File " << __FILE__ << ", line " << __LINE__ << "." << std::endl;                                                                             \
    std::cout << "Should never be reached." << std::endl;                                                                                                      \
    print_trace();                                                                                                                                   \
    exit(1);                                                                                                                                         \
  } while (0)

// Workaround
#ifndef __TIMESTAMP__
#define __TIMESTAMP__ ""
#endif

#ifndef ARRAYLENGTH
#define ARRAYLENGTH(x) (sizeof(x) / sizeof(*(x)))
#endif

// *** Memory usage

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach_init.h>
#include <mach/task.h>

// Source: http://blog.kuriositaet.de/?p=257
inline int getmem(unsigned int *rss, unsigned int *vs) {
  struct task_basic_info t_info{};
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) { return -1; }
  *rss = t_info.resident_size;
  *vs  = t_info.virtual_size;
  return 0;
}

inline int memoryused() {
  unsigned int rss;
  unsigned int vs;

  getmem(&rss, &vs);

  return rss / 1024; // result in kB block of memory!
}
#define HAS_MEMORY_USAGE
#endif

#if defined(__linux)
#include <string.h>
#include <stdio.h>
// Returns the number of kB blocks of memory used by this process (self). 
inline int memoryused() {
  int used = 0;

  // See /usr/src/linux/Documentation/filesystems/proc.txt
  if (FILE *info = fopen("/proc/self/status", "r")) {
    char buf[512];
    int val;

    while (fgets(buf, sizeof(buf), info)) {
//      if (!strncmp(buf, "VmSize:", 7) && sscanf(buf + 7, "%d", &val)) used = val;
      if (!strncmp(buf, "VmPeak:", 7) && sscanf(buf + 7, "%d", &val)) used = val; // peak usage is more useful information!!
    }
    fclose(info);
  }

  return used;
}
#define HAS_MEMORY_USAGE
#endif

// Fall-back
#if !defined(HAS_MEMORY_USAGE)
inline int memoryused() { return 0; }
#endif

inline ofstream safe_open(string filename, bool binary = false) {
  my_assert(filename != "");
  ios::openmode flag = ios::out;
  if (binary) flag |= ios::binary;
  ofstream F(filename, flag);
  if (!F) throw std::runtime_error(fmt::format("Can't open {} for writing", filename));
  return F;
}

inline ifstream safe_open_for_reading(string filename, bool binary = false) {
  my_assert(filename != "");
  ios::openmode flag = ios::in;
  if (binary) flag |= ios::binary;
  ifstream F(filename, flag);
  if (!F) throw std::runtime_error(fmt::format("Can't open {} for reading", filename));
  return F;
}

inline int remove(const string &filename) { return remove(filename.c_str()); } // C function

#if defined(__clang__) || defined (__GNUC__)
# define ATTRIBUTE_NO_SANITIZE_ADDRESS __attribute__((no_sanitize_address))
#else
# define ATTRIBUTE_NO_SANITIZE_ADDRESS
#endif

#if defined(__clang__) || defined (__GNUC__)
# define ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO __attribute__((no_sanitize("float-divide-by-zero")))
#else
# define ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO
#endif

#endif
