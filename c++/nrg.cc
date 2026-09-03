#if defined(NRGLJUBLJANA_ENABLE_FP_TRAPS) && NRGLJUBLJANA_ENABLE_FP_TRAPS
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif

#define NRG_EXECUTABLE

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <boost/mpi/collectives.hpp>
#if defined(NRGLJUBLJANA_ENABLE_FP_TRAPS) && NRGLJUBLJANA_ENABLE_FP_TRAPS
#include <cstdio>
#include <fenv.h>
#endif

#include "nrg-general.hpp" // common
#include "nrg-lib.hpp"     // exposed in library
#include "nrg.hpp"         // specific to executable
#include "openmp.hpp"      // report_openMP() called from main()
#include "workdir.hpp"
#include <common/version.hpp>

#ifndef NRG_BUILD_CXX_COMPILER_ID
 #define NRG_BUILD_CXX_COMPILER_ID ""
#endif
#ifndef NRG_BUILD_CXX_COMPILER_VERSION
 #define NRG_BUILD_CXX_COMPILER_VERSION ""
#endif
#ifndef NRG_BUILD_SYSTEM_PROCESSOR
 #define NRG_BUILD_SYSTEM_PROCESSOR ""
#endif
#ifndef NRG_BUILD_SYSTEM_NAME
 #define NRG_BUILD_SYSTEM_NAME ""
#endif
#ifndef NRG_BUILD_SYSTEM_VERSION
 #define NRG_BUILD_SYSTEM_VERSION ""
#endif

#define NRG_STRINGIFY_IMPL(x) #x
#define NRG_STRINGIFY(x) NRG_STRINGIFY_IMPL(x)

#if defined(NRGLJUBLJANA_ENABLE_FP_TRAPS) && NRGLJUBLJANA_ENABLE_FP_TRAPS
#pragma STDC FENV_ACCESS ON
#endif

using namespace NRG;

#ifndef NRG_HAS_ADDRESS_SANITIZER
#if defined(__SANITIZE_ADDRESS__)
#define NRG_HAS_ADDRESS_SANITIZER 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define NRG_HAS_ADDRESS_SANITIZER 1
#endif
#endif
#endif

#ifndef NRG_HAS_ADDRESS_SANITIZER
#define NRG_HAS_ADDRESS_SANITIZER 0
#endif

namespace {

auto has_text(const char *value) -> bool { return value != nullptr && value[0] != '\0'; }

auto fallback_compiler_name() -> const char * {
#if defined(__INTEL_LLVM_COMPILER)
  return "IntelLLVM";
#elif defined(__INTEL_COMPILER)
  return "Intel";
#elif defined(__apple_build_version__) && defined(__clang__)
  return "AppleClang";
#elif defined(__clang__)
  return "Clang";
#elif defined(__GNUC__)
  return "GCC";
#elif defined(_MSC_VER)
  return "MSVC";
#else
  return "unknown";
#endif
}

auto fallback_compiler_version() -> const char * {
#if defined(__INTEL_LLVM_COMPILER)
  return NRG_STRINGIFY(__INTEL_LLVM_COMPILER);
#elif defined(__INTEL_COMPILER)
  return NRG_STRINGIFY(__INTEL_COMPILER);
#elif defined(__clang_version__)
  return __clang_version__;
#elif defined(__VERSION__)
  return __VERSION__;
#elif defined(_MSC_FULL_VER)
  return NRG_STRINGIFY(_MSC_FULL_VER);
#elif defined(_MSC_VER)
  return NRG_STRINGIFY(_MSC_VER);
#else
  return "";
#endif
}

auto fallback_architecture() -> const char * {
#if defined(__x86_64__) || defined(_M_X64)
  return "x86_64";
#elif defined(__i386__) || defined(_M_IX86)
  return "x86";
#elif defined(__aarch64__) || defined(_M_ARM64)
  return "aarch64";
#elif defined(__arm__) || defined(_M_ARM)
  return "arm";
#elif defined(__powerpc64__) || defined(__ppc64__)
  return "ppc64";
#elif defined(__powerpc__) || defined(__ppc__) || defined(_M_PPC)
  return "ppc";
#elif defined(__riscv)
 #if defined(__riscv_xlen) && __riscv_xlen == 64
  return "riscv64";
 #elif defined(__riscv_xlen) && __riscv_xlen == 32
  return "riscv32";
 #else
  return "riscv";
 #endif
#else
  return "unknown";
#endif
}

auto fallback_os_name() -> const char * {
#if defined(__linux__)
  return "Linux";
#elif defined(__APPLE__) && defined(__MACH__)
  return "Darwin";
#elif defined(_WIN32)
  return "Windows";
#elif defined(__FreeBSD__)
  return "FreeBSD";
#elif defined(__OpenBSD__)
  return "OpenBSD";
#elif defined(__NetBSD__)
  return "NetBSD";
#elif defined(__unix__)
  return "Unix";
#else
  return "unknown";
#endif
}

auto compiler_info() {
  std::string compiler = has_text(NRG_BUILD_CXX_COMPILER_ID) ? NRG_BUILD_CXX_COMPILER_ID : fallback_compiler_name();
  if (compiler == "GNU") compiler = "GCC";
  const char *version = has_text(NRG_BUILD_CXX_COMPILER_VERSION) ? NRG_BUILD_CXX_COMPILER_VERSION : fallback_compiler_version();
  if (has_text(version)) {
    compiler += ' ';
    compiler += version;
  }
  return compiler;
}

auto architecture_info() {
  return has_text(NRG_BUILD_SYSTEM_PROCESSOR) ? std::string{NRG_BUILD_SYSTEM_PROCESSOR} : std::string{fallback_architecture()};
}

auto os_info() {
  std::string os = has_text(NRG_BUILD_SYSTEM_NAME) ? NRG_BUILD_SYSTEM_NAME : fallback_os_name();
  if (has_text(NRG_BUILD_SYSTEM_VERSION)) {
    os += ' ';
    os += NRG_BUILD_SYSTEM_VERSION;
  }
  return os;
}

void print_nrg_about_message() {
  fmt::print("NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print("Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print("Compiled on {} at {}\n", __DATE__, __TIME__);
  fmt::print("Compiler/platform: {}; arch={}; os={}\n\n", compiler_info(), architecture_info(), os_info());
}

void configure_asan_mpi_environment()
{
#if !defined(_WIN32) && NRG_HAS_ADDRESS_SANITIZER
  // Libfabric memhooks patch allocation functions during MPI_Init and conflict with ASan.
  if (std::getenv("FI_MR_CACHE_MONITOR") == nullptr) setenv("FI_MR_CACHE_MONITOR", "disabled", 0);
#endif
}

void enable_fp_traps()
{
#if defined(NRGLJUBLJANA_ENABLE_FP_TRAPS) && NRGLJUBLJANA_ENABLE_FP_TRAPS
  feclearexcept(FE_ALL_EXCEPT);
  if (feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) != -1) {
    std::puts("Floating point exception trapping enabled");
  }
#endif
}

}

inline auto help(int argc, char **argv, const std::string &help_message) -> bool
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (!args.empty() && args[0] == "-h") {
    std::cout << help_message << std::endl;
    return true;
  }
  return false;
}

struct WorkdirOptions {
  std::string dir;
  WorkdirMode mode;
};

auto parse_workdir_options(int argc, char **argv) {
  std::string dir = default_workdir; // defined in workdir.h
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  auto mode = WorkdirMode::unique_temporary;
  bool workdir_option_seen = false;

  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index]; // NOLINT
    if (option != "-w" && option != "--checkpoint-dir") {
      throw std::invalid_argument("Unknown argument '" + option + "'");
    }
    if (workdir_option_seen) {
      throw std::invalid_argument("Specify only one of -w DIR and --checkpoint-dir DIR");
    }
    if (index + 1 >= argc) throw std::invalid_argument("Option " + option + " requires DIR");

    const std::string value = argv[++index]; // NOLINT
    if (value.empty() || value == "-h" || value == "-w" || value == "--checkpoint-dir") {
      throw std::invalid_argument("Option " + option + " requires DIR");
    }
    dir = value;
    mode = option == "--checkpoint-dir" ? WorkdirMode::persistent_exact : WorkdirMode::unique_temporary;
    workdir_option_seen = true;
  }
  return WorkdirOptions{dir, mode};
}

int main(int argc, char **argv) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "nrg")) return EXIT_SUCCESS;
  configure_asan_mpi_environment();
  boost::mpi::environment mpienv(argc, argv);
  boost::mpi::communicator mpiw;
  if (mpiw.rank() == 0) print_nrg_about_message();
  enable_fp_traps();
  if (!prepare_parallel_runtime(std::cerr, mpiw.rank() == 0)) return 1;

  constexpr int startup_continue = -1;
  const std::string help_message =
    "Usage: nrg [-h] [-V|--version] [-w DIR | --checkpoint-dir DIR]\n"
    "  -w DIR                 create a unique temporary workdir below DIR\n"
    "  --checkpoint-dir DIR   create or reopen exactly DIR and preserve it";
  int startup_status = startup_continue;
  std::unique_ptr<Workdir> workdir;
  if (mpiw.rank() == 0) {
    if (help(argc, argv, help_message)) {
      startup_status = EXIT_SUCCESS;
    } else {
      try {
        const auto options = parse_workdir_options(argc, argv);
        workdir = std::make_unique<Workdir>(options.dir, options.mode, true);
        clear_done_marker();
        if (!file_exists("data")) {
          std::cout << "Input file 'data' does not exist. Terminating." << std::endl;
          startup_status = EXIT_FAILURE;
        } else if (!file_exists("param")) {
          std::cout << "Input file 'param' does not exist. Terminating." << std::endl;
          startup_status = EXIT_FAILURE;
        }
      } catch (const std::invalid_argument &error) {
        std::cerr << "nrg: " << error.what() << '\n' << help_message << std::endl;
        startup_status = EXIT_FAILURE;
      } catch (const std::runtime_error &error) {
        std::cerr << "nrg: " << error.what() << std::endl;
        startup_status = EXIT_FAILURE;
      }
    }
  }

  boost::mpi::broadcast(mpiw, startup_status, 0);
  if (startup_status != startup_continue) return startup_status;

  if (mpiw.rank() == 0) {
    std::cout << "MPI job running on " << mpiw.size() << " processors." << std::endl << std::endl;
    report_openMP(std::cout, mpiw.size());
    std::cout << "workdir=" << workdir->get() << std::endl << std::endl;
    try {
      run_nrg_master(mpienv, mpiw, std::move(workdir));
    } catch (const std::exception &error) {
      std::cerr << "nrg: " << error.what() << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    run_nrg_slave(mpienv, mpiw); // slaves do no disk I/O to workdir
  }
  return 0;
}
