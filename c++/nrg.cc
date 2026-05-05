#define NRG_EXECUTABLE

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "nrg-general.hpp" // common
#include "nrg-lib.hpp"     // exposed in library
#include "nrg.hpp"         // specific to executable
#include "openmp.hpp"      // report_openMP() called from main()
#include "workdir.hpp"

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

void configure_asan_mpi_environment()
{
#if !defined(_WIN32) && NRG_HAS_ADDRESS_SANITIZER
  // Libfabric memhooks patch allocation functions during MPI_Init and conflict with ASan.
  if (std::getenv("FI_MR_CACHE_MONITOR") == nullptr) setenv("FI_MR_CACHE_MONITOR", "disabled", 0);
#endif
}

}

inline void help(int argc, char **argv, std::string help_message)
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() >= 1 && args[0] == "-h") {
    std::cout << help_message << std::endl;
    exit(EXIT_SUCCESS);
  }
}

auto set_workdir(int argc, char **argv) { // not inline!
  std::string dir = default_workdir; // defined in workdir.h
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() == 2 && args[0] == "-w") dir = args[1];
  return std::make_unique<Workdir>(dir);
}

int main(int argc, char **argv) {
  configure_asan_mpi_environment();
  print_about_message();
  boost::mpi::environment mpienv(argc, argv);
  boost::mpi::communicator mpiw;
  if (!prepare_parallel_runtime(std::cerr, mpiw.rank() == 0)) return 1;
  if (mpiw.rank() == 0) {
    help(argc, argv, "Usage: nrg [-h] [-w workdir]");
    if (!file_exists("data")) {
      std::cout << "Input file 'data' does not exist. Terminating." << std::endl;
      return 1;
    }
    if (!file_exists("param")) {
      std::cout << "Input file 'param' does not exist. Terminating." << std::endl;
      return 1;
    }
    std::cout << "MPI job running on " << mpiw.size() << " processors." << std::endl << std::endl;
    report_openMP(std::cout, mpiw.size());
    auto workdir = set_workdir(argc, argv);
    run_nrg_master(mpienv, mpiw, std::move(workdir));
  } else {
    run_nrg_slave(mpienv, mpiw); // slaves do no disk I/O to workdir
  }
  return 0;
}
