#ifndef _nrg_lib_h_
#define _nrg_lib_h_

// This header is included for both executable and library.

#define NRG_COMMON

#include <string>

void run_nrg_master();
void run_nrg_slave(); // note: only defined if compled using NRG_MPI
void set_workdir(const std::string &workdir);
void print_about_message(std::ostream &s = std::cout);

namespace time_mem {
  void timing_report();
  void memory_report();
}

#endif
