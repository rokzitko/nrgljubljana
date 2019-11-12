#ifndef _nrg_lib_h_
#define _nrg_lib_h_

#include <string>

void run_nrg_master();
void run_nrg_slave(); // note: only defined if compled using NRG_MPI
void set_workdir(std::string workdir);

#endif
