#include "nrg-general.h" // common
#include "nrg-lib.h"     // exposed in library
#include "nrg.h"         // specific to executable
#include "openmp.h"

#ifdef NRG_MPI
mpi::environment *mpienv;
mpi::communicator *mpiw;
#endif

int main(int argc, char **argv) {
#ifdef NRG_MPI
  mpi::environment env(argc, argv);
  mpi::communicator world;
  mpienv     = &env;
  mpiw       = &world;
  if (myrank() == 0) {
    std::cout << "Parallelization using MPI: Running on " << mpiw->size() << " processors." << std::endl << std::endl;
#else
    std::cout << "No MPI: single node calculation." << std::endl << std::endl;
#endif
    print_about_message();
    report_openMP();
    set_workdir(argc, argv);
    run_nrg_master();
    time_mem::memory_report();
    time_mem::timing_report();
#ifdef NRG_MPI
  } else
    run_nrg_slave();
#endif
}
