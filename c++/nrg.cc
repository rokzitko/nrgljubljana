#include "nrg-general.h" // common
#include "nrg-lib.h"     // exposed in library
#include "nrg.h"         // specific to executable

int main(int argc, char *argv[]) {
#ifdef NRG_MPI
  mpi::environment env(argc, argv);
  mpi::communicator world;
  mpienv     = &env;
  mpiw       = &world;
  int myrank = mpiw->rank();
  if (myrank == 0) {
    set_workdir(argc, argv);
    std::cout << "Parallelization using MPI: Running on " << mpiw->size() << " processors." << std::endl << std::endl;
    run_nrg_master();
    memory_report();
    timing_report();
  } else
    run_nrg_slave();
#else
  set_workdir(argc, argv);
  std::cout << "No MPI: single node calculation." << std::endl << std::endl;
  run_nrg_master();
  memory_report();
  timing_report();
#endif // NRG_MPI
}
