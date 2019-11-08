#include "nrg-general.h"
#include "nrg-lib.h"
#include "nrg.h"

int main(int argc, char *argv[])
{
#ifdef NRG_MPI
   mpi::environment env(argc, argv);
   mpi::communicator world;
   mpienv = &env;
   mpiw = &world;
   myrank = mpiw->rank();
   if (myrank == 0) {
      set_workdir(argc, argv);
      std::cout << "Parallelization using MPI: Running on " << mpiw->size() << " processors." << std::endl << std::endl;
      run_nrg_master();
   } else
      run_nrg_slave();
#else
   set_workdir(argc, argv);
   std::cout << "No MPI: single node calculation." << std::endl << std::endl;
   run_nrg_master();
#endif  // NRG_MPI
}
