#define NRG_EXECUTABLE

#include "nrg-general.h" // common
#include "nrg-lib.h"     // exposed in library
#include "nrg.h"         // specific to executable
#include "openmp.h"      // report_openMP() called from main()
#include "workdir.h"

#ifdef NRG_MPI
mpi::environment *mpienv;
mpi::communicator *mpiw;
#endif

inline void help(int argc, char **argv, std::string help_message)
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() >= 1 && args[0] == "-h") {
    std::cout << help_message << std::endl;
    exit(EXIT_SUCCESS);
  }
}

Workdir set_workdir(int argc, char **argv) { // not inline!
  std::string dir = default_workdir; // defined in workdir.h
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() == 2 && args[0] == "-w") dir = args[1];
  return Workdir(dir);
}

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
    help(argc, argv, "Usage: nrg [-h] [-w workdir]");
    auto workdir = set_workdir(argc, argv);
    const bool embedded = false;
    run_nrg_master(workdir, embedded);
#ifdef NRG_MPI
  } else
    run_nrg_slave(); // slaves do no disk I/O
#endif
}
