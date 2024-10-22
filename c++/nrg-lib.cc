#include "nrg-general.hpp"

// This is compiled n the library only. Should not be used if a cblas library is available.
#ifdef CBLAS_WORKAROUND
 #define ADD_
 #include "cblas_globals.c"
 #include "cblas_dgemm.c"
 #include "cblas_zgemm.c"
 #include "cblas_xerbla.c"
#endif

namespace NRG {

void print_about_message() {
  fmt::print("NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print("Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print("Compiled on {} at {}\n\n", __DATE__, __TIME__);
}

// Called from the NRG stand-alone executable
void run_nrg_master(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, std::unique_ptr<Workdir> workdir) {
//  MPI_diag mpi(mpienv, mpiw);
  const bool embedded = false;
  if (complex_data())
    NRG_calculation<std::complex<double>> calc(std::move(workdir), embedded);
  else
    NRG_calculation<double> calc(std::move(workdir), embedded);
//  mpi.done();
}

// Called from a third-party application
void run_nrg_master(const std::string &dir) {
//  boost::mpi::environment mpienv;
//  boost::mpi::communicator mpiw;
//  MPI_diag mpi(mpienv, mpiw);
  auto workdir = set_workdir(dir);
  const bool embedded = true;
  if (complex_data())
    NRG_calculation<std::complex<double>> calc(std::move(workdir), embedded);
  else
    NRG_calculation<double> calc(std::move(workdir), embedded);
}

template<scalar S>
void run_nrg_slave_impl(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw) {
  DiagMPI<S> eng(mpienv, mpiw);
  constexpr auto master = 0;
  DiagParams DP;
  for (;;) {
    if (mpiw.iprobe(master, boost::mpi::any_tag)) { // message can be received.
      int task = 0;
      const auto status = mpiw.recv(master, boost::mpi::any_tag, task);
      mpilog("Slave " << mpiw.rank() << " received message with tag " << status.tag());
      check_status(status);
      switch (status.tag()) {
      case TAG_SYNC:
        DP = eng.receive_params();
        break;
      case TAG_DIAG:
        eng.slave_diag(master, DP);
        break;
      case TAG_EXIT:
        return; // exit from run_slave()
      default:
        std::cout << "MPI error: unknown tag on " << mpiw.rank() << std::endl;
        break;
      }
    } else usleep(100); // sleep to reduce the load on the computer. (OpenMPI "feature" workaround)
  }
}

void run_nrg_slave(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw) {
  if (complex_data())
    run_nrg_slave_impl<std::complex<double>>(mpienv, mpiw);
  else
    run_nrg_slave_impl<double>(mpienv, mpiw);
}
    
} // namespace
