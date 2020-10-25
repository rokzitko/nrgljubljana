#include "nrg-general.hpp"

namespace NRG {

void print_about_message() {
  fmt::print(fmt::emphasis::bold, "NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print(fmt::emphasis::bold, "Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print(fmt::emphasis::bold, "Compiled on {} at {}\n\n", __DATE__, __TIME__);
}

// Called from the NRG stand-alone executable
void run_nrg_master(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, const Workdir &workdir) {
  MPI_diag mpi(mpienv, mpiw);
  const bool embedded = false;
  if (complex_data())
    NRG_calculation<std::complex<double>> calc(mpi, workdir, embedded);
  else
    NRG_calculation<double> calc(mpi, workdir, embedded);
  mpi.done();
}

// Called from a third-party application
void run_nrg_master(const std::string &dir) {
  boost::mpi::environment mpienv;
  boost::mpi::communicator mpiw;
  MPI_diag mpi(mpienv, mpiw);
  auto workdir = set_workdir(dir);
  const bool embedded = true;
  if (complex_data())
    NRG_calculation<std::complex<double>> calc(mpi, workdir, embedded);
  else
    NRG_calculation<double> calc(mpi, workdir, embedded);
}

void run_nrg_slave(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw) {
  MPI_diag mpi(mpienv, mpiw);
  constexpr auto master = 0;
  DiagParams DP;
  for (;;) {
    if (mpiw.iprobe(master, boost::mpi::any_tag)) { // message can be received.
      int task;
      const auto status = mpiw.recv(master, boost::mpi::any_tag, task);
      mpilog("Slave " << mpiw.rank() << " received message with tag " << status.tag());
      mpi.check_status(status);
      switch (status.tag()) {
      case TAG_SYNC:
        DP = mpi.receive_params();
        break;
      case TAG_DIAG_DBL:
        mpi.slave_diag<double>(master, DP);
        break;
      case TAG_DIAG_CMPL:
        mpi.slave_diag<std::complex<double>>(master, DP);
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

} // namespace