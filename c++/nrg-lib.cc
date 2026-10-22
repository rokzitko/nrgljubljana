#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <system_error>
#include <stdexcept>
#include <exception>
#include <boost/io/ios_state.hpp>
#include <gmp.h>

#include "nrg-general.hpp"

namespace NRG {

void print_about_message() {
  fmt::print("NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print("Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print("Compiled on {} at {}\n\n", __DATE__, __TIME__);
}

void clear_done_marker() {
  std::error_code ec;
  static_cast<void>(std::filesystem::remove("DONE", ec));
  if (ec) throw std::runtime_error("Can't remove stale DONE: " + ec.message());
}

template <scalar S>
void run_nrg_master_impl(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, std::unique_ptr<Workdir> workdir) {
  std::shared_ptr<DiagEngine<S>> eng = std::make_shared<DiagMPI<S>>(mpienv, mpiw);
  clear_done_marker();
  const bool embedded = false;
  NRG_calculation<S> calc(std::move(workdir), eng, embedded);
}

// Called from the NRG stand-alone executable (MPI version)
void run_nrg_master(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, std::unique_ptr<Workdir> workdir) {
  if (complex_data())
    run_nrg_master_impl<std::complex<double>>(mpienv, mpiw, std::move(workdir));
  else
    run_nrg_master_impl<double>(mpienv, mpiw, std::move(workdir));
}

template <scalar S>
void run_nrg_master_impl(std::unique_ptr<Workdir> workdir) {
  std::shared_ptr<DiagEngine<S>> eng = std::make_shared<DiagOpenMP<S>>();
  const bool embedded = true;
  NRG_calculation<S> calc(std::move(workdir), eng, embedded);
}

// Called from a third-party application
void run_nrg_master(const std::string &dir) {
  run_nrg_master(dir, WorkdirMode::unique_temporary);
}

void run_nrg_master(const std::string &dir, const WorkdirMode mode) {
  std::ostringstream diagnostics;
  if (!prepare_parallel_runtime(diagnostics))
    throw std::runtime_error("NRG numerical runtime validation failed:\n" + diagnostics.str());
  auto workdir = set_workdir(dir, mode);
  clear_done_marker();
  if (complex_data())
    run_nrg_master_impl<std::complex<double>>(std::move(workdir));
  else
    run_nrg_master_impl<double>(std::move(workdir));
}

void run_nrg(const std::string &calculation_directory, const std::string &scratch_directory) {
  static std::mutex invocation_mutex;
  const std::lock_guard lock(invocation_mutex);
  if (calculation_directory.empty()) throw std::invalid_argument("NRG calculation_directory must not be empty");
  if (!std::filesystem::path(scratch_directory).is_absolute())
    throw std::invalid_argument("NRG scratch_directory must be an absolute path");
  if (!std::filesystem::is_directory(scratch_directory))
    throw std::invalid_argument("NRG scratch_directory must be an existing directory");

  struct SavedState {
    const mp_bitcnt_t precision = mpf_get_default_prec();
    boost::io::ios_base_all_saver cout_format{std::cout}, cerr_format{std::cerr};
    boost::io::ios_fill_saver cout_fill{std::cout}, cerr_fill{std::cerr};
    ~SavedState() { mpf_set_default_prec(precision); }
  } saved_state;

  const auto original_directory = std::filesystem::current_path();
  try {
    std::filesystem::current_path(calculation_directory);
    run_nrg_master(scratch_directory);
  } catch (...) {
    std::error_code ec;
    std::filesystem::current_path(original_directory, ec);
    if (ec)
      std::throw_with_nested(std::filesystem::filesystem_error("Can't restore caller's current directory", original_directory, ec));
    throw;
  }
  std::filesystem::current_path(original_directory);
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
