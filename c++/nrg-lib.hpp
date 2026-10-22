#ifndef _nrg_lib_hpp_
#define _nrg_lib_hpp_

// This header is included for both executable and library.

#define NRG_COMMON

#include <memory>
#include <string>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "workdir.hpp"

namespace NRG {

void print_about_message();
void clear_done_marker();

// Sequential, non-MPI embedding entry point. Read param/data and write solver
// outputs (including DONE) in calculation_directory, which must already exist.
// A relative calculation_directory is resolved from the caller's current directory.
// scratch_directory must be an existing absolute directory: a unique temporary
// child is created there and removed on return or exception, not used for outputs.
// Restore the caller's current directory, GMP default precision, and std::cout /
// std::cerr formatting (flags, precision, width, fill) on success and failure.
// The caller's original directory must remain accessible for restoration.
// Calls to this API are serialized, but the solver is NOT generally thread-safe:
// do not overlap other solver APIs or host operations using/changing the current
// directory, GMP defaults, or these streams. In particular, the lock cannot protect
// arbitrary host relative I/O. Numerical runtime setup follows the CLI, including
// initialization of an unset MKL_THREADING_LAYER from the build configuration.
void run_nrg(const std::string &calculation_directory, const std::string &scratch_directory);

// Legacy entry points: workdir selects scratch/checkpoint storage; input/output
// paths remain relative to the current directory. No embedding state isolation.
void run_nrg_master(const std::string &workdir);
void run_nrg_master(const std::string &workdir, WorkdirMode mode);

void run_nrg_master(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, std::unique_ptr<Workdir> workdir);
void run_nrg_slave(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw);

} // namespace

#endif
