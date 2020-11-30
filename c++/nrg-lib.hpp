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
void run_nrg_master(const std::string &workdir);

void run_nrg_master(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw, std::unique_ptr<Workdir> workdir);
void run_nrg_slave(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw);

} // namespace

#endif
