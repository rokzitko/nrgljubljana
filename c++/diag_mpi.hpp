#ifndef _diag_mpi_hpp_
#define _diag_mpi_hpp_

#include <iostream>
#include <list>
#include <deque>
#include <algorithm>

#ifdef OMPI_SKIP_MPICXX // workaround to avoid warnings for for redefinition in mpi/environment.hpp
 #undef OMPI_SKIP_MPICXX
#endif
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp> // broadcast

#include "traits.hpp"
#include "invar.hpp"
#include "step.hpp"
#include "operators.hpp"
#include "symmetry.hpp"
#include "params.hpp"
#include "diag.hpp"
#include "misc.hpp"
#include "core.hpp"

namespace NRG {

enum TAG : int { TAG_EXIT = 1, TAG_DIAG_DBL, TAG_DIAG_CMPL, TAG_SYNC, TAG_MATRIX, TAG_INVAR,
                 TAG_MATRIX_SIZE, TAG_MATRIX_LINE, TAG_EIGEN_INT, TAG_EIGEN_VEC };

class MPI_diag {
 private:
   boost::mpi::environment &mpienv;
   boost::mpi::communicator &mpiw;

 public:
   MPI_diag(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw) : mpienv(mpienv), mpiw(mpiw) {}
   auto myrank() { return mpiw.rank(); } // used in diag.h, time_mem.h
   void send_params(const DiagParams &DP) {
     mpilog("Sending diag parameters " << DP.diag << " " << DP.diagratio);
     for (auto i = 1; i < mpiw.size(); i++) mpiw.send(i, TAG_SYNC, 0);
     auto DPcopy = DP;
     boost::mpi::broadcast(mpiw, DPcopy, 0);
   }
   DiagParams receive_params() {
     DiagParams DP;
     boost::mpi::broadcast(mpiw, DP, 0);
     mpilog("Received diag parameters " << DP.diag << " " << DP.diagratio);
     return DP;
   }
   void check_status(const boost::mpi::status &status) {
     if (status.error()) {
       std::cout << "MPI communication error. rank=" << mpiw.rank() << std::endl;
       mpienv.abort(1);
     }
   }
   // NOTE: MPI is limited to message size of 2GB (or 4GB). For big problems we thus need to send objects line by line.
    template<scalar S>
    void send_matrix(const int dest, const EigenMatrix<S> &m) {
      const auto size1 = NRG::size1(m);
      mpiw.send(dest, TAG_MATRIX_SIZE, size1);
      const auto size2 = NRG::size2(m);
      mpiw.send(dest, TAG_MATRIX_SIZE, size2);
      mpilog("Sending matrix of size " << size1 << " x " << size2 << " line by line to " << dest);
      for (const auto i: range0(size1)) {
        mpiw.send(dest, TAG_MATRIX_LINE, m.row(i));
      }
   }
   template<scalar S> auto receive_matrix(const int source) {
      size_t size1;
      check_status(mpiw.recv(source, TAG_MATRIX_SIZE, size1));
      size_t size2;
      check_status(mpiw.recv(source, TAG_MATRIX_SIZE, size2));
      EigenMatrix<S> m(size1, size2);
      mpilog("Receiving matrix of size " << size1 << " x " << size2 << " line by line from " << source);
      for (const auto i: range0(size1)) {
        EigenVector<S> vec;
        check_status(mpiw.recv(source, TAG_MATRIX_LINE, vec));
        my_assert(vec.size() == size2);
        m.row(i) = vec;
      }
      return m;
   }
   template<scalar S> void send_raweigen(const int dest, const RawEigen<S> &eig) {
     mpilog("Sending eigen from " << mpiw.rank() << " to " << dest);
     mpiw.send(dest, TAG_EIGEN_VEC, eig.val);
     send_matrix<S>(dest, eig.vec);
   }
   template<scalar S> auto receive_raweigen(const int source) {
     mpilog("Receiving eigen from " << source << " on " << mpiw.rank());
     RawEigen<S> eig;
     check_status(mpiw.recv(source, TAG_EIGEN_VEC, eig.val));
     eig.vec = receive_matrix<S>(source);
     return eig;
   } 
   // Read results from a slave process.
   template<scalar S> std::pair<Invar, RawEigen<S>> read_from(const int source) {
     mpilog("Reading results from " << source);
     const auto eig = receive_raweigen<S>(source);
     Invar Irecv;
     check_status(mpiw.recv(source, TAG_INVAR, Irecv));
     mpilog("Received results for subspace " << Irecv << " [nr=" << eig.getnrcomputed() << ", dim=" << eig.getdim() << "]");
     return {Irecv, eig};
   }
   // Handle a diagonalisation request
   template<scalar S> void slave_diag(const int master, const DiagParams &DP) {
     // 1. receive the matrix and the subspace identification
     auto m = receive_matrix<S>(master);
     Invar I;
     check_status(mpiw.recv(master, TAG_INVAR, I));
     // 2. preform the diagonalisation
     const auto eig = diagonalise(m, DP, myrank());
     // 3. send back the results
     send_raweigen<S>(master, eig);
     mpiw.send(master, TAG_INVAR, I);
   }
   template<scalar S>
   DiagInfo<S> diagonalisations_MPI(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                                    const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) {
       DiagInfo<S> diagnew;
       send_params(DP);                                         // Synchronise parameters
       std::list<Invar> tasks_todo(tasks.begin(), tasks.end());
       std::list<Invar> tasks_done;
       std::deque<int> nodes_available(mpiw.size());            // Available nodes including the master, which is always at the head of the deque
       std::iota(nodes_available.begin(), nodes_available.end(), 0);
       nrglog('M', "nrtasks=" << tasks_todo.size() << " nrnodes=" << nodes_available.size());
       while (!tasks_todo.empty()) {
         my_assert(!nodes_available.empty());
         // i is the node to which the next job will be scheduled. (If a single task is left undone, do it on the master
         // node to avoid the unnecessary network copying.)
         const auto i = tasks_todo.size() != 1 ? get_back(nodes_available) : 0;
         // On master, we take short jobs from the end. On slaves, we take long jobs from the beginning.
         const Invar I = i == 0 ? get_back(tasks_todo) : get_front(tasks_todo);
         auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P); // non-const
         nrglog('M', "Scheduler: job " << I << " (dim=" << dim(h) << ")" << " on node " << i);
         if (i == 0) {
           // On master, diagonalize immediately.
           auto e = diagonalise(h, DP, myrank());
           diagnew[I] = Eigen<S>(std::move(e), step);
           tasks_done.push_back(I);
           nodes_available.push_back(0);
         } else {
           mpiw.send(i, std::is_same_v<S, double> ? TAG_DIAG_DBL : TAG_DIAG_CMPL, 0);
           send_matrix<S>(i, h);
           mpiw.send(i, TAG_INVAR, I);
         }
         // Check for terminated jobs
         while (auto status = mpiw.iprobe(boost::mpi::any_source, TAG_EIGEN_VEC)) {
           nrglog('M', "Receiveing results from " << status->source());
           auto [Irecv, eig] = read_from<S>(status->source());
           diagnew[Irecv] = Eigen<S>(std::move(eig), step);
           tasks_done.push_back(Irecv);
           // The node is now available for new tasks!
           nodes_available.push_back(status->source());
         }
       }
       // Keep reading results sent from the slave processes until all tasks have been completed.
       while (tasks_done.size() != tasks.size()) {
         const auto status = mpiw.probe(boost::mpi::any_source, TAG_EIGEN_VEC);
         auto [Irecv, eig]  = read_from<S>(status.source());
         diagnew[Irecv] = Eigen<S>(std::move(eig), step);
         tasks_done.push_back(Irecv);
       }
       return diagnew;
     }
   void done() {
     for (auto i = 1; i < mpiw.size(); i++) mpiw.send(i, TAG_EXIT, 0); // notify slaves we are done
   }
};

} // namespace

#endif
