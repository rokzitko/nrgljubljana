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
#include "diagengine.hpp"

namespace NRG {

enum TAG : int { TAG_EXIT = 1, TAG_DIAG, TAG_SYNC, TAG_INVAR, TAG_MATRIX, TAG_MATRIX_PART, TAG_MATRIX_SIZE, TAG_VEC };

template <scalar S>
class DiagMPI : public DiagEngine<S>{
 private:
   boost::mpi::environment &mpienv;
   boost::mpi::communicator &mpiw;

 public:
   DiagMPI(boost::mpi::environment &mpienv, boost::mpi::communicator &mpiw) : mpienv(mpienv), mpiw(mpiw) {}
   ~DiagMPI() {
     for (auto i = 1; i < mpiw.size(); i++) mpiw.send(i, TAG_EXIT, 0); // notify slaves we are done
   }
   void send_params(const DiagParams &DP) {
     mpilog("Sending diag parameters: diag=" << DP.diag << " diagratio=" << DP.diagratio);
     for (auto i = 1; i < mpiw.size(); i++) mpiw.send(i, TAG_SYNC, 0);
     auto DPcopy = DP;
     boost::mpi::broadcast(mpiw, DPcopy, 0);
   }
   DiagParams receive_params() {
     DiagParams DP;
     boost::mpi::broadcast(mpiw, DP, 0);
     mpilog("Received diag parameters: diag=" << DP.diag << " diagratio=" << DP.diagratio);
     return DP;
   }
   void send_matrix(const int dest, const EigenMatrix<S> &m) {
     mpilog("send_matrix() caled, dest=" << dest);
     const size_t size1 = NRG::size1(m);
     mpiw.send(dest, TAG_MATRIX_SIZE, size1);
     const size_t size2 = NRG::size2(m);
     mpiw.send(dest, TAG_MATRIX_SIZE, size2);
     mpilog("Sending matrix of size " << size1 << " x " << size2 << " to " << dest);
//   NOTE: MPI is limited to message size of 2GB (or 4GB). For big problems we thus need to send objects line by line.
//   mpiw.send(dest, TAG_MATRIX, (S*)m.data(), size1*size2);
     for (size_t i = 0; i < size1; i++)
       mpiw.send(dest, TAG_MATRIX_PART, (S*)m.data() + i*size2, size2);
   }
   auto receive_matrix(const int source) {
     mpilog("receive_matrix() called, source=" << source);
     size_t size1;
     mpiw.recv(source, TAG_MATRIX_SIZE, size1);
     size_t size2;
     mpiw.recv(source, TAG_MATRIX_SIZE, size2);
     mpilog("Receiving matrix of size " << size1 << " x " << size2 << " from " << source);
     EigenMatrix<S> m(size1, size2);
//     mpiw.recv(source, TAG_MATRIX, (S*)m.data(), size1*size2);
     for (size_t i = 0; i < size1; i++)
       mpiw.recv(source, TAG_MATRIX_PART, (S*)m.data() + i*size2, size2);
     return m;
   }
   void send_raweigen(const int dest, const RawEigen<S> &eig) {
     mpilog("Sending eigen from " << mpiw.rank() << " to " << dest);
     mpiw.send(dest, TAG_VEC, eig.val);
     send_matrix(dest, eig.vec);
   }
   auto receive_raweigen(const int source) {
     mpilog("Receiving eigen from " << source << " on " << mpiw.rank());
     RawEigen<S> eig;
     mpiw.recv(source, TAG_VEC, eig.val);
     eig.vec = receive_matrix(source);
     return eig;
   } 
   // Read results from a slave process.
   std::pair<Invar, RawEigen<S>> read_from(const int source) {
     mpilog("Reading results from " << source);
     const auto eig = receive_raweigen(source);
     Invar Irecv;
     mpiw.recv(source, TAG_INVAR, Irecv);
     mpilog("Received results for subspace " << Irecv << " [nr=" << eig.getnrcomputed() << ", dim=" << eig.getdim() << "]");
     return {Irecv, eig};
   }
   auto myrank() { return mpiw.rank(); }
   // Handle a diagonalisation request
   void slave_diag(const int master, const DiagParams &DP) {
     // 1. receive the matrix and the subspace identification
     mpilog("slave_diag() called, master=" << master);
     Invar I;
     mpiw.recv(master, TAG_INVAR, I);
     mpilog("Received I=" << I);
     auto m = receive_matrix(master);
     // 2. preform the diagonalisation
     const auto eig = diagonalise(m, DP, myrank());
     // 3. send back the results
     send_raweigen(master, eig);
     mpiw.send(master, TAG_INVAR, I);
   }
   DiagInfo<S> diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
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
           mpiw.send(i, TAG_DIAG, 0);
           mpiw.send(i, TAG_INVAR, I);
           send_matrix(i, h);
         }
         // Check for terminated jobs
         while (auto status = mpiw.iprobe(boost::mpi::any_source, TAG_VEC)) {
           nrglog('M', "Receiveing results from " << status->source());
           auto [Irecv, eig] = read_from(status->source());
           diagnew[Irecv] = Eigen<S>(std::move(eig), step);
           tasks_done.push_back(Irecv);
           // The node is now available for new tasks!
           nodes_available.push_back(status->source());
         }
       }
       // Keep reading results sent from the slave processes until all tasks have been completed.
       while (tasks_done.size() != tasks.size()) {
         const auto status = mpiw.probe(boost::mpi::any_source, TAG_VEC);
         auto [Irecv, eig]  = read_from(status.source());
         diagnew[Irecv] = Eigen<S>(std::move(eig), step);
         tasks_done.push_back(Irecv);
       }
       return diagnew;
   }
};

} // namespace

#endif
