#ifndef _diag_openmp_hpp_
#define _diag_openmp_hpp_

#include "traits.hpp"
#include "step.hpp"
#include "operators.hpp"
#include "coef.hpp"
#include "eigen.hpp"
#include "output.hpp"
#include "invar.hpp"
#include "params.hpp"
#include "symmetry.hpp"

namespace NRG {

template<scalar S>
auto diagonalisations_OpenMP(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                             const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) {
  DiagInfo<S> diagnew;
  const auto nr = tasks.size();
  size_t itask = 0;
  // cppcheck-suppress unreadVariable symbolName=nth
  const int nth = P.diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
  for (itask = 0; itask < nr; itask++) {
    const Invar I  = tasks[itask];
    auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P); // non-const, consumed by diagonalise()
    const int thid = omp_get_thread_num();
#pragma omp critical
    { nrglog('(', "Diagonalizing " << I << " dim=" << dim(h) << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
    auto e = diagonalise(h, DP, -1); // -1 = not using MPI
#pragma omp critical
    { diagnew[I] = Eigen(std::move(e), step); }
  }
  return diagnew;
}

}

#endif
