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
#include "diagengine.hpp"

#ifndef NRG_ENABLE_APP_OPENMP
#define NRG_ENABLE_APP_OPENMP 0
#endif

#if NRG_ENABLE_APP_OPENMP
#include <omp.h>
#endif

namespace NRG {

template<scalar S>
class DiagOpenMP : public DiagEngine<S> {
public:
   DiagInfo<S> diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                                const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) {
      DiagInfo<S> diagnew;
      const auto nr = tasks.size();
      size_t itask = 0;
#if NRG_ENABLE_APP_OPENMP
      // cppcheck-suppress unreadVariable symbolName=nth
      const int nth = P.diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
#endif
      for (itask = 0; itask < nr; itask++) {
        const Invar I  = tasks[itask];
        auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P); // non-const, consumed by diagonalise()
#if NRG_ENABLE_APP_OPENMP
        const int thid = omp_get_thread_num();
#pragma omp critical
        { nrglog('(', "[OpenMP] Diagonalizing " << I << " dim=" << dim(h) << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
#else
        nrglog('(', "[serial scheduler] Diagonalizing " << I << " dim=" << dim(h) << " (task " << itask + 1 << "/" << nr << ")");
#endif
        auto e = diagonalise(h, DP, -1); // -1 = not using MPI
#if NRG_ENABLE_APP_OPENMP
#pragma omp critical
#endif
        { diagnew[I] = Eigen(std::move(e), step); }
      }
      return diagnew;
   }
};

}

#endif
