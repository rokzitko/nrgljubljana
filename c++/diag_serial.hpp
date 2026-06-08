#ifndef _diag_serial_hpp_
#define _diag_serial_hpp_

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
#include "subspaces.hpp"

namespace NRG {

template<scalar S>
class DiagSerial : public DiagEngine<S> {
public:
   DiagInfo<S> diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                                const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) {
     DiagInfo<S> diagnew;
     const auto tasks_by_size = tasks_descending_by_subspace_dimension(tasks, diagprev, Sym);
     for (const auto &[dim, I]: tasks_by_size) {
       auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P);
       nrglog('(', "[serial] Diagonalizing " << I << " dim=" << dim);
       auto e = diagonalise(h, DP, -1); // -1 = not using MPI
       diagnew[I] = Eigen(std::move(e), step);
     }
     return diagnew;
   }
};

}

#endif
