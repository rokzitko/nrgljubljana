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

namespace NRG {

template<scalar S>
class DiagSerial : public DiagEngine<S> {
public:
   DiagInfo<S> diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                                const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) {
     DiagInfo<S> diagnew;
     for (const auto &I: tasks) {
       auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P);
       nrglog('(', "[serial] Diagonalizing " << I << " dim=" << dim(h));
       auto e = diagonalise(h, DP, -1); // -1 = not using MPI
       diagnew[I] = Eigen(std::move(e), step);
     }
     return diagnew;
   }
};

}

#endif
