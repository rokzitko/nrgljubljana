#ifndef _diagengine_hpp_
#define _diagengine_hpp_

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

template <scalar S>
class DiagEngine {
 public:
   virtual DiagInfo<S> diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                                        const std::vector<Invar> &tasks, const DiagParams &DP, const Symmetry<S> *Sym, const Params &P) = 0;
   virtual ~DiagEngine() = default;
};

}

#endif
