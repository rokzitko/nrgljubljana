#include <string>
#include <vector>

#include "abelian-cran.hpp"

namespace NRG {

// No symmetry: a single invariant subspace. Any number of channels; the coefficient tables are generated at run time
// (abelian-cran.hpp). In the combination order of NONE the annihilation part of the hopping does not contribute to the
// upper triangle of the Hamiltonian, so only the creation-operator blocks of <||f||> are used and recalculated.
template<typename SC>
class SymmetryNONE : public SymmetryCRAN<SC> {
  public:
   SymmetryNONE(const Params &P) : SymmetryCRAN<SC>(P, AbelianShell::Kind::NONE, Invar(0), Invar(0), false) {
     initInvar({
       {"x", additive} // dummy quantum number
     });
   }

};

}
