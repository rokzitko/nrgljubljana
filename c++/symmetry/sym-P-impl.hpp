#include <string>
#include <vector>

#include "abelian-cran.hpp"

namespace NRG {

// Fermion parity: Z_2 symmetry for problems with superconducting pairing and arbitrary spin interactions. Any number
// of channels; the coefficient tables are generated at run time (abelian-cran.hpp).
template<typename SC>
class SymmetryP : public SymmetryCRAN<SC> {
  public:
   SymmetryP(const Params &P) : SymmetryCRAN<SC>(P, AbelianShell::Kind::P, Invar(1), Invar(-1), true) {
     initInvar({
       {"P", multiplicative} // fermion parity
     });
   }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override { return z2_equality(I1.get("P"), I2.get("P"), I3.get("P")); }
};

}
