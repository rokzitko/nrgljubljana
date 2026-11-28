#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetrySU2 : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySU2(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Q^2>", "II", 4)}, Invar(1)) {
       initInvar({
         {"II", additive} // isospin
       });
     }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"));
  }

  // Multiplicity of the I=(II) subspace = (2I+1) = II.
  size_t mult(const Invar &I) const override {
    return I.get("II"); // isospin multiplicity
  }

  // We always must have I >= 0.
  bool Invar_allowed(const Invar &I) const override { return I.get("II") > 0; }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_doublet_factor(Ip.get("II"), I1.get("II"));
  }

  void load() override {
    switch (P.channels) {
      case 1:
#include "su2/su2-1ch-In2.dat"
#include "su2/su2-1ch-QN.dat"
        break;
      case 2:
#include "su2/su2-2ch-In2.dat"
#include "su2/su2-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

// For SU2, the <||f||> depend on the "type"!
// OFFDIAG_1 corresponds to the first combination, OFFDIAG_2 to the
// second combination of operators. Each of them has contributions
// in each channel.

template<typename SC>
void SymmetrySU2<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ii = I.get("II");
  int NN   = step.getnn();
  switch (P.channels) {
    case 1:
#include "su2/su2-1ch-offdiag-1.dat"
#include "su2/su2-1ch-offdiag-2.dat"
      break;
    case 2:
#include "su2/su2-2ch-offdiag-1.dat"
#include "su2/su2-2ch-offdiag-2.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-SU2.hpp"
