#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;
using std::sqrt;

template<typename SC>
class SymmetryDBLISOSZ : public SymField<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryDBLISOSZ(const Params &P) : SymField<SC>(P, {td::sz_square("<Sz^2>"), td::sz("<Sz>"), td::su2_z_square("<Q1^2>", "II1", 4), td::su2_z_square("<Q2^2>", "II2", 4)}, Invar(1, 1, 0)) {
     initInvar({
        {"II1", additive}, // isospin 1
        {"II2", additive}, // isospin 2
        {"SSZ", additive}  // spin projection
     });
   }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     return su2_triangle_inequality(I1.get("II1"), I2.get("II1"), I3.get("II1"))
       && su2_triangle_inequality(I1.get("II2"), I2.get("II2"), I3.get("II2")) && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"));
   }

  // Multiplicity of the I=(II1,II2) subspace = (2I1+1)(2I2+1) = II1 II2.
  size_t mult(const Invar &I) const override { return I.get("II1") * I.get("II2"); }

  // We always must have I1 >= 0 and I2 >= 0.
  bool Invar_allowed(const Invar &I) const override { return (I.get("II1") > 0) && (I.get("II2") > 0); }

  // This function supports isospin doublets with respect to either first or the
  // second isospin quantum number. Since the doublet operators are either doublets
  // with respect to one or another, this function should always return a non-zero
  // value!
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_abs_diff(Ip, I1, "SSZ", 1);
    return two_su2_doublet_factor(Ip, I1, "II1", "II2");
  }

  void load() override {
    switch (P.channels) {
      case 2:
#include "dblisosz/dblisosz-2ch-In2.dat"
#include "dblisosz/dblisosz-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL; // HAS_TRIPLET;
};

template<typename SC>
void SymmetryDBLISOSZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  switch (P.channels) {
    case 2:
#include "dblisosz/dblisosz-2ch-offdiag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-DBLISOSZ.hpp"
