#include <cmath>
#include <cstddef>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;
using std::sqrt;

template<typename SC>
class SymmetryISOSZLR : public SymFieldLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISOSZLR(const Params &P) : SymFieldLR<SC>(P, {td::sz_square("<Sz^2>"), td::sz("<Sz>"), td::su2_z_square("<Q^2>", "II", 4)}, Invar(1,0,1)) {
     initInvar({
        {"II", additive},     // isospin
        {"SSZ", additive},    // spin projection
        {"P", multiplicative} // parity
     });
   }
  // Multiplicity of the I=(II,SSZ) subspace = (2I+1) = II.
  size_t mult(const Invar &I) const override { return I.get("II"); }
  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }
  void load() override {
    my_assert(P.channels == 2);
#include "isoszlr/isoszlr-2ch-In2.dat"
#include "isoszlr/isoszlr-2ch-QN.dat"
  }
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_abs_diff(Ip, I1, "SSZ", 1);
    return su2_doublet_factor(Ip.get("II"), I1.get("II"));
  }
  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

template<typename SC>
void SymmetryISOSZLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ii = I.get("II");
#include "isoszlr/isoszlr-2ch-offdiag.dat"
}

}

#include "nrg-recalc-ISOSZLR.hpp"
