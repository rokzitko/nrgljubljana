#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetryQSLR : public SymLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSLR(const Params &P) : SymLR<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1,1)) {
     initInvar({
        {"Q", additive},      // charge
        {"SS", additive},     // spin
        {"P", multiplicative} // parity
     });
   }

  // Multiplicity of the I=(Q,SS,P) subspace = 2S+1 = SS.
  size_t mult(const Invar &I) const override {
    return I.get("SS"); // spin multiplicity
  }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "qslr/qslr-2ch-In2.dat"
#include "qslr/qslr-2ch-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

template<typename SC>
void SymmetryQSLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ss = I.get("SS");
#include "qslr/qslr-2ch-diag.dat"
#include "qslr/qslr-2ch-offdiag.dat"
}

}

#include "nrg-recalc-QSLR.hpp"
