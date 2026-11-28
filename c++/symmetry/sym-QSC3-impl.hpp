#include <complex>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

template<typename SC>
class SymmetryQSC3 : public SymC3<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSC3(const Params &P) : SymC3<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1,0)) {
     initInvar({
        {"Q", additive},  // charge
        {"SS", additive}, // spin
        {"P", mod3}       // C_3 rep
     });
   }

  // Multiplicity of the I=(Q,SS,P) subspace = 2S+1 = SS.
  size_t mult(const Invar &I) const override {
    return I.get("SS"); // spin multiplicity
  }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && c3_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 3);
#include "qsc3/qsc3-In2.dat"
#include "qsc3/qsc3-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "P", 0);
    return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    check_diff(Ip, I1, "P", 0); // only P=0 implemented
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
  }

  DECL;
};

template<typename SC>
void SymmetryQSC3<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  my_assert(P.channels == 3);
  int ss = I.get("SS");
#undef Complex
#define Complex(x, y) cmpl(x, y)
#define sqrt(x) csqrt(x)
#include "qsc3/qsc3-diag.dat"
#include "qsc3/qsc3-offdiag.dat"
#undef sqrt
}

}

#include "nrg-recalc-QSC3.hpp"
