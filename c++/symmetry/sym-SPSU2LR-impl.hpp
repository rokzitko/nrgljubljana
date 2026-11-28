#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetrySPSU2LR : public SymLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPSU2LR(const Params &P) : SymLR<SC>(P, {td::su2_z_square("<Sz^2>", "SS")}, Invar(1,1)) {
     initInvar({
        {"SS", additive},     // spin
        {"P", multiplicative} // parity
     });
   }

  size_t mult(const Invar &I) const override { return I.get("SS"); }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "spsu2lr/spsu2lr-2ch-In2.dat"
#include "spsu2lr/spsu2lr-2ch-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
  }

  DECL;
  bool uses_pairing() const override { return true; }
  bool uses_anomalous_hopping() const override { return true; }
  HAS_DOUBLET;
  HAS_TRIPLET;
};

template<typename SC>
void SymmetrySPSU2LR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  my_assert(P.channels == 2);
  int ss = I.get("SS");
#include "spsu2lr/spsu2lr-2ch-diag.dat"
#include "spsu2lr/spsu2lr-2ch-offdiag.dat"
#include "spsu2lr/spsu2lr-2ch-anomalous.dat"
#include "spsu2lr/spsu2lr-2ch-isospinx.dat"
}

}

#include "nrg-recalc-SPSU2LR.hpp"
