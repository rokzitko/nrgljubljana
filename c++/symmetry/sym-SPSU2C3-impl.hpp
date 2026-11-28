#include <complex>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

template<typename SC>
class SymmetrySPSU2C3 : public SymC3<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPSU2C3(const Params &P) : SymC3<SC>(P, {td::su2_z_square("<Sz^2>", "SS")}, Invar(1,0)) {
     initInvar({
        {"SS", additive}, // spin
        {"P", mod3}       // C_3 rep
     });
   }

  size_t mult(const Invar &I) const override { return I.get("SS"); }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) && c3_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 3);
#include "spsu2c3/spsu2c3-In2.dat"
#include "spsu2c3/spsu2c3-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
  }

  DECL;
  bool uses_pairing() const override { return true; }
};

template<typename SC>
void SymmetrySPSU2C3<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  my_assert(P.channels == 3);
  int ss = I.get("SS");
#undef Complex
#define Complex(x, y) cmpl(x, y)
#define sqrt(x) csqrt(x)
#include "spsu2c3/spsu2c3-diag.dat"
#include "spsu2c3/spsu2c3-offdiag.dat"
#include "spsu2c3/spsu2c3-isospinx.dat"
#undef sqrt
}

}

#include "nrg-recalc-SPSU2C3.hpp"
