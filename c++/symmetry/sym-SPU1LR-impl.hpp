#include <cmath>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;
using std::sqrt;

template<typename SC>
class SymmetrySPU1LR : public SymFieldLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPU1LR(const Params &P) : SymFieldLR<SC>(P, {td::sz_square("<Sz^2>"), td::sz("<Sz>")}, Invar(0,1)) {
     initInvar({
        {"SSZ", additive},    // spin projection
        {"P", multiplicative} // parity
     });
   }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    switch (P.channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-In2.dat"
#include "spu1lr/spu1lr-1ch-QN.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-In2.dat"
#include "spu1lr/spu1lr-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  bool uses_pairing() const override { return true; }
  bool uses_anomalous_hopping() const override { return true; }
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_GLOBAL;
};

template<typename SC>
void SymmetrySPU1LR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  switch (P.channels) {
    case 1:
#include "spu1lr/spu1lr-1ch-offdiag.dat"
#include "spu1lr/spu1lr-1ch-anomalous.dat"
#include "spu1lr/spu1lr-1ch-diag.dat"
#include "spu1lr/spu1lr-1ch-isospinx.dat"
      break;
    case 2:
#include "spu1lr/spu1lr-2ch-diag.dat"
#include "spu1lr/spu1lr-2ch-offdiag.dat"
#include "spu1lr/spu1lr-2ch-anomalous.dat"
#include "spu1lr/spu1lr-2ch-isospinx.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-SPU1LR.hpp"
