#include <cmath>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;
using std::sqrt;

template<typename SC>
class SymmetryQSZLR : public SymFieldLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSZLR(const Params &P) : SymFieldLR<SC>(P, {td::sz_square("<Sz^2>"), td::sz("<Sz>"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,0,1)) {
     initInvar({
        {"Q", additive},      // charge
        {"SSZ", additive},    // spin projection
        {"P", multiplicative} // parity
     });
   }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "qszlr/qszlr-2ch-In2.dat"
#include "qszlr/qszlr-2ch-QN.dat"
  }

  DECL;
};

template<typename SC>
void SymmetryQSZLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
#include "qszlr/qszlr-2ch-offdiag.dat"
#include "qszlr/qszlr-2ch-diag.dat"
}

}

#include "nrg-recalc-QSZLR.hpp"
