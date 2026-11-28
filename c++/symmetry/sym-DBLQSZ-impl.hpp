#include <cmath>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;

template<typename SC>
class SymmetryDBLQSZ : public SymField<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryDBLQSZ(const Params &P) : SymField<SC>(P, {td::sz_square("<Sz^2>"), td::sz("<Sz>"), td::mean("<Q1>", "Q1"), td::square("<Q1^2>", "Q1"), td::mean("<Q2>", "Q2"), td::square("<Q2^2>", "Q2")}, Invar(0, 0, 0)) {
     initInvar({
        {"Q1", additive}, // charge 1
        {"Q2", additive}, // charge 2
        {"SSZ", additive}  // spin projection
     });
   }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     return u1_equality(I1.get("Q1"),  I2.get("Q1"),  I3.get("Q1"))
         && u1_equality(I1.get("Q2"),  I2.get("Q2"),  I3.get("Q2"))
         && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"));
   }

  void load() override {
    switch (P.channels) {
      case 2:
#include "dblqsz/dblqsz-2ch-In2.dat"
#include "dblqsz/dblqsz-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

template<typename SC>
void SymmetryDBLQSZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  switch (P.channels) {
    case 2:
#include "dblqsz/dblqsz-2ch-offdiag.dat"
#include "dblqsz/dblqsz-2ch-diag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-DBLQSZ.hpp"
