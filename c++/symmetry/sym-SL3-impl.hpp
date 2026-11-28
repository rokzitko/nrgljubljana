#include <cmath>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::pow;

template<typename SC>
class SymmetrySL3 : public Symmetry<SC> {
  private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySL3(const Params &P) : Symmetry<SC>(P, {td::mean("<Q1>", "Q1"), td::square("<Q1^2>", "Q1"), td::variance("<sQ1^2>", "Q1"), td::mean("<Q2>", "Q2"), td::square("<Q2^2>", "Q2"), td::variance("<sQ2^2>", "Q2"), td::mean("<Q3>", "Q3"), td::square("<Q3^2>", "Q3"), td::variance("<sQ3^2>", "Q3")}, Invar(0,0,0)) {
     initInvar({
        {"Q1", additive}, // charge in channel 1
        {"Q2", additive}, // charge in channel 2
        {"Q3", additive}  // charge in channel 3
     });
   }

  void load() override {
    switch (P.channels) {
      case 3:
#include "sl3/sl3-3ch-In2.dat"
#include "sl3/sl3-3ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

template<typename SC>
void SymmetrySL3<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  switch (P.channels) {
    case 3:
#include "sl3/sl3-3ch-offdiag.dat"
#include "sl3/sl3-3ch-diag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-SL3.hpp"
