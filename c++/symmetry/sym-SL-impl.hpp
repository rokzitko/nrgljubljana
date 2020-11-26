namespace NRG {

template<typename SC>
class SymmetrySL : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySL(const Params &P) : Symmetry<SC>(P, std::vector{"<Q>", "<Q^2>", "<sQ^2>"}, Invar(0)) {
     initInvar({
        {"Q", additive} // charge
     });
   }

  void load() override {
    switch (P.channels) {
      case 1:
#include "sl/sl-1ch-In2.dat"
#include "sl/sl-1ch-QN.dat"
        break;

      case 2:
#include "sl/sl-2ch-In2.dat"
#include "sl/sl-2ch-QN.dat"
        break;

      case 3:
#include "sl/sl-3ch-In2.dat"
#include "sl/sl-3ch-QN.dat"
        break;

      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const int q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
    }
    stats.td.set("<Q>",    trQ / stats.Z);
    stats.td.set("<Q^2>",  trQ2 / stats.Z);
    // charge fluctuations -> susceptibility
    stats.td.set("<sQ^2>", (trQ2 / stats.Z) - pow(trQ / stats.Z, 2));
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetrySL<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  switch (P.channels) {
    case 1:
#include "sl/sl-1ch-offdiag.dat"
#include "sl/sl-1ch-diag.dat"
      break;
    case 2:
#include "sl/sl-2ch-offdiag.dat"
#include "sl/sl-2ch-diag.dat"
      break;
    case 3:
#include "sl/sl-3ch-offdiag.dat"
#include "sl/sl-3ch-diag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-SL.hpp"
