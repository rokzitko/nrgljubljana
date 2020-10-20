template<typename SC>
class SymmetrySL_tmpl : public Symmetry_tmpl<SC> {
 private:
   outfield Q, Q2, sQ2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySL_tmpl(const Params &P, Allfields &allfields) : Symmetry_tmpl<SC>(P),
     Q(P, allfields, "<Q>", 1), Q2(P, allfields, "<Q^2>", 2), sQ2(P, allfields, "<sQ^2>", 3) {
       initInvar({
         {"Q", additive} // charge
       });
       this->InvarSinglet = Invar(0);
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

  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats_tmpl<SC> &stats, const double factor) override {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const Number q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
    }
    Q  = trQ / stats.Z;
    Q2 = trQ2 / stats.Z;
    // charge fluctuations -> susceptibility
    sQ2 = (trQ2 / stats.Z) - pow(trQ / stats.Z, 2);
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
void SymmetrySL_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
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

#include "nrg-recalc-SL.cc"
