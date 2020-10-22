template<typename SC>
class SymmetrySL3 : public Symmetry<SC> {
  private:
   outfield Q1, Q12, sQ12, Q2, Q22, sQ22, Q3, Q32, sQ32;
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySL3(const Params &P, Allfields &allfields) : Symmetry<SC>(P),
     Q1(P, allfields, "<Q1>", 1), Q12(P, allfields, "<Q1^2>", 2), sQ12(P, allfields, "<sQ1^2>", 3),
     Q2(P, allfields, "<Q2>", 4), Q22(P, allfields, "<Q2^2>", 5), sQ22(P, allfields, "<sQ2^2>", 6),
     Q3(P, allfields, "<Q3>", 7), Q32(P, allfields, "<Q3^2>", 8), sQ32(P, allfields, "<sQ3^2>", 9) {
       initInvar({
         {"Q1", additive}, // charge in channel 1
         {"Q2", additive}, // charge in channel 2
         {"Q3", additive}  // charge in channel 3
       });
       this->InvarSinglet = Invar(0, 0, 0);
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

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, const Stats<SC> &stats, const double factor) override {
    bucket trQ1, trQ12; // Tr[Q], Tr[Q^2]
    bucket trQ2, trQ22;
    bucket trQ3, trQ32;
    for (const auto &[I, eig]: diag) {
      const Number q1   = I.get("Q1");
      const Number q2   = I.get("Q2");
      const Number q3   = I.get("Q3");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ1 += sumZ * q1;
      trQ12 += sumZ * q1 * q1;
      trQ2 += sumZ * q2;
      trQ22 += sumZ * q2 * q2;
      trQ3 += sumZ * q3;
      trQ32 += sumZ * q3 * q3;
    }
    Q1  = trQ1 / stats.Z;
    Q12 = trQ12 / stats.Z;
    // charge fluctuations -> susceptibility
    sQ12 = (trQ12 / stats.Z) - pow(trQ1 / stats.Z, 2);
    Q2  = trQ2 / stats.Z;
    Q22 = trQ22 / stats.Z;
    // charge fluctuations -> susceptibility
    sQ22 = (trQ22 / stats.Z) - pow(trQ2 / stats.Z, 2);
    Q3  = trQ3 / stats.Z;
    Q32 = trQ32 / stats.Z;
    // charge fluctuations -> susceptibility
    sQ32 = (trQ32 / stats.Z) - pow(trQ3 / stats.Z, 2);
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
void SymmetrySL3<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) {
  switch (P.channels) {
    case 3:
#include "sl3/sl3-3ch-offdiag.dat"
#include "sl3/sl3-3ch-diag.dat"
      break;
    default: my_assert_not_reached();
  }
}

#include "nrg-recalc-SL3.h"
