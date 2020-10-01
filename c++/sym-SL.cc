class SymmetrySL : public Symmetry {
  private:
  outfield Q, Q2, sQ2;

  public:
  SymmetrySL() : Symmetry() { all_syms["SL"] = this; }

  void init() override {
    Q.set("<Q>", 1);
    Q2.set("<Q^2>", 2);
    sQ2.set("<sQ^2>", 3);
    InvarStructure InvStruc[] = {
       {"Q", additive} // charge
    };
    initInvar(InvStruc, ARRAYLENGTH(InvStruc));
    InvarSinglet = Invar(0);
  }

  void load() override {
    switch (channels) {
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

  void calculate_TD(const Step &step, const DiagInfo &diag, double factor) override {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]

    for (const auto &[I, eig]: diag) {
      const Number q    = I.get("Q");
      const double sumZ = calculate_Z(I, eig, factor);

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

Symmetry *SymSL = new SymmetrySL;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(step.N() + 1, ch), h, qq)

void SymmetrySL::makematrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  switch (channels) {
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
