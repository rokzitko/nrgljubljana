class SymmetryPP : public Symmetry {
  public:
  SymmetryPP() : Symmetry() { all_syms["PP"] = this; }

  void init() override {
    InvarStructure InvStruc[] = {
       {"Pa", multiplicative}, // fermion parity in channel a
       {"Pb", multiplicative}  // fermion parity in channel b
    };
    initInvar(InvStruc, ARRAYLENGTH(InvStruc));
    InvarSinglet = Invar(1, 1);
  }

  void load() override {
    switch (channels) {

      case 2:
#include "pp/pp-2ch-In2.dat"
#include "pp/pp-2ch-QN.dat"
        break;

      default: my_assert_not_reached();
    }
  }

  void makematrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);
  void makematrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);

  void calculate_TD(const Step &step, const DiagInfo &diag, double factor) override{};

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) override {
    return z2_equality(I1.get("Pa"), I2.get("Pa"), I3.get("Pa")) && z2_equality(I1.get("Pb"), I2.get("Pb"), I3.get("Pb"));
  }

  DECL;
  HAS_GLOBAL;
};

Symmetry *SymPP = new SymmetryPP;

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP
#undef OFFDIAG_AN_DO
#undef OFFDIAG_AN_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(i, j, ch, 1, t_matel(factor) * xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_AN_DO(i, j, ch, factor) offdiag_function(i, j, ch, 2, t_matel(factor) * xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_AN_UP(i, j, ch, factor) offdiag_function(i, j, ch, 3, t_matel(factor) * xi(step.N(), ch), h, qq, In, opch)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(step.N() + 1, ch), h, qq)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(step.N() + 1, ch), h, qq)

void SymmetryPP::makematrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  switch (channels) {

    case 2:
#include "pp/pp-2ch-offdiag-CR-UP.dat"
#include "pp/pp-2ch-offdiag-CR-DO.dat"
#include "pp/pp-2ch-offdiag-AN-UP.dat"
#include "pp/pp-2ch-offdiag-AN-DO.dat"
#include "pp/pp-2ch-diag.dat"
#include "pp/pp-2ch-Ixtot.dat"
      break;

    default: my_assert_not_reached();
  }
}

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP
#undef OFFDIAG_AN_DO
#undef OFFDIAG_AN_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * xiDOWN(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(i, j, ch, 1, t_matel(factor) * xiUP(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_AN_DO(i, j, ch, factor) offdiag_function(i, j, ch, 2, t_matel(factor) * xiDOWN(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_AN_UP(i, j, ch, factor) offdiag_function(i, j, ch, 3, t_matel(factor) * xiUP(step.N(), ch), h, qq, In, opch)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(step.N() + 1, ch), h, qq)

#undef DIAG_UP
#define DIAG_UP(i, j, ch, number) diag_function(i, ch, number, zetaUP(step.N() + 1, ch), h, qq)

#undef DIAG_DOWN
#define DIAG_DOWN(i, j, ch, number) diag_function(i, ch, number, zetaDOWN(step.N() + 1, ch), h, qq)

void SymmetryPP::makematrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  switch (channels) {

    case 2:
#include "pp/pp-2ch-offdiag-CR-UP.dat"
#include "pp/pp-2ch-offdiag-CR-DO.dat"
#include "pp/pp-2ch-offdiag-AN-UP.dat"
#include "pp/pp-2ch-offdiag-AN-DO.dat"
#include "pp/pp-2ch-diag-UP.dat"
#include "pp/pp-2ch-diag-DOWN.dat"
#include "pp/pp-2ch-Ixtot.dat"
      break;

    default: my_assert_not_reached();
  }
}

void SymmetryPP::makematrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (P.polarized)
    makematrix_polarized(h, step, qq, I, In, opch);
  else
    makematrix_nonpolarized(h, step, qq, I, In, opch);
}

#include "nrg-recalc-PP.cc"
