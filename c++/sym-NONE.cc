class SymmetryNONE : public Symmetry {
  public:
  SymmetryNONE() : Symmetry() { all_syms["NONE"] = this; }

  void init() override {
    InvarStructure InvStruc[] = {
       {"x", additive} // dummy quantum number
    };
    initInvar(InvStruc, ARRAYLENGTH(InvStruc));
    InvarSinglet = Invar(0);
  }

  void load() override {
    switch (channels) {
      case 1:
#include "none/none-1ch-In2.dat"
#include "none/none-1ch-QN.dat"
        break;

      case 2:
#include "none/none-2ch-In2.dat"
#include "none/none-2ch-QN.dat"
        break;

      default: my_assert_not_reached();
    }
  }

  void makematrix_polarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In);
  void makematrix_nonpolarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In);

  void calculate_TD(const DiagInfo &diag, double factor) override{};

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

Symmetry *SymNONE = new SymmetryNONE;

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * xi(STAT::N, ch), h, qq, In)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(i, j, ch, 1, t_matel(factor) * xi(STAT::N, ch), h, qq, In)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(STAT::N + 1, ch), h, qq)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(STAT::N + 1, ch), h, qq)

void SymmetryNONE::makematrix_nonpolarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In) {
  switch (channels) {
    case 1:
#include "none/none-1ch-offdiag-CR-UP.dat"
#include "none/none-1ch-offdiag-CR-DO.dat"
#include "none/none-1ch-diag.dat"
#include "none/none-1ch-Ixtot.dat"
      break;

    case 2:
#include "none/none-2ch-offdiag-CR-UP.dat"
#include "none/none-2ch-offdiag-CR-DO.dat"
#include "none/none-2ch-diag.dat"
#include "none/none-2ch-Ixtot.dat"
      break;

    default: my_assert_not_reached();
  }
}

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * xiDOWN(STAT::N, ch), h, qq, In)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(i, j, ch, 1, t_matel(factor) * xiUP(STAT::N, ch), h, qq, In)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(STAT::N + 1, ch), h, qq)

#undef DIAG_UP
#define DIAG_UP(i, j, ch, number) diag_function(i, ch, number, zetaUP(STAT::N + 1, ch), h, qq)

#undef DIAG_DOWN
#define DIAG_DOWN(i, j, ch, number) diag_function(i, ch, number, zetaDOWN(STAT::N + 1, ch), h, qq)

void SymmetryNONE::makematrix_polarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In) {
  switch (channels) {
    case 1:
#include "none/none-1ch-offdiag-CR-UP.dat"
#include "none/none-1ch-offdiag-CR-DO.dat"
#include "none/none-1ch-diag-UP.dat"
#include "none/none-1ch-diag-DOWN.dat"
#include "none/none-1ch-Ixtot.dat"
      break;

    case 2:
#include "none/none-2ch-offdiag-CR-UP.dat"
#include "none/none-2ch-offdiag-CR-DO.dat"
#include "none/none-2ch-diag-UP.dat"
#include "none/none-2ch-diag-DOWN.dat"
#include "none/none-2ch-Ixtot.dat"
      break;

    default: my_assert_not_reached();
  }
}

void SymmetryNONE::makematrix(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In) {
  if (P::polarized) {
    makematrix_polarized(h, qq, I, In);
  } else {
    makematrix_nonpolarized(h, qq, I, In);
  }
}

#include "nrg-recalc-NONE.cc"
