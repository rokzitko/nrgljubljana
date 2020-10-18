template<typename SC>
class SymmetryNONE_tmpl : public Symmetry_tmpl<SC> {
  private:
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   SymmetryNONE_tmpl(const Params &P, Allfields &allfields) : Symmetry_tmpl<SC>(P) {
     initInvar({
       {"x", additive} // dummy quantum number
     });
     this->InvarSinglet = Invar(0);
   }

  void load() override {
    switch (P.channels) {
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

  void make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef);
  void make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef);
  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats &stats, const double factor) override {};

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(step, i, j, ch, 1, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, ch, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryNONE_tmpl<SC>::make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  switch (P.channels) {
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

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xiDOWN(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(step, i, j, ch, 1, t_matel(factor) * coef.xiUP(step.N(), ch), h, qq, In, opch)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, ch, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef DIAG_UP
#define DIAG_UP(i, j, ch, number) this->diag_function(step, i, ch, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#undef DIAG_DOWN
#define DIAG_DOWN(i, j, ch, number) this->diag_function(step, i, ch, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryNONE_tmpl<SC>::make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  switch (P.channels) {
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

template<typename SC>
void SymmetryNONE_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  if (P.polarized)
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  else
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
}

#include "nrg-recalc-NONE.cc"
