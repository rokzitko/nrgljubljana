namespace NRG {

template<typename SC>
class SymmetryNONE : public Symmetry<SC> {
  private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryNONE(const Params &P) : Symmetry<SC>(P, std::vector<std::string>{}, Invar(0)) {
     initInvar({
       {"x", additive} // dummy quantum number
     });
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
  void make_matrix_polarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const;
  void make_matrix_nonpolarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const;
  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {};
  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

#undef OFFDIAG_CR_DO
#undef OFFDIAG_CR_UP

#define OFFDIAG_CR_DO(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_CR_UP(i, j, ch, factor) offdiag_function(step, i, j, ch, 1, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryNONE<SC>::make_matrix_nonpolarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
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
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef DIAG_UP
#define DIAG_UP(i, j, ch, number) this->diag_function(step, i, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#undef DIAG_DOWN
#define DIAG_DOWN(i, j, ch, number) this->diag_function(step, i, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryNONE<SC>::make_matrix_polarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
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
void SymmetryNONE<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  if (P.polarized)
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  else
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
}

}

#include "nrg-recalc-NONE.hpp"
