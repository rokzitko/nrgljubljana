class SymmetryU1 : public Symmetry {
 private:
   outfield Q, Q2;

 public:
   template<typename ... Args> SymmetryU1(Args&& ... args) : Symmetry(std::forward<Args>(args)...),
     Q(P, allfields, "<Q>", 1), Q2(P, allfields, "<Q^2>", 2) {
       InvarStructure InvStruc[] = {
         {"Q", additive} // charge
       };
       initInvar(InvStruc, ARRAYLENGTH(InvStruc));
       InvarSinglet = Invar(0);
     }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) override { return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")); }

  void load() override {
    switch (P.channels) {
      case 1:
#include "u1/u1-1ch-In2.dat"
#include "u1/u1-1ch-QN.dat"
        break;
      case 2:
#include "u1/u1-2ch-In2.dat"
#include "u1/u1-2ch-QN.dat"
        break;
      case 3:
#include "u1/u1-3ch-In2.dat"
#include "u1/u1-3ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const Step &step, const DiagInfo &diag, const Stats &stats, double factor) override {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const Number q    = I.get("Q");
      const double sumZ = calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
    }
    Q  = trQ / stats.Z;
    Q2 = trQ2 / stats.Z;
  }

  void make_matrix_pol2x2(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef);
  void make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef);
  void make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef);

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

#undef DIAG
#define DIAG(i, ch, number) diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

#undef OFFDIAG_UP
#undef OFFDIAG_DO
#undef OFFDIAG_UPDO
#undef OFFDIAG_DOUP
#undef DIAG_UP
#undef DIAG_DOWN
#undef DIAG_DOUP

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DO(i, j, ch, factor0) offdiag_function(step, i, j, ch, 1, t_matel(factor0) * coef.xiDOWN(step.N(), ch), h, qq, In, opch)

// UPDO -> <f> from previous site for spin UP (index fnr=0)
#define OFFDIAG_UPDO(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiUPDO(step.N(), ch), h, qq, In, opch)

// DOUP -> <f> from previous site for spin DO (index fnr=1)
#define OFFDIAG_DOUP(i, j, ch, factor0) offdiag_function(step, i, j, ch, 1, t_matel(factor0) * coef.xiDOUP(step.N(), ch), h, qq, In, opch)

// Note the _half !!
#define DIAG_UP(i, j, ch, number) diag_function_half(step, i, ch, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) diag_function_half(step, i, ch, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

// Compare with ISOSPINX for symtype=SPSU2 case
// See also coefnew/u1/u1.m
#define DIAG_DOUP(i, j, ch, factor) diag_offdiag_function(step, i, j, ch, t_matel(factor) * coef.zetaDOUP(step.N() + 1, ch), h, qq)

#undef SPINZ
#define SPINZ(i, j, ch, factor) spinz_function(step, i, j, ch, t_matel(factor), h, qq)
#undef SPINX
#define SPINX(i, j, ch, factor) spinx_function(step, i, j, ch, t_matel(factor), h, qq)

void SymmetryU1::make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  switch (P.channels) {
    case 1:
#include "u1/u1-1ch-offdiag-UP.dat"
#include "u1/u1-1ch-offdiag-DO.dat"
#include "u1/u1-1ch-diag-UP.dat"
#include "u1/u1-1ch-diag-DOWN.dat"
#include "u1/u1-1ch-spinz.dat" // P.globalB
#include "u1/u1-1ch-spinx.dat" // P.globalBx
      break;
    case 2:
#include "u1/u1-2ch-offdiag-UP.dat"
#include "u1/u1-2ch-offdiag-DO.dat"
#include "u1/u1-2ch-diag-UP.dat"
#include "u1/u1-2ch-diag-DOWN.dat"
#include "u1/u1-2ch-spinz.dat" // P.globalB
#include "u1/u1-2ch-spinx.dat" // P.globalBx
      break;
    case 3:
#include "u1/u1-3ch-offdiag-UP.dat"
#include "u1/u1-3ch-offdiag-DO.dat"
#include "u1/u1-3ch-diag-UP.dat"
#include "u1/u1-3ch-diag-DOWN.dat"
#include "u1/u1-3ch-spinz.dat" // P.globalB
#include "u1/u1-3ch-spinx.dat" // P.globalBx
      break;
    default: my_assert_not_reached();
  }
}

// Full 2x2 spin matrix structure. Added 10.9.2012
void SymmetryU1::make_matrix_pol2x2(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  switch (P.channels) {
    case 1:
#include "u1/u1-1ch-offdiag-UP.dat"
#include "u1/u1-1ch-offdiag-DO.dat"
#include "u1/u1-1ch-offdiag-UPDO.dat"
#include "u1/u1-1ch-offdiag-DOUP.dat"
#include "u1/u1-1ch-diag-UP.dat"
#include "u1/u1-1ch-diag-DOWN.dat"
#include "u1/u1-1ch-diag-DOUP.dat"
#include "u1/u1-1ch-spinz.dat" // P.globalB
#include "u1/u1-1ch-spinx.dat" // P.globalBx
      break;
    case 2:
#include "u1/u1-2ch-offdiag-UP.dat"
#include "u1/u1-2ch-offdiag-DO.dat"
#include "u1/u1-2ch-offdiag-UPDO.dat"
#include "u1/u1-2ch-offdiag-DOUP.dat"
#include "u1/u1-2ch-diag-UP.dat"
#include "u1/u1-2ch-diag-DOWN.dat"
#include "u1/u1-2ch-diag-DOUP.dat"
#include "u1/u1-2ch-spinz.dat" // P.globalB
#include "u1/u1-2ch-spinx.dat" // P.globalBx
      break;
    case 3:
#include "u1/u1-3ch-offdiag-UP.dat"
#include "u1/u1-3ch-offdiag-DO.dat"
#include "u1/u1-3ch-offdiag-UPDO.dat"
#include "u1/u1-3ch-offdiag-DOUP.dat"
#include "u1/u1-3ch-diag-UP.dat"
#include "u1/u1-3ch-diag-DOWN.dat"
#include "u1/u1-3ch-diag-DOUP.dat"
#include "u1/u1-3ch-spinz.dat" // P.globalB
#include "u1/u1-3ch-spinx.dat" // P.globalBx
      break;
    default: my_assert_not_reached();
  }
}

#undef OFFDIAG_DO
#undef OFFDIAG_UP

#define OFFDIAG_DO(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)
#define OFFDIAG_UP(i, j, ch, factor) offdiag_function(step, i, j, ch, 1, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)

void SymmetryU1::make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  switch (P.channels) {
    case 1:
#include "u1/u1-1ch-offdiag-UP.dat"
#include "u1/u1-1ch-offdiag-DO.dat"
#include "u1/u1-1ch-diag.dat"
#include "u1/u1-1ch-spinz.dat" // P.globalB
#include "u1/u1-1ch-spinx.dat" // P.globalBx
      break;
    case 2:
#include "u1/u1-2ch-offdiag-UP.dat"
#include "u1/u1-2ch-offdiag-DO.dat"
#include "u1/u1-2ch-diag.dat"
#include "u1/u1-2ch-spinz.dat" // P.globalB
#include "u1/u1-2ch-spinx.dat" // P.globalBx
      break;
    case 3:
#include "u1/u1-3ch-offdiag-UP.dat"
#include "u1/u1-3ch-offdiag-DO.dat"
#include "u1/u1-3ch-diag.dat"
#include "u1/u1-3ch-spinz.dat" // P.globalB
#include "u1/u1-3ch-spinx.dat" // P.globalBx
      break;
    default: my_assert_not_reached();
  }
}

void SymmetryU1::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  if (P.pol2x2) {
    make_matrix_pol2x2(h, step, qq, I, In, opch, coef);
  } else if (P.polarized) {
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  } else {
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
  }
}

#include "nrg-recalc-U1.cc"
