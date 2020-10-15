// *** WARNING!!! Modify nrg-recalc-NONE.cc.m4, not nrg-recalc-NONE.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, April 2010
// This file pertains to the case with no symmetry.

include(recalc-macros.m4)

// Driver routine for recalc_f()
Opch SymmetryNONE::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1 = Invar();

    ONETWO(`RECALC_F_TAB_N("none/none-1ch-a-CR-DO.dat", 0, 0);
            RECALC_F_TAB_N("none/none-1ch-a-CR-UP.dat", 0, 1);',

           `RECALC_F_TAB_N("none/none-2ch-a-CR-DO.dat", 0, 0);
	          RECALC_F_TAB_N("none/none-2ch-b-CR-DO.dat", 1, 0);
            RECALC_F_TAB_N("none/none-2ch-a-CR-UP.dat", 0, 1);
            RECALC_F_TAB_N("none/none-2ch-b-CR-UP.dat", 1, 1)');
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryNONE::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Invar Ip = Invar();

    ONETWO(`RECALC_TAB("none/none-1ch-doublet.dat", Invar())',
           `RECALC_TAB("none/none-2ch-doublet.dat", Invar())');
  }
  return cnew;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

// Isospin operator need an appropriate phase factor (bipartite
// sublattice index)
#define USEISOFACTOR

#if defined(USEISOFACTOR)
#define ISOFACTOR psgn(step.getnn() + 1)
#else
#define ISOFACTOR 1
#endif

#ifdef NRG_COMPLEX
#undef SPINY
#define SPINY(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINY
#define ISOSPINY(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *complex<double>(ISOFACTOR))

#undef Complex
#define Complex(x, y) cmpl(x, y)
#endif // NRG_COMPLEX

#undef CHARGE
#define CHARGE(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

void SymmetryNONE::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-spinz.dat"
          break;
        case 2:
#include "none/none-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

#ifdef NRG_COMPLEX
  if (name == "SYtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-spiny.dat"
          break;
        case 2:
#include "none/none-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-spinx.dat"
          break;
        case 2:
#include "none/none-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Qtot.dat"
          break;
        case 2:
#include "none/none-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iztot.dat"
          break;
        case 2:
#include "none/none-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Ixtot.dat"
          break;
        case 2:
#include "none/none-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

#ifdef NRG_COMPLEX
  if (name == "Iytot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iytot.dat"
          break;
        case 2:
#include "none/none-2ch-Iytot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iptot.dat"
          break;
        case 2:
#include "none/none-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Imtot.dat"
          break;
        case 2:
#include "none/none-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  
  my_assert_not_reached();
}
