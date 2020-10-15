// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Dec 2007
// This file pertains to (S) subspaces

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetrySPSU2::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-doubletp.dat", Invar(2))',
          `RECALC_TAB("spsu2/spsu2-2ch-doubletp.dat", Invar(2))',
          `RECALC_TAB("spsu2/spsu2-3ch-doubletp.dat", Invar(2))');

    Ip = Invar(ss1-1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-doubletm.dat", Invar(2))',
          `RECALC_TAB("spsu2/spsu2-2ch-doubletm.dat", Invar(2))',
          `RECALC_TAB("spsu2/spsu2-3ch-doubletm.dat", Invar(2))');
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
      RECALC_TAB("spsu2/spsu2-1ch-doubletp.dat", Invar(2));

      Ip = Invar(ss1 - 1);
      RECALC_TAB("spsu2/spsu2-1ch-doubletm.dat", Invar(2));
    }
  }
  return cnew;
}

// Driver routine for recalc_f()
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO
Opch SymmetrySPSU2::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  my_assert(!P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    ONE23(`RECALC_F_TAB("spsu2/spsu2-1ch-spinupa.dat", 0)',

          `RECALC_F_TAB("spsu2/spsu2-2ch-spinupa.dat", 0);
	         RECALC_F_TAB("spsu2/spsu2-2ch-spinupb.dat", 1)',
           
          `RECALC_F_TAB("spsu2/spsu2-3ch-spinupa.dat", 0);
       	   RECALC_F_TAB("spsu2/spsu2-3ch-spinupb.dat", 1);
	         RECALC_F_TAB("spsu2/spsu2-3ch-spinupc.dat", 2)');

    I1 = Invar(ssp-1);
    ONE23(`RECALC_F_TAB("spsu2/spsu2-1ch-spindowna.dat", 0)',

          `RECALC_F_TAB("spsu2/spsu2-2ch-spindowna.dat", 0);
           RECALC_F_TAB("spsu2/spsu2-2ch-spindownb.dat", 1)',

          `RECALC_F_TAB("spsu2/spsu2-3ch-spindowna.dat", 0);
           RECALC_F_TAB("spsu2/spsu2-3ch-spindownb.dat", 1);
           RECALC_F_TAB("spsu2/spsu2-3ch-spindownc.dat", 2)');

    // Note: for 3ch cases, the lengths in the three channels are not the same!
    // The same thing occurs for all SU(2)_spin cases, for instance for symtype=QS.
    // RZ, oct 2015
  }
  return opch;
}

// Driver routine for recalc_f()
OpchChannel SymmetrySPSU2::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  my_assert(P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    RECALC_F_TAB("spsu2/spsu2-1ch-spinupa.dat", M);

    I1 = Invar(ssp - 1);
    RECALC_F_TAB("spsu2/spsu2-1ch-spindowna.dat", M);
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetrySPSU2::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-triplets.dat", Invar(3))',
          `RECALC_TAB("spsu2/spsu2-2ch-triplets.dat", Invar(3))',
          `RECALC_TAB("spsu2/spsu2-3ch-triplets.dat", Invar(3))');

    Ip = Invar(ss1+2);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-tripletp.dat", Invar(3))',
          `RECALC_TAB("spsu2/spsu2-2ch-tripletp.dat", Invar(3))',
          `RECALC_TAB("spsu2/spsu2-3ch-tripletp.dat", Invar(3))');

    Ip = Invar(ss1-2);
    ONETWO(`RECALC_TAB("spsu2/spsu2-1ch-tripletm.dat", Invar(3))',
           `RECALC_TAB("spsu2/spsu2-2ch-tripletm.dat", Invar(3))',
           `RECALC_TAB("spsu2/spsu2-3ch-tripletm.dat", Invar(3))');
    }
  } else my_assert_not_reached();
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

void SymmetrySPSU2::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  // NOTE: none of these are implemented for substeps==true.

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Qtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-qdiff.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iztot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Ixtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iptot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Imtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
