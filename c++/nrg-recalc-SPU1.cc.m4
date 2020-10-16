// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2008.
// This file pertains to (S) subspaces

include(recalc-macros.m4)

MatrixElements SymmetrySPU1::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1 + 1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-doubletp.dat", Invar(-1))',
           `RECALC_TAB("spu1/spu1-2ch-doubletp.dat", Invar(-1))');

     Ip = Invar(ssz1-1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-doubletm.dat", Invar(+1))',
            `RECALC_TAB("spu1/spu1-2ch-doubletm.dat", Invar(+1))');
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1 + 1);
      RECALC_TAB("spu1/spu1-1ch-doubletp.dat", Invar(-1));

      Ip = Invar(ssz1 - 1);
      RECALC_TAB("spu1/spu1-1ch-doubletm.dat", Invar(+1));
    }
  }
  return cnew;
}

Opch SymmetrySPU1::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  my_assert(!P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
     ONETWO(`RECALC_F_TAB("spu1/spu1-1ch-spinupa.dat", 0)',

            `RECALC_F_TAB("spu1/spu1-2ch-spinupa.dat", 0);
 	           RECALC_F_TAB("spu1/spu1-2ch-spinupb.dat", 1)');

     I1 = Invar(sszp-1);
     ONETWO(`RECALC_F_TAB("spu1/spu1-1ch-spindowna.dat", 0)',

            `RECALC_F_TAB("spu1/spu1-2ch-spindowna.dat", 0);
             RECALC_F_TAB("spu1/spu1-2ch-spindownb.dat", 1)');
  }
  return opch;
}

OpchChannel SymmetrySPU1::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  my_assert(P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
    RECALC_F_TAB("spu1/spu1-1ch-spinupa.dat", M);

    I1 = Invar(sszp - 1);
    RECALC_F_TAB("spu1/spu1-1ch-spindowna.dat", M);
  }
  return opch[M];
}

MatrixElements SymmetrySPU1::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-triplets.dat", Invar(0))',
           `RECALC_TAB("spu1/spu1-2ch-triplets.dat", Invar(0))');

     Ip = Invar(ssz1+2);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-tripletp.dat", Invar(-2))',
           `RECALC_TAB("spu1/spu1-2ch-tripletp.dat", Invar(-2))');

     Ip = Invar(ssz1-2);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-tripletm.dat", Invar(+2))',
            `RECALC_TAB("spu1/spu1-2ch-tripletm.dat", Invar(+2))');
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1);
      RECALC_TAB("spu1/spu1-1ch-triplets.dat", Invar(0));

      Ip = Invar(ssz1 + 2);
      RECALC_TAB("spu1/spu1-1ch-tripletp.dat", Invar(-2));

      Ip = Invar(ssz1 - 2);
      RECALC_TAB("spu1/spu1-1ch-tripletm.dat", Invar(+2));
    }
  }
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

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

// 2-channel only

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1UP
#define Q1UP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2UP
#define Q2UP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1DO
#define Q1DO(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2DO
#define Q2DO(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

void SymmetrySPU1::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1/spu1-1ch-Qtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Qtot.dat"
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
#include "spu1/spu1-1ch-Iztot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Iztot.dat"
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
#include "spu1/spu1-1ch-Ixtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Ixtot.dat"
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
#include "spu1/spu1-1ch-Iptot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Iptot.dat"
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
#include "spu1/spu1-1ch-Imtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Imtot.dat"
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
#include "spu1/spu1-2ch-qdiff.dat"
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
#include "spu1/spu1-2ch-q1.dat"
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
#include "spu1/spu1-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1UP") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spu1/spu1-2ch-q1up.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2UP") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spu1/spu1-2ch-q2up.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1DO") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spu1/spu1-2ch-q1do.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2DO") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spu1/spu1-2ch-q2do.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

// XXX: use m4 macros rather than C++ macros
#undef Q1
#undef Q2
#undef Q3
