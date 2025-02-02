namespace NRG {

// *** WARNING!!! Modify nrg-recalc-QS.cc.m4, not nrg-recalc-QS.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (Q,S) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetryQS<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
    ONE234(`RECALC_TAB("qs/qs-1ch-doubletp.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-2ch-doubletp.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-3ch-doubletp.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-4ch-doubletp.dat", Invar(1, 2))');

    Ip = Invar(q1-1, ss1-1);
    ONE234(`RECALC_TAB("qs/qs-1ch-doubletm.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-2ch-doubletm.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-3ch-doubletm.dat", Invar(1, 2))',
           `RECALC_TAB("qs/qs-4ch-doubletm.dat", Invar(1, 2))');
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
      RECALC_TAB("qs/qs-1ch-doubletp.dat", Invar(1, 2));

      Ip = Invar(q1 - 1, ss1 - 1);
      RECALC_TAB("qs/qs-1ch-doubletm.dat", Invar(1, 2));
    }
  }
  return cnew;
}

// (QS): Two calls of recalc_f() are necessary (for S+1/2 and S-1/2)
// for each channel.
// See Krishna-Murthy p. 1034, equation (B10).

// <I1 | f^dag | Ip>, hence q1=qp+1

// Driver routine for recalc_f()
template<typename SC>
Opch<SC> SymmetryQS<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    I1 = Invar(qp + 1, ssp + 1);
    ONE234(`RECALC_F_TAB("qs/qs-1ch-spinupa.dat", 0)',

           `RECALC_F_TAB("qs/qs-2ch-spinupa.dat", 0);
	          RECALC_F_TAB("qs/qs-2ch-spinupb.dat", 1)',

           `RECALC_F_TAB("qs/qs-3ch-spinupa.dat", 0);
	          RECALC_F_TAB("qs/qs-3ch-spinupb.dat", 1);
	          RECALC_F_TAB("qs/qs-3ch-spinupc.dat", 2)',

           `RECALC_F_TAB("qs/qs-4ch-spinupa.dat", 0);
	          RECALC_F_TAB("qs/qs-4ch-spinupb.dat", 1);
	          RECALC_F_TAB("qs/qs-4ch-spinupc.dat", 2);
	          RECALC_F_TAB("qs/qs-4ch-spinupd.dat", 3)');

    I1 = Invar(qp+1, ssp-1);
    ONE234(`RECALC_F_TAB("qs/qs-1ch-spindowna.dat", 0)',

           `RECALC_F_TAB("qs/qs-2ch-spindowna.dat", 0);
            RECALC_F_TAB("qs/qs-2ch-spindownb.dat", 1)',

           `RECALC_F_TAB("qs/qs-3ch-spindowna.dat", 0);
	          RECALC_F_TAB("qs/qs-3ch-spindownb.dat", 1);
	          RECALC_F_TAB("qs/qs-3ch-spindownc.dat", 2)',

           `RECALC_F_TAB("qs/qs-4ch-spindowna.dat", 0);
	          RECALC_F_TAB("qs/qs-4ch-spindownb.dat", 1);
	          RECALC_F_TAB("qs/qs-4ch-spindownc.dat", 2);
	          RECALC_F_TAB("qs/qs-4ch-spindownd.dat", 3)');
  }
  return opch;
}
 
// Driver routine for recalc_f() for substeps=true, i.e., chain by chain diagonalisations
template<typename SC>
OpchChannel<SC> SymmetryQS<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const int M) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(qp + 1, ssp + 1);
    RECALC_F_TAB("qs/qs-1ch-spinupa.dat", M);

    I1 = Invar(qp + 1, ssp - 1);
    RECALC_F_TAB("qs/qs-1ch-spindowna.dat", M);
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
template<typename SC>
MatrixElements<SC> SymmetryQS<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1, ss1);
    ONE234(`RECALC_TAB("qs/qs-1ch-triplets.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-2ch-triplets.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-3ch-triplets.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-4ch-triplets.dat", Invar(0, 3))');

    Ip = Invar(q1, ss1+2);
    ONE234(`RECALC_TAB("qs/qs-1ch-tripletp.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-2ch-tripletp.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-3ch-tripletp.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-4ch-tripletp.dat", Invar(0, 3))');

    Ip = Invar(q1, ss1-2);
    ONE234(`RECALC_TAB("qs/qs-1ch-tripletm.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-2ch-tripletm.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-3ch-tripletm.dat", Invar(0, 3))',
           `RECALC_TAB("qs/qs-4ch-tripletm.dat", Invar(0, 3))');
    }
  } else my_assert_not_reached();
  return cnew;
}

#undef QDIFF
#define QDIFF(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef QTOT
#define QTOT(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryQS<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, 
                                        std::string name, MatrixElements<SC> &cnew) const {
  // AAA: m4 macros!! RECALC_GLOBAL("Qdiff", ...)
  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-qdiff.dat"
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
#include "qs/qs-2ch-q1.dat"
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
#include "qs/qs-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

#undef Q2

}
