namespace NRG {

// *** WARNING!!! Modify nrg-recalc-QSZ.cc.m4, not nrg-recalc-QSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, 2006-2020
// This file pertains to (Q,SZ) subspaces

include(recalc-macros.m4)

// NOTE: p is ket (right side), 1 is bra (left side). OP is sandwiched in between. Thus Q[p] + Q[op] = Q[1].

template<typename SC>
MatrixElements<SC> SymmetryQSZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1   = I1.get("Q");
      int ssz1 = I1.get("SSZ");
      Invar Ip;

      // In the case of (Q,S_z) basis, spin up and spin down are not
      // equivalent. The distinction appears in this recalculation code, but
      // also during the evaluation of the spectral densities, where two
      // spin-diagonal spectral densities (and also spin-flip spectral
      // density) can be defined.

      Ip = Invar(q1 - 1, ssz1 + 1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-doubletp.dat", Invar(1, -1))',
          `RECALC_TAB("qsz/qsz-2ch-doubletp.dat", Invar(1, -1))',
          `RECALC_TAB("qsz/qsz-3ch-doubletp.dat", Invar(1, -1))');

    Ip = Invar(q1-1, ssz1-1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-doubletm.dat", Invar(1, +1))',
          `RECALC_TAB("qsz/qsz-2ch-doubletm.dat", Invar(1, +1))',
          `RECALC_TAB("qsz/qsz-3ch-doubletm.dat", Invar(1, +1))');
    }      // loop
  } else { // substeps
    for(const auto &[I1, eig]: diag) {
      int q1   = I1.get("Q");
      int ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(q1 - 1, ssz1 + 1);
      RECALC_TAB("qsz/qsz-1ch-doubletp.dat", Invar(1, -1));

      Ip = Invar(q1 - 1, ssz1 - 1);
      RECALC_TAB("qsz/qsz-1ch-doubletm.dat", Invar(1, +1));
    } // loop
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryQSZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int sszp = Ip.get("SSZ");
    Invar I1;

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels.

    I1 = Invar(qp + 1, sszp + 1);
    ONE23(`RECALC_F_TAB("qsz/qsz-1ch-spinupa.dat", 0)',

          `RECALC_F_TAB("qsz/qsz-2ch-spinupa.dat", 0);
      	   RECALC_F_TAB("qsz/qsz-2ch-spinupb.dat", 1)',

          `RECALC_F_TAB("qsz/qsz-3ch-spinupa.dat", 0);
           RECALC_F_TAB("qsz/qsz-3ch-spinupb.dat", 1);
      	   RECALC_F_TAB("qsz/qsz-3ch-spinupc.dat", 2)');

    I1 = Invar(qp+1, sszp-1);
    ONE23(`RECALC_F_TAB("qsz/qsz-1ch-spindowna.dat", 0)',

          `RECALC_F_TAB("qsz/qsz-2ch-spindowna.dat", 0);
	         RECALC_F_TAB("qsz/qsz-2ch-spindownb.dat", 1)',

          `RECALC_F_TAB("qsz/qsz-3ch-spindowna.dat", 0);
           RECALC_F_TAB("qsz/qsz-3ch-spindownb.dat", 1);
      	   RECALC_F_TAB("qsz/qsz-3ch-spindownc.dat", 2)');
  } // loop
  return opch;
}

template<typename SC>
OpchChannel<SC> SymmetryQSZ<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const int M) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1);
    RECALC_F_TAB("qsz/qsz-1ch-spinupa.dat", M);

    I1 = Invar(qp + 1, sszp - 1);
    RECALC_F_TAB("qsz/qsz-1ch-spindowna.dat", M);
  } // loop
  return opch[M];
}

template<typename SC>
MatrixElements<SC> SymmetryQSZ<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1   = I1.get("Q");
      int ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(q1, ssz1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-triplets.dat", Invar(0, 0))',
          `RECALC_TAB("qsz/qsz-2ch-triplets.dat", Invar(0, 0))',
          `RECALC_TAB("qsz/qsz-3ch-triplets.dat", Invar(0, 0))');

    Ip = Invar(q1, ssz1+2);
    ONETWO(`RECALC_TAB("qsz/qsz-1ch-tripletp.dat", Invar(0, -2))',
           `RECALC_TAB("qsz/qsz-2ch-tripletp.dat", Invar(0, -2))',
           `RECALC_TAB("qsz/qsz-3ch-tripletp.dat", Invar(0, -2))');

    Ip = Invar(q1, ssz1-2);
    ONETWO(`RECALC_TAB("qsz/qsz-1ch-tripletm.dat", Invar(0, +2))',
           `RECALC_TAB("qsz/qsz-2ch-tripletm.dat", Invar(0, +2))',
           `RECALC_TAB("qsz/qsz-3ch-tripletm.dat", Invar(0, +2))');
    }      // loop
  } else { // substeps
    my_assert_not_reached();
  }
  return cnew;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)
#undef Q1U
#define Q1U(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)
#undef Q1D
#define Q1D(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryQSZ<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-spinz.dat"
          break;
        case 2:
#include "qsz/qsz-2ch-spinz.dat"
          break;
        case 3:
#include "qsz/qsz-3ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1u") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-q1u.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1d") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-q1d.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
