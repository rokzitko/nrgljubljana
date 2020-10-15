// *** WARNING!!! Modify nrg-recalc-ISOSZLR.cc.m4, not nrg-recalc-ISOSZLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2009
// This file pertains to (I,Sz,P) subspaces

include(recalc-macros.m4)

// Driver routine for recalc_f()
Opch SymmetryISOSZLR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    Ispin iip   = Ip.get("II");
    SZspin sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");

    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip + 1, sszp + 1, pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isoupa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isoupb.dat", 1);

    I1 = Invar(iip + 1, sszp - 1, pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isoupa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isoupb.dat", 1);

    I1 = Invar(iip - 1, sszp + 1, pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isodowna.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isodownb.dat", 1);

    I1 = Invar(iip - 1, sszp - 1, pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isodowna.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isodownb.dat", 1);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip + 1, sszp + 1, -pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isoupdiffa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isoupdiffb.dat", 1);

    I1 = Invar(iip + 1, sszp - 1, -pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isoupdiffa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isoupdiffb.dat", 1);

    I1 = Invar(iip - 1, sszp + 1, -pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isodowndiffa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spinup-isodowndiffb.dat", 1);

    I1 = Invar(iip - 1, sszp - 1, -pp);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isodowndiffa.dat", 0);
    RECALC_F_TAB("isoszlr/isoszlr-2ch-spindown-isodowndiffb.dat", 1);
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
MatrixElements SymmetryISOSZLR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1   = I1.get("II");
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1 - 1, ssz1 + 1, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-doubletmp.dat", Invar(2, -1, +1));

    Ip = Invar(ii1 - 1, ssz1 - 1, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-doubletmm.dat", Invar(2, +1, +1));

    Ip = Invar(ii1 + 1, ssz1 + 1, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-doubletpp.dat", Invar(2, -1, +1));

    Ip = Invar(ii1 + 1, ssz1 - 1, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-doubletpm.dat", Invar(2, +1, +1));
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
MatrixElements SymmetryISOSZLR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1   = I1.get("II");
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ssz1, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-triplets.dat", Invar(1, 0, +1));

    Ip = Invar(ii1, ssz1 + 2, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-tripletp.dat", Invar(1, -2, +1));

    Ip = Invar(ii1, ssz1 - 2, p1);
    RECALC_TAB("isoszlr/isoszlr-2ch-tripletm.dat", Invar(1, +2, +1));
  }
  return cnew;
}
