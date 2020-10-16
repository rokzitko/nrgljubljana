// *** WARNING!!! Modify nrg-recalc-SPSU2LR.cc.m4, not nrg-recalc-SPSU2LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Sep 2015
// This file pertains to (S,P) subspaces

include(recalc-macros.m4)

MatrixElements SymmetrySPSU2LR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1 + 1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-doubletp.dat", Invar(-1, 1));

    Ip = Invar(ss1 - 1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-doubletm.dat", Invar(+1, 1));
  }
  return cnew;
}

Opch SymmetrySPSU2LR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    int pp    = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(ssp + 1, pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupa.dat", 0);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupb.dat", 1);

    I1 = Invar(ssp - 1, pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowna.dat", 0);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindownb.dat", 1);

    // CASE II: OPPOSITE PARITY

    I1 = Invar(ssp + 1, -pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupdiffa.dat", 0);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupdiffb.dat", 1);

    I1 = Invar(ssp - 1, -pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowndiffa.dat", 0);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowndiffb.dat", 1);
  }
  return opch;
}

MatrixElements SymmetrySPSU2LR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-triplets.dat", Invar(0, 1));

    Ip = Invar(ss1 + 2, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-tripletp.dat", Invar(-2, 1));

    Ip = Invar(ss1 - 2, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-tripletm.dat", Invar(+2, 1));
  }
  return cnew;
}
