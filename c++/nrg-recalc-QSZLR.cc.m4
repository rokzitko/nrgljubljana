// *** WARNING!!! Modify nrg-recalc-QSZLR.cc.m4, not nrg-recalc-QSZLR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006
// This file pertains to (Q,Sz,P) subspaces

include(recalc-macros.m4)

template<typename SC>
Opch_tmpl<SC> SymmetryQSZLR_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp   = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels.

    // ****** CASE I: SAME PARITY ******

    Invar I1;
    I1 = Invar(qp + 1, sszp + 1, pp);
    RECALC_F_TAB("qszlr/qszlr-2ch-spinupa.dat", 0);
    RECALC_F_TAB("qszlr/qszlr-2ch-spinupb.dat", 1);

    I1 = Invar(qp + 1, sszp - 1, pp);
    RECALC_F_TAB("qszlr/qszlr-2ch-spindowna.dat", 0);
    RECALC_F_TAB("qszlr/qszlr-2ch-spindownb.dat", 1);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, sszp + 1, -pp);
    RECALC_F_TAB("qszlr/qszlr-2ch-spinupdiffa.dat", 0);
    RECALC_F_TAB("qszlr/qszlr-2ch-spinupdiffb.dat", 1);

    I1 = Invar(qp + 1, sszp - 1, -pp);
    RECALC_F_TAB("qszlr/qszlr-2ch-spindowndiffa.dat", 0);
    RECALC_F_TAB("qszlr/qszlr-2ch-spindowndiffb.dat", 1);
  }
  return opch;
}
