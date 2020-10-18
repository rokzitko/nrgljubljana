// *** WARNING!!! Modify nrg-recalc-QSLR.cc.m4, not nrg-recalc-QSLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

include(recalc-macros.m4)

template<typename SC>
Opch_tmpl<SC> SymmetryQSLR_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    // IMPORTANT NEW ELEMENT: the parity is important here!!
    int lrp = Ip.get("P");

    // Both parities yield non-zero <Q+1, S+-1/2, P| a^\dag_\nu
    // |Q,S,P'>.  Coefficients *DO* depend on P,P', or more
    // accurately, on whether or not P and P' are the same.
    //
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(qp + 1, ssp + 1, lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spinupa.dat", 0);
    RECALC_F_TAB("qslr/qslr-2ch-spinupb.dat", 1);

    I1 = Invar(qp + 1, ssp - 1, lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spindowna.dat", 0);
    RECALC_F_TAB("qslr/qslr-2ch-spindownb.dat", 1);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, ssp + 1, -lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spinupdiffa.dat", 0);
    RECALC_F_TAB("qslr/qslr-2ch-spinupdiffb.dat", 1);

    I1 = Invar(qp + 1, ssp - 1, -lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spindowndiffa.dat", 0);
    RECALC_F_TAB("qslr/qslr-2ch-spindowndiffb.dat", 1);
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements_tmpl<SC> SymmetryQSLR_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1 - 1, ss1 + 1, p1);
    RECALC_TAB("qslr/qslr-2ch-doubletp.dat", Invar(1, 2, +1));

    Ip = Invar(q1 - 1, ss1 - 1, p1);
    RECALC_TAB("qslr/qslr-2ch-doubletm.dat", Invar(1, 2, +1));
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
template<typename SC>
MatrixElements_tmpl<SC> SymmetryQSLR_tmpl<SC>::recalc_triplet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    RECALC_TAB("qslr/qslr-2ch-triplets.dat", Invar(0, 3, +1));

    Ip = Invar(q1, ss1 + 2, p1);
    RECALC_TAB("qslr/qslr-2ch-tripletp.dat", Invar(0, 3, +1));

    Ip = Invar(q1, ss1 - 2, p1);
    RECALC_TAB("qslr/qslr-2ch-tripletm.dat", Invar(0, 3, +1));
  }
  return cnew;
}
