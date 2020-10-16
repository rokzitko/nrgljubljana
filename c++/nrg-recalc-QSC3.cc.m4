// *** WARNING!!! Modify nrg-recalc-QSC3.cc.m4, not nrg-recalc-QSC3.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (Q,S,P) subspaces, P modulo 3

include(recalc-macros.m4)

#define xRECALC_F_TAB(a, b, c) 0;

// Driver routine for recalc_f()
Opch SymmetryQSC3::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
#ifdef NRG_COMPLEX

  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    int p     = Ip.get("P");

    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(qp + 1, ssp + 1, (p + 0) % 3);
    RECALC_F_TAB("qsc3/qsc3-spinup0-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spinup0-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spinup0-c.dat", 2);

    I1 = Invar(qp + 1, ssp - 1, (p + 0) % 3);
    RECALC_F_TAB("qsc3/qsc3-spindown0-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spindown0-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spindown0-c.dat", 2);

    I1 = Invar(qp + 1, ssp + 1, (p + 1) % 3);
    RECALC_F_TAB("qsc3/qsc3-spinup1-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spinup1-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spinup1-c.dat", 2);

    I1 = Invar(qp + 1, ssp - 1, (p + 1) % 3);
    RECALC_F_TAB("qsc3/qsc3-spindown1-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spindown1-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spindown1-c.dat", 2);

    I1 = Invar(qp + 1, ssp + 1, (p + 2) % 3);
    RECALC_F_TAB("qsc3/qsc3-spinup2-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spinup2-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spinup2-c.dat", 2);

    I1 = Invar(qp + 1, ssp - 1, (p + 2) % 3);
    RECALC_F_TAB("qsc3/qsc3-spindown2-a.dat", 0);
    RECALC_F_TAB("qsc3/qsc3-spindown2-b.dat", 1);
    RECALC_F_TAB("qsc3/qsc3-spindown2-c.dat", 2);
#undef Power
#undef sqrt
  }
#endif
  return opch;
}
