// *** WARNING!!! Modify nrg-recalc-QST.cc.m4, not nrg-recalc-QST.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (Q,S,T) subspaces

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryQST::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double T   = t1; // trick!
    double S   = (ss1 - 1.) / 2.;
    Invar Ip;

    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,1) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ss1 + 1, t1 - 1);
    RECALC_TAB("qst/qst-doubletp-1.dat", Invar(1, 2, 1));

    Ip = Invar(q1 - 1, ss1 - 1, t1 - 1);
    RECALC_TAB("qst/qst-doubletm-1.dat", Invar(1, 2, 1));

    Ip = Invar(q1 - 1, ss1 + 1, t1);
    RECALC_TAB("qst/qst-doubletp0.dat", Invar(1, 2, 1));

    Ip = Invar(q1 - 1, ss1 - 1, t1);
    RECALC_TAB("qst/qst-doubletm0.dat", Invar(1, 2, 1));

    Ip = Invar(q1 - 1, ss1 + 1, t1 + 1);
    RECALC_TAB("qst/qst-doubletp+1.dat", Invar(1, 2, 1));

    Ip = Invar(q1 - 1, ss1 - 1, t1 + 1);
    RECALC_TAB("qst/qst-doubletm+1.dat", Invar(1, 2, 1));
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
Opch SymmetryQST::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp  = Ip.get("Q");
    Sspin ssp  = Ip.get("SS");
    Tangmom tp = Ip.get("T");
    double T   = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    I1 = Invar(qp + 1, ssp + 1, tp + 1);
    RECALC_F_TAB("qst/qst-spinup+1.dat", 0);

    I1 = Invar(qp + 1, ssp + 1, tp);
    RECALC_F_TAB("qst/qst-spinup0.dat", 0);

    I1 = Invar(qp + 1, ssp + 1, tp - 1);
    RECALC_F_TAB("qst/qst-spinup-1.dat", 0);

    I1 = Invar(qp + 1, ssp - 1, tp + 1);
    RECALC_F_TAB("qst/qst-spindo+1.dat", 0);

    I1 = Invar(qp + 1, ssp - 1, tp);
    RECALC_F_TAB("qst/qst-spindo0.dat", 0);

    I1 = Invar(qp + 1, ssp - 1, tp - 1);
    RECALC_F_TAB("qst/qst-spindo-1.dat", 0);
  }
  return opch;
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetryQST::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    Ip = Invar(q1, ss1, t1);
    RECALC_TAB("qst/qst-triplets.dat", Invar(0, 3, 0));

    Ip = Invar(q1, ss1 + 2, t1);
    RECALC_TAB("qst/qst-tripletp.dat", Invar(0, 3, 0));

    Ip = Invar(q1, ss1 - 2, t1);
    RECALC_TAB("qst/qst-tripletm.dat", Invar(0, 3, 0));
  }
  return cnew;
}

// Recalculate matrix elements of an orbital triplet tenzor operator
MatrixElements SymmetryQST::recalc_orb_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    // 0 = chargeless
    // 1 = spin singlet (deg=2S+1=1)
    // 1 = spin triplet (T=1)

    Ip = Invar(q1, ss1, t1);
    RECALC_TAB("qst/qst-orb-triplets.dat", Invar(0, 1, 1));

    Ip = Invar(q1, ss1, t1 + 1);
    RECALC_TAB("qst/qst-orb-tripletp.dat", Invar(0, 1, 1));

    Ip = Invar(q1, ss1, t1 - 1);
    RECALC_TAB("qst/qst-orb-tripletm.dat", Invar(0, 1, 1));
  }
  return cnew;
}
