// *** WARNING!!! Modify nrg-recalc-QST.cc.m4, not nrg-recalc-QST.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (Q,S,T) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





   namespace QST {
#include "qst/qst-def.dat"
}

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
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp-1.dat" << ", len=" << QST::LENGTH_D_3CH_a << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletp-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_a);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm-1.dat" << ", len=" << QST::LENGTH_D_3CH_a << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletm-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_a);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
  }
};

    Ip = Invar(q1 - 1, ss1 + 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp0.dat" << ", len=" << QST::LENGTH_D_3CH_b << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletp0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_b);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_b, Invar(1, 2, 1));
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm0.dat" << ", len=" << QST::LENGTH_D_3CH_b << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletm0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_b);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_b, Invar(1, 2, 1));
  }
};

    Ip = Invar(q1 - 1, ss1 + 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp+1.dat" << ", len=" << QST::LENGTH_D_3CH_a << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletp+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_a);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm+1.dat" << ", len=" << QST::LENGTH_D_3CH_a << ", Iop=" << Invar(1, 2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-doubletm+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_D_3CH_a);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
  }
};
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
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup+1.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_0 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spinup+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_0);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_0);
  }
};

    I1 = Invar(qp + 1, ssp + 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup0.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_1 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spinup0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_1);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_1);
  }
};

    I1 = Invar(qp + 1, ssp + 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup-1.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_2 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spinup-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_2);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_2);
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo+1.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_0 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spindo+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_0);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_0);
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo0.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_1 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spindo0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_1);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_1);
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo-1.dat" << ", ch=" << 0 << ", len=" << QST::LENGTH_I_3CH_2 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qst/qst-spindo-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_I_3CH_2);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QST::LENGTH_I_3CH_2);
  }
};
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
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-triplets.dat" << ", len=" << QST::LENGTH_T0_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_T0_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_T0_3CH, Invar(0, 3, 0));
  }
};

    Ip = Invar(q1, ss1 + 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-tripletp.dat" << ", len=" << QST::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_Tpm_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_Tpm_3CH, Invar(0, 3, 0));
  }
};

    Ip = Invar(q1, ss1 - 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-tripletm.dat" << ", len=" << QST::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_Tpm_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_Tpm_3CH, Invar(0, 3, 0));
  }
};
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tenzor operator
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
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-triplets.dat" << ", len=" << QST::LENGTH_OT0_3CH << ", Iop=" << Invar(0, 1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-orb-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_OT0_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_OT0_3CH, Invar(0, 1, 1));
  }
};

    Ip = Invar(q1, ss1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-tripletp.dat" << ", len=" << QST::LENGTH_OTpm_3CH << ", Iop=" << Invar(0, 1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-orb-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_OTpm_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_OTpm_3CH, Invar(0, 1, 1));
  }
};

    Ip = Invar(q1, ss1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-tripletm.dat" << ", len=" << QST::LENGTH_OTpm_3CH << ", Iop=" << Invar(0, 1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qst/qst-orb-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QST::LENGTH_OTpm_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, QST::LENGTH_OTpm_3CH, Invar(0, 1, 1));
  }
};
  }
  return cnew;
}
