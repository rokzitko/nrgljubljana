// *** WARNING!!! Modify nrg-recalc-QSZTZ.cc.m4, not nrg-recalc-QSZTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace QSZTZ {
#include "qsztz/qsztz-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryQSZTZ::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1   = I1.get("Q");
    Sspin ssz1  = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp-1.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, -1, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, -1, +1));
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm-1.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, +1, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, +1, +1));
  }
};

    Ip = Invar(q1 - 1, ssz1 + 1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp0.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, -1, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, -1, 0));
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm0.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, +1, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, +1, 0));
  }
};

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp+1.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, -1, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, -1, -1));
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm+1.dat" << ", len=" << QSZTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, +1, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_D_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_D_3CH, Invar(1, +1, -1));
  }
};
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
Opch SymmetryQSZTZ::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp   = Ip.get("Q");
    Sspin sszp  = Ip.get("SZ");
    Tangmom tzp = Ip.get("TZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1, tzp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup+1.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, sszp + 1, tzp);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup0.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, sszp + 1, tzp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup-1.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo+1.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo+1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo0.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo0.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo-1.dat" << ", ch=" << 0 << ", len=" << QSZTZ::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo-1.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_I_3CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, QSZTZ::LENGTH_I_3CH);
  }
};
  }
  return opch;
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetryQSZTZ::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1   = I1.get("Q");
    Sspin ssz1  = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ssz1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-triplets.dat" << ", len=" << QSZTZ::LENGTH_T0_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_T0_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_T0_3CH, Invar(0, 3, 0));
  }
};

    Ip = Invar(q1, ssz1 + 2, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-tripletp.dat" << ", len=" << QSZTZ::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_Tpm_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
  }
};

    Ip = Invar(q1, ssz1 - 2, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-tripletm.dat" << ", len=" << QSZTZ::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "qsztz/qsztz-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZTZ::LENGTH_Tpm_3CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QSZTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
  }
};
  }
  return cnew;
}
