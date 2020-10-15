// *** WARNING!!! Modify nrg-recalc-QJ.cc.m4, not nrg-recalc-QJ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,J) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace QJ {
#include "qj/qj-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryQJ::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-doubletp.dat" << ", len=" << QJ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_D_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_D_3CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-doubletm.dat" << ", len=" << QJ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_D_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_D_3CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
  }
  return cnew;
}

#undef If
#define If(cond, a, b) (cond ? a : b)

// Recalculate matrix elements of a quadruplet tensor operator
MatrixElements SymmetryQJ::recalc_quadruplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
    //    double J = (jj1-1.0)/2.0;
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 3);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad1.dat" << ", len=" << QJ::LENGTH_Q1_3CH << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-quad1.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q1_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_Q1_3CH, Invar(1, 4));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1 - 1, jj1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad2.dat" << ", len=" << QJ::LENGTH_Q2_3CH << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-quad2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q2_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_Q2_3CH, Invar(1, 4));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad3.dat" << ", len=" << QJ::LENGTH_Q2_3CH << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-quad3.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q2_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_Q2_3CH, Invar(1, 4));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 3);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad4.dat" << ", len=" << QJ::LENGTH_Q1_3CH << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qj/qj-quad4.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q1_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QJ::LENGTH_Q1_3CH, Invar(1, 4));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
  }
  return cnew;
}

Opch SymmetryQJ::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin jjp = Ip.get("JJ");
    double j  = J(jjp);
    Invar I1;

    I1 = Invar(qp + 1, jjp + 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz3_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jz3_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j3_2-jz3_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jz3_2);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j3_2_jz3_2);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};

    I1 = Invar(qp + 1, jjp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz1_2.dat" << ", ch=" << 0 << ", len=" << QJ::LENGTH_I_3CH_j1_2_jz1_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j1_2-jz1_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j1_2_jz1_2);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j1_2_jz1_2);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz1_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jz1_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j3_2-jz1_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jz1_2);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j3_2_jz1_2);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};

    I1 = Invar(qp + 1, jjp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz-1_2.dat" << ", ch=" << 0 << ", len=" << QJ::LENGTH_I_3CH_j1_2_jzM1_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j1_2-jz-1_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j1_2_jzM1_2);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j1_2_jzM1_2);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-1_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jzM1_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j3_2-jz-1_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jzM1_2);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j3_2_jzM1_2);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};

    I1 = Invar(qp + 1, jjp - 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-3_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jzM3_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qj/qj-spin_j3_2-jz-3_2.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jzM3_2);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QJ::LENGTH_I_3CH_j3_2_jzM3_2);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};
  }
  return opch;
}
