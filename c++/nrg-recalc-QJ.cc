// *** WARNING!!! Modify nrg-recalc-QJ.cc.m4, not nrg-recalc-QJ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,J) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace QJ {
#include "qj/qj-def.dat"
}


// Recalculate matrix elements of a doublet tensor operator
void SymmetryQJ::recalc_doublet(DiagInfo &diag,
                                MatrixElements &cold,
                                MatrixElements &cnew)
{
  nrglog('f', "QJ::recalc_doublet() called");
 
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
    Invar Ip;

    Ip = Invar(q1-1, jj1+1);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-doubletp.dat" << ", len=" << QJ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_D_3CH, Invar(1, 2));
 }
};

    Ip = Invar(q1-1, jj1-1);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-doubletm.dat" << ", len=" << QJ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_D_3CH, Invar(1, 2));
 }
};
  }

  nrglog('f', "QJ::recalc_doublet() end");
}

#undef If
#define If(cond, a, b) (cond ? a : b)

// Recalculate matrix elements of a quadruplet tensor operator
void SymmetryQJ::recalc_quadruplet(DiagInfo &diag,
                                   MatrixElements &cold,
                                   MatrixElements &cnew)
{
  nrglog('f', "QJ::recalc_quadruplet() called");
 
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
//    double J = (jj1-1.0)/2.0;
    Invar Ip;

    Ip = Invar(q1-1, jj1+3);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-quad1.dat" << ", len=" << QJ::LENGTH_Q1_3CH << ", Iop=" << Invar(1, 4) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-quad1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q1_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_Q1_3CH, Invar(1, 4));
 }
};

    Ip = Invar(q1-1, jj1+1);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-quad2.dat" << ", len=" << QJ::LENGTH_Q2_3CH << ", Iop=" << Invar(1, 4) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-quad2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q2_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_Q2_3CH, Invar(1, 4));
 }
};

    Ip = Invar(q1-1, jj1-1);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-quad3.dat" << ", len=" << QJ::LENGTH_Q2_3CH << ", Iop=" << Invar(1, 4) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-quad3.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q2_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_Q2_3CH, Invar(1, 4));
 }
};

    Ip = Invar(q1-1, jj1-3);
    {
 nrglog('f', "RECALC(fn=" << "qj/qj-quad4.dat" << ", len=" << QJ::LENGTH_Q1_3CH << ", Iop=" << Invar(1, 4) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qj/qj-quad4.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_Q1_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QJ::LENGTH_Q1_3CH, Invar(1, 4));
 }
};
  }

  nrglog('f', "QJ::recalc_quadruplet() end");
}

void SymmetryQJ::recalc_irreduc(const DiagInfo &diag)
{
  nrglog('f', "SymmetryQJ::recalc_irreduc() start");

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin jjp = Ip.get("JJ");
    double j = J(jjp);
    Invar I1;

    I1 = Invar(qp+1, jjp+3);
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz3_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jz3_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j3_2-jz3_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jz3_2);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j3_2_jz3_2);
 }
};

    I1 = Invar(qp+1, jjp+1);
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz1_2.dat" << ", ch=" << 0 << ", len=" << QJ::LENGTH_I_3CH_j1_2_jz1_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j1_2-jz1_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j1_2_jz1_2);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j1_2_jz1_2);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz1_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jz1_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j3_2-jz1_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jz1_2);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j3_2_jz1_2);
 }
};

    I1 = Invar(qp+1, jjp-1);
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz-1_2.dat" << ", ch=" << 0 << ", len=" << QJ::LENGTH_I_3CH_j1_2_jzM1_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j1_2-jz-1_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j1_2_jzM1_2);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j1_2_jzM1_2);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-1_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jzM1_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j3_2-jz-1_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jzM1_2);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j3_2_jzM1_2);
 }
};

    I1 = Invar(qp+1, jjp-3);
    {
 nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-3_2.dat" << ", ch=" << 1 << ", len=" << QJ::LENGTH_I_3CH_j3_2_jzM3_2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qj/qj-spin_j3_2-jz-3_2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QJ::LENGTH_I_3CH_j3_2_jzM3_2);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QJ::LENGTH_I_3CH_j3_2_jzM3_2);
 }
};
  }
  
  nrglog('f', "SymmetryQJ::recalc_irreduc() end");
}

