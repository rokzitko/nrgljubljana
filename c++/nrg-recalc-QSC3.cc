// *** WARNING!!! Modify nrg-recalc-QSC3.cc.m4, not nrg-recalc-QSC3.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (Q,S,P) subspaces, P modulo 3

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace QSC3 {
#include "qsc3/qsc3-def.dat"
}

#define xRECALC_F_TAB(a,b,c) 0;

// Driver routine for recalc_f()
void SymmetryQSC3::recalc_irreduc(const DiagInfo &diag)
{
#ifdef NRG_COMPLEX
  // CONVENTION: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    int p = Ip.get("P");
    
    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(qp+1, ssp+1, (p+0)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup0-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup0-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup0-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(qp+1, ssp-1, (p+0)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown0-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown0-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown0-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};

    I1 = Invar(qp+1, ssp+1, (p+1)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup1-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup1-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup1-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(qp+1, ssp-1, (p+1)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown1-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown1-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown1-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};

    I1 = Invar(qp+1, ssp+1, (p+2)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup2-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup2-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spinup2-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(qp+1, ssp-1, (p+2)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown2-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown2-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsc3/qsc3-spindown2-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSC3::LENGTH_I_3CH);
 }
};
#undef Power
#undef sqrt
}
#endif
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSC3::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
#ifdef NRG_COMPLEX
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1-1, ss1+1, p1);
//    xRECALC_TAB("qsc3/qsc3-doubletp.dat", QSC3::LENGTH_D_3CH, 
//    					     Invar(1, 2, 0));

    Ip = Invar(q1-1, ss1-1, p1);
//    xRECALC_TAB("qsc3/qsc3-doubletm.dat", QSC3::LENGTH_D_3CH, 
//    					     Invar(1, 2, 0));
  }
#endif
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSC3::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
#ifdef NRG_COMPLEX
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
//    xRECALC_TAB("qsc3/qsc3-triplets.dat", QSC3::LENGTH_T0_3CH,
//    					     Invar(0, 3, 0));

    Ip = Invar(q1, ss1+2, p1);
//    xRECALC_TAB("qsc3/qsc3-tripletp.dat", QSC3::LENGTH_Tpm_3CH,
//    					     Invar(0, 3, 0));

    Ip = Invar(q1, ss1-2, p1);
//    xRECALC_TAB("qsc3/qsc3-tripletm.dat", QSC3::LENGTH_Tpm_3CH, 
//    					     Invar(0, 3, 0));
  }
#endif
}
