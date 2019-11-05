// *** WARNING!!! Modify nrg-recalc-QSTZ.cc.m4, not nrg-recalc-QSTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2016
// This file pertains to (Q,S,Tz) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace QSTZ {
#include "qstz/qstz-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSTZ::recalc_doublet(DiagInfo &diag,
                                  MatrixElements &cold,
                                  MatrixElements &cnew)
{
  nrglog('f', "QSTZ::recalc_doublet() called");
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;
    
    nrglog('f', "I1=" << I1);
    
    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet, 
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1-1, ss1+1, tz1-1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletp-1.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletp-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, +1));
 }
};
      
    Ip = Invar(q1-1, ss1-1, tz1-1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletm-1.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletm-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, +1));
 }
};

    Ip = Invar(q1-1, ss1+1, tz1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletp0.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletp0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, 0));
 }
};
      
    Ip = Invar(q1-1, ss1-1, tz1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletm0.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletm0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, 0));
 }
};

    Ip = Invar(q1-1, ss1+1, tz1+1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletp+1.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletp+1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, -1));
 }
};
      
    Ip = Invar(q1-1, ss1-1, tz1+1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-doubletm+1.dat" << ", len=" << QSTZ::LENGTH_D_3CH << ", Iop=" << Invar(1, 2, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-doubletm+1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_D_3CH, Invar(1, 2, -1));
 }
};
 }
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
void SymmetryQSTZ::recalc_irreduc(const DiagInfo &diag)
{
  nrglog('f', "QSTZ::recalc_irreduc() called");
  my_assert(!substeps);

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Tangmom tzp = Ip.get("TZ");
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces 
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    nrglog('f', "qp=" << qp << " ssp=" << ssp << " tzp=" << tzp);
    nrglog('f', "spinup+1");
    I1 = Invar(qp+1, ssp+1, tzp+1);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spinup+1.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spinup+1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};
    
    nrglog('f', "spinup0");
    I1 = Invar(qp+1, ssp+1, tzp);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spinup0.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spinup0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};

    nrglog('f', "spinup-1");
    I1 = Invar(qp+1, ssp+1, tzp-1);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spinup-1.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spinup-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};
    
    nrglog('f', "spindo+1");
    I1 = Invar(qp+1, ssp-1, tzp+1);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spindo+1.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spindo+1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};

    nrglog('f', "spindo0");
    I1 = Invar(qp+1, ssp-1, tzp);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spindo0.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spindo0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};
    
    nrglog('f', "spindo-1");
    I1 = Invar(qp+1, ssp-1, tzp-1);
    {
 nrglog('f', "RECALC_F(fn=" << "qstz/qstz-spindo-1.dat" << ", ch=" << 0 << ", len=" << QSTZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qstz/qstz-spindo-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSTZ::LENGTH_I_3CH);
 }
};
  }
}


// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSTZ::recalc_triplet(DiagInfo &diag,
                                  MatrixElements &cold,
                                  MatrixElements &cnew)
{
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ss1, tz1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-triplets.dat" << ", len=" << QSTZ::LENGTH_T0_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_T0_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_T0_3CH, Invar(0, 3, 0));
 }
};
      
    Ip = Invar(q1, ss1+2, tz1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-tripletp.dat" << ", len=" << QSTZ::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_Tpm_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
 }
};
      
    Ip = Invar(q1, ss1-2, tz1);
    {
 nrglog('f', "RECALC(fn=" << "qstz/qstz-tripletm.dat" << ", len=" << QSTZ::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qstz/qstz-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSTZ::LENGTH_Tpm_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
 }
};
   }
}


