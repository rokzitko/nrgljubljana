// *** WARNING!!! Modify nrg-recalc-QSLR.cc.m4, not nrg-recalc-QSLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace QSLR {
#include "qslr/qslr-2ch-def.dat"
}

// Driver routine for recalc_f()
void SymmetryQSLR::recalc_irreduc(const DiagInfo &diag)
{
  // CONVENTION: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
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

    I1 = Invar(qp+1, ssp+1, lrp);
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spinupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(qp+1, ssp-1, lrp);
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spindownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp+1, ssp+1, -lrp);
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spinupdiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spinupdiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(qp+1, ssp-1, -lrp);
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spindowndiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ", len=" << QSLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qslr/qslr-2ch-spindowndiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSLR::LENGTH_I_2CH);
 }
};
  }
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
void SymmetryQSLR::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1-1, ss1+1, p1);
    {
 nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletp.dat" << ", len=" << QSLR::LENGTH_D_2CH << ", Iop=" << Invar(1, 2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qslr/qslr-2ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSLR::LENGTH_D_2CH, Invar(1, 2, +1));
 }
};

    Ip = Invar(q1-1, ss1-1, p1);
    {
 nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletm.dat" << ", len=" << QSLR::LENGTH_D_2CH << ", Iop=" << Invar(1, 2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qslr/qslr-2ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSLR::LENGTH_D_2CH, Invar(1, 2, +1));
 }
};
  }
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
void SymmetryQSLR::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    {
 nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-triplets.dat" << ", len=" << QSLR::LENGTH_T0_2CH << ", Iop=" << Invar(0, 3, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qslr/qslr-2ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_T0_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSLR::LENGTH_T0_2CH, Invar(0, 3, +1));
 }
};

    Ip = Invar(q1, ss1+2, p1);
    {
 nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletp.dat" << ", len=" << QSLR::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, 3, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qslr/qslr-2ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSLR::LENGTH_Tpm_2CH, Invar(0, 3, +1));
 }
};

    Ip = Invar(q1, ss1-2, p1);
    {
 nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletm.dat" << ", len=" << QSLR::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, 3, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qslr/qslr-2ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSLR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSLR::LENGTH_Tpm_2CH, Invar(0, 3, +1));
 }
};
  }
}
