// *** WARNING!!! Modify nrg-recalc-ISOSZLR.cc.m4, not nrg-recalc-ISOSZLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2009
// This file pertains to (I,Sz,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace ISOSZLR {
#include "isoszlr/isoszlr-2ch-def.dat"
}

// Driver routine for recalc_f()
void SymmetryISOSZLR::recalc_irreduc(const DiagInfo &diag)
{
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;
    
    Ispin iip = Ip.get("II");
    SZspin sszp = Ip.get("SSZ");
    int pp = Ip.get("P");
       
    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip+1, sszp+1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isoupa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isoupb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isoupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(iip+1, sszp-1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isoupa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isoupb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isoupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};

    I1 = Invar(iip-1, sszp+1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isodowna.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isodownb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isodownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};

    I1 = Invar(iip-1, sszp-1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isodowna.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isodownb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isodownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip+1, sszp+1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isoupdiffa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isoupdiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isoupdiffb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isoupdiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(iip+1, sszp-1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isoupdiffa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isoupdiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isoupdiffb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isoupdiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
      
    I1 = Invar(iip-1, sszp+1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isodowndiffa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isodowndiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spinup-isodowndiffb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spinup-isodowndiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};

    I1 = Invar(iip-1, sszp-1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isodowndiffa.dat" << ", ch=" << 0 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isodowndiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "isoszlr/isoszlr-2ch-spindown-isodowndiffb.dat" << ", ch=" << 1 << ", len=" << ISOSZLR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "isoszlr/isoszlr-2ch-spindown-isodowndiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISOSZLR::LENGTH_I_2CH);
 }
};
  }
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
void SymmetryISOSZLR::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    SZspin ssz1 = I1.get("SSZ");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1-1, ssz1+1, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-doubletmp.dat" << ", len=" << ISOSZLR::LENGTH_D_2CH << ", Iop=" << Invar(2, -1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-doubletmp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_D_2CH, Invar(2, -1, +1));
 }
};

    Ip = Invar(ii1-1, ssz1-1, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-doubletmm.dat" << ", len=" << ISOSZLR::LENGTH_D_2CH << ", Iop=" << Invar(2, +1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-doubletmm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_D_2CH, Invar(2, +1, +1));
 }
};

    Ip = Invar(ii1+1, ssz1+1, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-doubletpp.dat" << ", len=" << ISOSZLR::LENGTH_D_2CH << ", Iop=" << Invar(2, -1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-doubletpp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_D_2CH, Invar(2, -1, +1));
 }
};

    Ip = Invar(ii1+1, ssz1-1, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-doubletpm.dat" << ", len=" << ISOSZLR::LENGTH_D_2CH << ", Iop=" << Invar(2, +1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-doubletpm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_D_2CH, Invar(2, +1, +1));
 }
};
  }
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
void SymmetryISOSZLR::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    SZspin ssz1 = I1.get("SSZ");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ssz1, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-triplets.dat" << ", len=" << ISOSZLR::LENGTH_T0_2CH << ", Iop=" << Invar(1, 0, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_T0_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_T0_2CH, Invar(1, 0, +1));
 }
};

    Ip = Invar(ii1, ssz1+2, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-tripletp.dat" << ", len=" << ISOSZLR::LENGTH_Tpm_2CH << ", Iop=" << Invar(1, -2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_Tpm_2CH, Invar(1, -2, +1));
 }
};

    Ip = Invar(ii1, ssz1-2, p1);
    {
 nrglog('f', "RECALC(fn=" << "isoszlr/isoszlr-2ch-tripletm.dat" << ", len=" << ISOSZLR::LENGTH_Tpm_2CH << ", Iop=" << Invar(1, +2, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "isoszlr/isoszlr-2ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISOSZLR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISOSZLR::LENGTH_Tpm_2CH, Invar(1, +2, +1));
 }
};
  }
}
