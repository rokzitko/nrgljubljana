// *** WARNING!!! Modify nrg-recalc-SPSU2LR.cc.m4, not nrg-recalc-SPSU2LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Sep 2015
// This file pertains to (S,P) subspaces

namespace SPSU2LR {
#include "spsu2lr/spsu2lr-2ch-def.dat"
}

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















// Recalculate matrix elements of a doublet tensor operator
void SymmetrySPSU2LR::recalc_doublet(DiagInfo &diag,
                                     MatrixElements &cold,
				     MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1+1, p1);
    {
 nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-doubletp.dat" << ", len=" << SPSU2LR::LENGTH_D_2CH << ", Iop=" << Invar(-1, 1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SPSU2LR::LENGTH_D_2CH, Invar(-1, 1));
 }
};
      
    Ip = Invar(ss1-1, p1);
    {
 nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-doubletm.dat" << ", len=" << SPSU2LR::LENGTH_D_2CH << ", Iop=" << Invar(+1, 1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SPSU2LR::LENGTH_D_2CH, Invar(+1, 1));
 }
};
  }
}

// Driver routine for recalc_f()
void SymmetrySPSU2LR::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Sspin ssp = Ip.get("SS");
    int pp = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(ssp+1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spinupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(ssp-1, pp);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spindownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};

    // CASE II: OPPOSITE PARITY

    I1 = Invar(ssp+1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spinupdiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spinupdiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    
    I1 = Invar(ssp-1, -pp);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spindowndiffa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ", len=" << SPSU2LR::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-spindowndiffb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2LR::LENGTH_I_2CH);
 }
};
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetrySPSU2LR::recalc_triplet(DiagInfo &diag,
                                     MatrixElements &cold,
                                     MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1, p1);
    {
 nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-triplets.dat" << ", len=" << SPSU2LR::LENGTH_T0_2CH << ", Iop=" << Invar(0, 1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_T0_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SPSU2LR::LENGTH_T0_2CH, Invar(0, 1));
 }
};
      
    Ip = Invar(ss1+2, p1);
    {
 nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-tripletp.dat" << ", len=" << SPSU2LR::LENGTH_Tpm_2CH << ", Iop=" << Invar(-2, 1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SPSU2LR::LENGTH_Tpm_2CH, Invar(-2, 1));
 }
};
      
    Ip = Invar(ss1-2, p1);
    {
 nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-tripletm.dat" << ", len=" << SPSU2LR::LENGTH_Tpm_2CH << ", Iop=" << Invar(+2, 1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "spsu2lr/spsu2lr-2ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2LR::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SPSU2LR::LENGTH_Tpm_2CH, Invar(+2, 1));
 }
};
  }
}
