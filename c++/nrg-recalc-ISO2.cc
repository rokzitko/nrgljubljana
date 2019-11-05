// *** WARNING!!! Modify nrg-recalc-ISO2.cc.m4, not nrg-recalc-ISO2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, May 2008
// This file pertains to (I,S) subspaces
// Version for EVEN number of impurities.

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace ISO2 {
#include "iso2/iso2-1ch-def.dat"
#include "iso2/iso2-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryISO2::recalc_doublet(DiagInfo &diag,
                          MatrixElements &cold,
                          MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1-1, ss1+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-doubletmp.dat" << ", len=" << ISO2::LENGTH_D_1CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-doubletmp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_1CH, Invar(2, 2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-doubletmp.dat" << ", len=" << ISO2::LENGTH_D_2CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-doubletmp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_2CH, Invar(2, 2));
 }
}
}
break;
default:
my_assert_not_reached();
};
      
    Ip = Invar(ii1-1, ss1-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-doubletmm.dat" << ", len=" << ISO2::LENGTH_D_1CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-doubletmm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_1CH, Invar(2, 2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-doubletmm.dat" << ", len=" << ISO2::LENGTH_D_2CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-doubletmm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_2CH, Invar(2, 2));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(ii1+1, ss1+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-doubletpp.dat" << ", len=" << ISO2::LENGTH_D_1CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-doubletpp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_1CH, Invar(2, 2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-doubletpp.dat" << ", len=" << ISO2::LENGTH_D_2CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-doubletpp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_2CH, Invar(2, 2));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(ii1+1, ss1-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-doubletpm.dat" << ", len=" << ISO2::LENGTH_D_1CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-doubletpm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_1CH, Invar(2, 2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-doubletpm.dat" << ", len=" << ISO2::LENGTH_D_2CH << ", Iop=" << Invar(2, 2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-doubletpm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_D_2CH, Invar(2, 2));
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}


// (ISO): Four calls of recalc_f() are necessary for each channel.

// Driver routine for recalc_f()
void SymmetryISO2::recalc_irreduc(const DiagInfo &diag)
{
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels. 

    Ispin iip = Ip.get("II");
    Sspin ssp = Ip.get("SS");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    I1 = Invar(iip+1, ssp+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-1ch-spinup-isoupa.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-1ch-spinup-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spinup-isoupa.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spinup-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
};
     	    {
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spinup-isoupb.dat" << ", ch=" << 1 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spinup-isoupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
}
}
break;
default:
my_assert_not_reached();
};
    
    I1 = Invar(iip+1, ssp-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-1ch-spindown-isoupa.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-1ch-spindown-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spindown-isoupa.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spindown-isoupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
};
	    {
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spindown-isoupb.dat" << ", ch=" << 1 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spindown-isoupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
}
}
break;
default:
my_assert_not_reached();
};

    I1 = Invar(iip-1, ssp+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-1ch-spinup-isodowna.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-1ch-spinup-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spinup-isodowna.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spinup-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
};
	    {
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spinup-isodownb.dat" << ", ch=" << 1 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spinup-isodownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
}
}
break;
default:
my_assert_not_reached();
};

    I1 = Invar(iip-1, ssp-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-1ch-spindown-isodowna.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-1ch-spindown-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spindown-isodowna.dat" << ", ch=" << 0 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spindown-isodowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
};
	    {
 nrglog('f', "RECALC_F(fn=" << "iso2/iso2-2ch-spindown-isodownb.dat" << ", ch=" << 1 << ", len=" << ISO2::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "iso2/iso2-2ch-spindown-isodownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, ISO2::LENGTH_I_2CH);
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryISO2::recalc_triplet(DiagInfo &diag,
                          MatrixElements &cold,
                          MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1, ss1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-triplets.dat" << ", len=" << ISO2::LENGTH_T0_1CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_T0_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_T0_1CH, Invar(1, 3));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-triplets.dat" << ", len=" << ISO2::LENGTH_T0_2CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_T0_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_T0_2CH, Invar(1, 3));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(ii1, ss1+2);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-tripletp.dat" << ", len=" << ISO2::LENGTH_Tpm_1CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_Tpm_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_Tpm_1CH, Invar(1, 3));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-tripletp.dat" << ", len=" << ISO2::LENGTH_Tpm_2CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_Tpm_2CH, Invar(1, 3));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(ii1, ss1-2);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-1ch-tripletm.dat" << ", len=" << ISO2::LENGTH_Tpm_1CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-1ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_Tpm_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_Tpm_1CH, Invar(1, 3));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "iso2/iso2-2ch-tripletm.dat" << ", len=" << ISO2::LENGTH_Tpm_2CH << ", Iop=" << Invar(1, 3) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "iso2/iso2-2ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ISO2::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ISO2::LENGTH_Tpm_2CH, Invar(1, 3));
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}
