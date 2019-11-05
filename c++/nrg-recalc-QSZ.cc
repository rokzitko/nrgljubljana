// *** WARNING!!! Modify nrg-recalc-QSZ.cc.m4, not nrg-recalc-QSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, 2006-2017
// This file pertains to (Q,SZ) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace QSZ {
#include "qsz/qsz-1ch-def.dat"
#include "qsz/qsz-2ch-def.dat"
#include "qsz/qsz-3ch-def.dat"
}

// NOTE: p is ket (right side), 1 is bra (left side). OP is sandwiched
// in between. Thus Q[p] + Q[op] = Q[1].

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSZ::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
 if (!substeps) {
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    // In the case of (Q,S_z) basis, spin up and spin down are not
    // equivalent. The distinction appears in this recalculation code, but
    // also during the evaluation of the spectral densities, where two
    // spin-diagonal spectral densities (and also spin-flip spectral
    // density) can be defined.

    Ip = Invar(q1-1, ssz1+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletp.dat" << ", len=" << QSZ::LENGTH_D_1CH << ", Iop=" << Invar(1, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_1CH, Invar(1, -1));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-doubletp.dat" << ", len=" << QSZ::LENGTH_D_2CH << ", Iop=" << Invar(1, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-2ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_2CH, Invar(1, -1));
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-doubletp.dat" << ", len=" << QSZ::LENGTH_D_3CH << ", Iop=" << Invar(1, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-3ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_3CH, Invar(1, -1));
 }
}
}
break;
default:
my_assert_not_reached();
};
      
    Ip = Invar(q1-1, ssz1-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletm.dat" << ", len=" << QSZ::LENGTH_D_1CH << ", Iop=" << Invar(1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_1CH, Invar(1, +1));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-doubletm.dat" << ", len=" << QSZ::LENGTH_D_2CH << ", Iop=" << Invar(1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-2ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_2CH, Invar(1, +1));
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-doubletm.dat" << ", len=" << QSZ::LENGTH_D_3CH << ", Iop=" << Invar(1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-3ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_3CH, Invar(1, +1));
 }
}
}
break;
default:
my_assert_not_reached();
};
  } // loop
 } else { // substeps
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q1-1, ssz1+1);
    {
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletp.dat" << ", len=" << QSZ::LENGTH_D_1CH << ", Iop=" << Invar(1, -1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-doubletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_1CH, Invar(1, -1));
 }
};

    Ip = Invar(q1-1, ssz1-1);
    {
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletm.dat" << ", len=" << QSZ::LENGTH_D_1CH << ", Iop=" << Invar(1, +1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-doubletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_D_1CH, Invar(1, +1));
 }
};
  } // loop
 } // if
}

// Driver routine for recalc_f()
void SymmetryQSZ::recalc_irreduc(const DiagInfo &diag)
{
  my_assert(!substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels. 

    I1 = Invar(qp+1, sszp+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-1ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-2ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_2CH);
 }
};
      	   {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QSZ::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-2ch-spinupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSZ::LENGTH_I_2CH);
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
};
           {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spinupb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
};
      	   {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupc.dat" << ", ch=" << 2 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spinupc.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
}
}
break;
default:
my_assert_not_reached();
};

    I1 = Invar(qp+1, sszp-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-1ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-2ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_2CH);
 }
};
	   {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QSZ::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-2ch-spindownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSZ::LENGTH_I_2CH);
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
};
           {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spindownb.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
};
      	   {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindownc.dat" << ", ch=" << 2 << ", len=" << QSZ::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-3ch-spindownc.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, QSZ::LENGTH_I_3CH);
 }
}
}
break;
default:
my_assert_not_reached();
};
  } // loop
}

void SymmetryQSZ::recalc_irreduc_substeps(const DiagInfo &diag, int M)
{
  my_assert(substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(qp+1, sszp+1);
    {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spinupa.dat" << ", ch=" << M << ", len=" << QSZ::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-1ch-spinupa.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_1CH);
        recalc_f(diag, a.opch[M][0], Ip, I1, recalc_table, QSZ::LENGTH_I_1CH);
 }
};

    I1 = Invar(qp+1, sszp-1);
    {
 nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spindowna.dat" << ", ch=" << M << ", len=" << QSZ::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "qsz/qsz-1ch-spindowna.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_I_1CH);
        recalc_f(diag, a.opch[M][0], Ip, I1, recalc_table, QSZ::LENGTH_I_1CH);
 }
};
  } // loop
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSZ::recalc_triplet(DiagInfo &diag,
                                 MatrixElements &cold,
				 MatrixElements &cnew)
{
 if (!substeps) {
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q1, ssz1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-triplets.dat" << ", len=" << QSZ::LENGTH_T0_1CH << ", Iop=" << Invar(0, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_T0_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_T0_1CH, Invar(0, 0));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-triplets.dat" << ", len=" << QSZ::LENGTH_T0_2CH << ", Iop=" << Invar(0, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-2ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_T0_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_T0_2CH, Invar(0, 0));
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-triplets.dat" << ", len=" << QSZ::LENGTH_T0_3CH << ", Iop=" << Invar(0, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-3ch-triplets.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_T0_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_T0_3CH, Invar(0, 0));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(q1, ssz1+2);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-tripletp.dat" << ", len=" << QSZ::LENGTH_Tpm_1CH << ", Iop=" << Invar(0, -2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_Tpm_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_Tpm_1CH, Invar(0, -2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-tripletp.dat" << ", len=" << QSZ::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, -2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-2ch-tripletp.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_Tpm_2CH, Invar(0, -2));
 }
}
}
break;
default:
my_assert_not_reached();
};

    Ip = Invar(q1, ssz1-2);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-tripletm.dat" << ", len=" << QSZ::LENGTH_Tpm_1CH << ", Iop=" << Invar(0, +2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-1ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_Tpm_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_Tpm_1CH, Invar(0, +2));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-tripletm.dat" << ", len=" << QSZ::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, +2) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "qsz/qsz-2ch-tripletm.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSZ::LENGTH_Tpm_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, QSZ::LENGTH_Tpm_2CH, Invar(0, +2));
 }
}
}
break;
default:
my_assert_not_reached();
};
  } // loop 
 } else { // substeps
  my_error("Not implemented.");
 }
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetryQSZ::recalc_global(DiagInfo &diag,
				string name,
				MatrixElements &cnew)
{
   if (name == "SZtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "qsz/qsz-1ch-spinz.dat"
	 break;
	case 2:
#include "qsz/qsz-2ch-spinz.dat"
	 break;
	case 3:
#include "qsz/qsz-3ch-spinz.dat"
	 break;
        default:
	 my_assert_not_reached();
       }
     } // LOOP
   }
}
