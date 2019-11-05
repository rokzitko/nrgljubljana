// *** WARNING!!! Modify nrg-recalc-ANYJ.cc.m4, not nrg-recalc-ANYJ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007
// This file pertains to (Q,SZ) subspaces for arbitrary spin

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace ANYJ {
const unsigned int LENGTH_D_1CH_N2 = 4;
const unsigned int LENGTH_D_1CH_N3 = 8;
const unsigned int LENGTH_I_1CH_N2 = 2;
const unsigned int LENGTH_I_1CH_N3 = 4;
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryANYJ::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    frac ssz1 = I1.get_frac("SSZ");
    Invar Ip;

    switch(P::spin) {
     case 2:
      {      
        Ip = Invar(q1 - 1, ssz1 - HALF); // Spin-up creation operator
        {
 nrglog('f', "RECALC(fn=" << "anyj/any-n2-doublet1.dat" << ", len=" << ANYJ::LENGTH_D_1CH_N2 << ", Iop=" << Invar(ONE, HALF) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "anyj/any-n2-doublet1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_D_1CH_N2);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ANYJ::LENGTH_D_1CH_N2, Invar(ONE, HALF));
 }
};
        Ip = Invar(q1 - 1, ssz1 + HALF); // Spin-down creation operator
        {
 nrglog('f', "RECALC(fn=" << "anyj/any-n2-doublet-1.dat" << ", len=" << ANYJ::LENGTH_D_1CH_N2 << ", Iop=" << Invar(ONE, -HALF) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "anyj/any-n2-doublet-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_D_1CH_N2);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ANYJ::LENGTH_D_1CH_N2, Invar(ONE, -HALF));
 }
};
      }
      break;
      
     case 3:
      {
        Ip = Invar(q1 - 1, ssz1 - 1); // Spin 1 creation operator
	{
 nrglog('f', "RECALC(fn=" << "anyj/any-n3-doublet2.dat" << ", len=" << ANYJ::LENGTH_D_1CH_N3 << ", Iop=" << Invar(ONE, ONE) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "anyj/any-n3-doublet2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_D_1CH_N3);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ANYJ::LENGTH_D_1CH_N3, Invar(ONE, ONE));
 }
};
	Ip = Invar(q1 - 1, ssz1); // Spin 0 creation operator
	{
 nrglog('f', "RECALC(fn=" << "anyj/any-n3-doublet0.dat" << ", len=" << ANYJ::LENGTH_D_1CH_N3 << ", Iop=" << Invar(ONE, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "anyj/any-n3-doublet0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_D_1CH_N3);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ANYJ::LENGTH_D_1CH_N3, Invar(ONE, 0));
 }
};
	Ip = Invar(q1 - 1, ssz1 + 1); // Spin -1 creation operator
	{
 nrglog('f', "RECALC(fn=" << "anyj/any-n3-doublet-2.dat" << ", len=" << ANYJ::LENGTH_D_1CH_N3 << ", Iop=" << Invar(ONE, -ONE) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "anyj/any-n3-doublet-2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_D_1CH_N3);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, ANYJ::LENGTH_D_1CH_N3, Invar(ONE, -ONE));
 }
};
      }
      break;

     default:
        my_error("Not implemented.");
    }
      
  } // loop over invariant subspaces
}

// Driver routine for recalc_f()
void SymmetryANYJ::recalc_irreduc(const DiagInfo &diag)
{
  my_assert(P::channels == 1);

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    frac sszp = Ip.get_frac("SSZ");
    Invar I1;

    switch (P::spin) {
     case 2:
      {
        I1 = Invar(qp + 1, sszp + HALF);
        {
 nrglog('f', "RECALC_F(fn=" << "anyj/any-n2-spin1.dat" << ", ch=" << 0 << ", len=" << ANYJ::LENGTH_I_1CH_N2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "anyj/any-n2-spin1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_I_1CH_N2);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ANYJ::LENGTH_I_1CH_N2);
 }
};
        I1 = Invar(qp + 1, sszp - HALF);
        {
 nrglog('f', "RECALC_F(fn=" << "anyj/any-n2-spin-1.dat" << ", ch=" << 0 << ", len=" << ANYJ::LENGTH_I_1CH_N2 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "anyj/any-n2-spin-1.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_I_1CH_N2);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ANYJ::LENGTH_I_1CH_N2);
 }
};
      }
      break;
      
     case 3:
      {
        I1 = Invar(qp + 1, sszp + 1);
	{
 nrglog('f', "RECALC_F(fn=" << "anyj/any-n3-spin2.dat" << ", ch=" << 0 << ", len=" << ANYJ::LENGTH_I_1CH_N3 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "anyj/any-n3-spin2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_I_1CH_N3);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ANYJ::LENGTH_I_1CH_N3);
 }
};
        I1 = Invar(qp + 1, sszp);
	{
 nrglog('f', "RECALC_F(fn=" << "anyj/any-n3-spin0.dat" << ", ch=" << 0 << ", len=" << ANYJ::LENGTH_I_1CH_N3 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "anyj/any-n3-spin0.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_I_1CH_N3);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ANYJ::LENGTH_I_1CH_N3);
 }
};
        I1 = Invar(qp + 1, sszp - 1);
	{
 nrglog('f', "RECALC_F(fn=" << "anyj/any-n3-spin-2.dat" << ", ch=" << 0 << ", len=" << ANYJ::LENGTH_I_1CH_N3 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "anyj/any-n3-spin-2.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == ANYJ::LENGTH_I_1CH_N3);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, ANYJ::LENGTH_I_1CH_N3);
 }
};
      }
      break;

     default:
        my_error("Not implemented.");
    }
  } // loop over isp
}

