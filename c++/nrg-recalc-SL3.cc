// *** WARNING!!! Modify nrg-recalc-SL3.cc.m4, not nrg-recalc-SL3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007, Oct 2010
// This file pertains to the spinless-fermions code.

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace SL3 {
#include "sl3/sl3-3ch-def.dat"
}

// Recalculate matrix elements of a "doublet" tensor operator
void SymmetrySL3::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q11 = I1.get("Q1");
    Number q21 = I1.get("Q2");
    Number q31 = I1.get("Q3");
    Invar Ip = Invar(q11-1, q21, q31); // This is a channel 1 operator
    {
 nrglog('f', "RECALC(fn=" << "sl3/sl3-3ch-doublet.dat" << ", len=" << SL::LENGTH_D_3CH << ", Iop=" << Invar(1, 0, 0) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "sl3/sl3-3ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SL::LENGTH_D_3CH, Invar(1, 0, 0));
 }
};
  }
}

// Driver routine for recalc_f()
void SymmetrySL3::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number q1p = Ip.get("Q1");
    Number q2p = Ip.get("Q2");
    Number q3p = Ip.get("Q3");
    
    Invar I1;
    
    I1 = Invar(q1p+1, q2p, q3p);
    {
 nrglog('f', "RECALC_F(fn=" << "sl3/sl3-3ch-a.dat" << ", ch=" << 0 << ", len=" << SL3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl3/sl3-3ch-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SL3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(q1p, q2p+1, q3p);
    {
 nrglog('f', "RECALC_F(fn=" << "sl3/sl3-3ch-b.dat" << ", ch=" << 1 << ", len=" << SL3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl3/sl3-3ch-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SL3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(q1p, q2p, q3p+1);
    {
 nrglog('f', "RECALC_F(fn=" << "sl3/sl3-3ch-c.dat" << ", ch=" << 2 << ", len=" << SL3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl3/sl3-3ch-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SL3::LENGTH_I_3CH);
 }
};
  }
}

#undef QTOT
#define QTOT(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N1
#define N1(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N2
#define N2(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N3
#define N3(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetrySL3::recalc_global(DiagInfo &diag,
                               string name,
                               MatrixElements &cnew)
{
  if (name == "Qtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
#include "sl3/sl3-3ch-qtot.dat"
     } // LOOP
   }

   if (name == "N1") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
#include "sl3/sl3-3ch-N1.dat"
     } // LOOP
   }

   if (name == "N2") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
#include "sl3/sl3-3ch-N2.dat"
     } // LOOP
   }
   
   if (name == "N3") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
#include "sl3/sl3-3ch-N3.dat"
     } // LOOP
   }
}
