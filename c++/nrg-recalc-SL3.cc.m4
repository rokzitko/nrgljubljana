// *** WARNING!!! Modify nrg-recalc-SL3.cc.m4, not nrg-recalc-SL3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007, Oct 2010
// This file pertains to the spinless-fermions code.

include(recalc-macros.m4)

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
    RECALC_TAB("sl3/sl3-3ch-doublet.dat", SL::LENGTH_D_3CH, Invar(1, 0, 0));
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
    RECALC_F_TAB("sl3/sl3-3ch-a.dat", 0, SL3::LENGTH_I_3CH);
    
    I1 = Invar(q1p, q2p+1, q3p);
    RECALC_F_TAB("sl3/sl3-3ch-b.dat", 1, SL3::LENGTH_I_3CH);
    
    I1 = Invar(q1p, q2p, q3p+1);
    RECALC_F_TAB("sl3/sl3-3ch-c.dat", 2, SL3::LENGTH_I_3CH);
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
