// *** WARNING!!! Modify nrg-recalc-SL.cc.m4, not nrg-recalc-SL.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007
// This file pertains to the spinless-fermions code.

include(recalc-macros.m4)

namespace SL {
#include "sl/sl-1ch-def.dat"
#include "sl/sl-2ch-def.dat"
#include "sl/sl-3ch-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySL::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Invar Ip = Invar(q1-1);
    ONE23(`RECALC_TAB("sl/sl-1ch-doublet.dat", SL::LENGTH_D_1CH, Invar(1))',
    	  `RECALC_TAB("sl/sl-2ch-doublet.dat", SL::LENGTH_D_2CH, Invar(1))',
          `RECALC_TAB("sl/sl-3ch-doublet.dat", SL::LENGTH_D_3CH, Invar(1))');
  }
}

// Driver routine for recalc_f()
void SymmetrySL::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Invar I1 = Invar(qp+1);
    ONE23(`RECALC_F_TAB("sl/sl-1ch-a.dat", 0, SL::LENGTH_I_1CH)',
    	  `RECALC_F_TAB("sl/sl-2ch-a.dat", 0, SL::LENGTH_I_2CH);
	   RECALC_F_TAB("sl/sl-2ch-b.dat", 1, SL::LENGTH_I_2CH)',
    	  `RECALC_F_TAB("sl/sl-3ch-a.dat", 0, SL::LENGTH_I_3CH);
    	   RECALC_F_TAB("sl/sl-3ch-b.dat", 1, SL::LENGTH_I_3CH);
	   RECALC_F_TAB("sl/sl-3ch-c.dat", 2, SL::LENGTH_I_3CH)'
	   );
  }
}

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef QTOT
#define QTOT(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N1
#define N1(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N2
#define N2(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N3
#define N3(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetrySL::recalc_global(DiagInfo &diag,
                               string name,
                               MatrixElements &cnew)
{
   if (name == "Qdiff") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 2:
#include "sl/sl-2ch-qdiff.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

  if (name == "Qtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 2:
#include "sl/sl-2ch-qtot.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

   if (name == "N1") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N1.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

   if (name == "N2") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N2.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
   
   if (name == "N3") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N3.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
}
