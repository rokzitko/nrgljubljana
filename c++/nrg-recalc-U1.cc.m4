// *** WARNING!!! Modify nrg-recalc-U1.cc.m4, not nrg-recalc-U1.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Oct 2012
// This file pertains to (Q) subspaces

namespace U1 {
#include "u1/u1-1ch-def.dat"
#include "u1/u1-2ch-def.dat"
#include "u1/u1-3ch-def.dat"
}

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
void SymmetryU1::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Invar Ip = Invar(q1-1);
    ONE23(`RECALC_TAB("u1/u1-1ch-doublet.dat", U1::LENGTH_D_1CH, Invar(1))',
    	  `RECALC_TAB("u1/u1-2ch-doublet.dat", U1::LENGTH_D_2CH, Invar(1))',
	  `RECALC_TAB("u1/u1-3ch-doublet.dat", U1::LENGTH_D_3CH, Invar(1))');
  }
}

// Override the recalc_f definition: we need to track the spin index of 
// the f-matrices.

define(`RECALC_F_TAB_U1',
{
 if (diag.count(I1)) {
   struct Recalc_f recalc_table[]={
   #include $1
           };
   BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == $4);
   recalc_f(diag, a.opch[$2][$3], Ip, I1, recalc_table, $4);
 }
})

// Driver routine for recalc_f()
void SymmetryU1::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Invar I1 = Invar(qp+1);
    ONE23(`RECALC_F_TAB_U1("u1/u1-1ch-a-DO.dat", 0, 1, U1::LENGTH_I_1CH);
           RECALC_F_TAB_U1("u1/u1-1ch-a-UP.dat", 0, 0, U1::LENGTH_I_1CH)',

          `RECALC_F_TAB_U1("u1/u1-2ch-a-DO.dat", 0, 1, U1::LENGTH_I_2CH);
	   RECALC_F_TAB_U1("u1/u1-2ch-b-DO.dat", 1, 1, U1::LENGTH_I_2CH);
	   RECALC_F_TAB_U1("u1/u1-2ch-a-UP.dat", 0, 0, U1::LENGTH_I_2CH);
	   RECALC_F_TAB_U1("u1/u1-2ch-b-UP.dat", 1, 0, U1::LENGTH_I_2CH)',

          `RECALC_F_TAB_U1("u1/u1-3ch-a-DO.dat", 0, 1, U1::LENGTH_I_3CH);
	   RECALC_F_TAB_U1("u1/u1-3ch-b-DO.dat", 1, 1, U1::LENGTH_I_3CH);
	   RECALC_F_TAB_U1("u1/u1-3ch-c-DO.dat", 2, 1, U1::LENGTH_I_3CH);
	   RECALC_F_TAB_U1("u1/u1-3ch-a-UP.dat", 0, 0, U1::LENGTH_I_3CH);
	   RECALC_F_TAB_U1("u1/u1-3ch-b-UP.dat", 1, 0, U1::LENGTH_I_3CH);
	   RECALC_F_TAB_U1("u1/u1-3ch-c-UP.dat", 2, 0, U1::LENGTH_I_3CH)');
  }
}

#undef SPINX
#define SPINX(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#ifdef NRG_COMPLEX
 #undef SPINY
 #define SPINY(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
 #undef Complex
 #define Complex(x,y) cmpl(x,y)
#endif // NRG_COMPLEX

void SymmetryU1::recalc_global(DiagInfo &diag,
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
#include "u1/u1-1ch-spinz.dat"
         break;
        case 2:
#include "u1/u1-2ch-spinz.dat"
         break;
        case 3:
#include "u1/u1-3ch-spinz.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

#ifdef NRG_COMPLEX
   if (name == "SYtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "u1/u1-1ch-spiny.dat"
         break;
        case 2:
#include "u1/u1-2ch-spiny.dat"
         break;
        case 3:
#include "u1/u1-3ch-spiny.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
#endif

   if (name == "SXtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "u1/u1-1ch-spinx.dat"
         break;
        case 2:
#include "u1/u1-2ch-spinx.dat"
         break;
        case 3:
#include "u1/u1-3ch-spinx.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
}
