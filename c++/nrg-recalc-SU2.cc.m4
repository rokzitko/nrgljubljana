// *** WARNING!!! Modify nrg-recalc-SU2.cc.m4, not nrg-recalc-SU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Sep 2009
// This file pertains to (I) subspaces

include(recalc-macros.m4)

namespace SU2 {
#include "su2/su2-1ch-def.dat"
#include "su2/su2-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
void SymmetrySU2::recalc_doublet(DiagInfo &diag,
                    	 MatrixElements &cold,
                    	 MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Invar Ip;

    Ip = Invar(ii1-1);
    ONETWO(`RECALC_TAB("su2/su2-1ch-doubletm.dat", SU2::LENGTH_D_1CH, Invar(2))',
           `RECALC_TAB("su2/su2-2ch-doubletm.dat", SU2::LENGTH_D_2CH, Invar(2))');
      
    Ip = Invar(ii1+1);
    ONETWO(`RECALC_TAB("su2/su2-1ch-doubletp.dat", SU2::LENGTH_D_1CH, Invar(2))',
    	   `RECALC_TAB("su2/su2-2ch-doubletp.dat", SU2::LENGTH_D_2CH, Invar(2))');
  }
}

// Override the recalc_f definition: we need to track the type (1 or 2) of
// the f-matrices.

define(`RECALC_F_TAB_SU2',
{
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include $1
                  };
		     BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == $4);
		        recalc_f(diag, a.opch[$2][$3], Ip, I1,
recalc_table, $4);
 }
 })
 

// Driver routine for recalc_f()
void SymmetrySU2::recalc_irreduc(const DiagInfo &diag)
{
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;

    Ispin iip = Ip.get("II");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(iip+1);
    ONETWO(`RECALC_F_TAB_SU2("su2/su2-1ch-type1-isoup-a.dat", 0, 0, SU2::LENGTH_I_1CH);
            RECALC_F_TAB_SU2("su2/su2-1ch-type2-isoup-a.dat", 0, 1, SU2::LENGTH_I_1CH)',

           `RECALC_F_TAB_SU2("su2/su2-2ch-type1-isoup-a.dat", 0, 0, SU2::LENGTH_I_2CH);
	    RECALC_F_TAB_SU2("su2/su2-2ch-type1-isoup-b.dat", 1, 0, SU2::LENGTH_I_2CH);
            RECALC_F_TAB_SU2("su2/su2-2ch-type2-isoup-a.dat", 0, 1, SU2::LENGTH_I_2CH);
	    RECALC_F_TAB_SU2("su2/su2-2ch-type2-isoup-b.dat", 1, 1, SU2::LENGTH_I_2CH)');
    
    I1 = Invar(iip-1);
    ONETWO(`RECALC_F_TAB_SU2("su2/su2-1ch-type1-isodown-a.dat", 0, 0, SU2::LENGTH_I_1CH);
            RECALC_F_TAB_SU2("su2/su2-1ch-type2-isodown-a.dat", 0, 1, SU2::LENGTH_I_1CH)',

           `RECALC_F_TAB_SU2("su2/su2-2ch-type1-isodown-a.dat", 0, 0, SU2::LENGTH_I_2CH);
	    RECALC_F_TAB_SU2("su2/su2-2ch-type1-isodown-b.dat", 1, 0, SU2::LENGTH_I_2CH);
            RECALC_F_TAB_SU2("su2/su2-2ch-type2-isodown-a.dat", 0, 1, SU2::LENGTH_I_2CH);
	    RECALC_F_TAB_SU2("su2/su2-2ch-type2-isodown-b.dat", 1, 1, SU2::LENGTH_I_2CH)');
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

void SymmetrySU2::recalc_global(DiagInfo &diag,
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
#include "su2/su2-1ch-spinz.dat"
         break;
        case 2:
#include "su2/su2-2ch-spinz.dat"
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
#include "su2/su2-1ch-spiny.dat"
         break;
        case 2:
#include "su2/su2-2ch-spiny.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
#endif // NRG_COMPLEX

   if (name == "SXtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "su2/su2-1ch-spinx.dat"
         break;
        case 2:
#include "su2/su2-2ch-spinx.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
}
