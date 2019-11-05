// *** WARNING!!! Modify nrg-recalc-ISO2.cc.m4, not nrg-recalc-ISO2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, May 2008
// This file pertains to (I,S) subspaces
// Version for EVEN number of impurities.

include(recalc-macros.m4)

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
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletmp.dat", ISO2::LENGTH_D_1CH, Invar(2, 2))',
           `RECALC_TAB("iso2/iso2-2ch-doubletmp.dat", ISO2::LENGTH_D_2CH, Invar(2, 2))');
      
    Ip = Invar(ii1-1, ss1-1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletmm.dat", ISO2::LENGTH_D_1CH, Invar(2, 2))',
    	   `RECALC_TAB("iso2/iso2-2ch-doubletmm.dat", ISO2::LENGTH_D_2CH, Invar(2, 2))');

    Ip = Invar(ii1+1, ss1+1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletpp.dat", ISO2::LENGTH_D_1CH, Invar(2, 2))',
    	   `RECALC_TAB("iso2/iso2-2ch-doubletpp.dat", ISO2::LENGTH_D_2CH, Invar(2, 2))');

    Ip = Invar(ii1+1, ss1-1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletpm.dat", ISO2::LENGTH_D_1CH, Invar(2, 2))',
   	   `RECALC_TAB("iso2/iso2-2ch-doubletpm.dat", ISO2::LENGTH_D_2CH, Invar(2, 2))');
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
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spinup-isoupa.dat", 0, ISO2::LENGTH_I_1CH)',

           `RECALC_F_TAB("iso2/iso2-2ch-spinup-isoupa.dat", 0, ISO2::LENGTH_I_2CH);
     	    RECALC_F_TAB("iso2/iso2-2ch-spinup-isoupb.dat", 1, ISO2::LENGTH_I_2CH)');
    
    I1 = Invar(iip+1, ssp-1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spindown-isoupa.dat", 0, ISO2::LENGTH_I_1CH)',

           `RECALC_F_TAB("iso2/iso2-2ch-spindown-isoupa.dat", 0, ISO2::LENGTH_I_2CH);
	    RECALC_F_TAB("iso2/iso2-2ch-spindown-isoupb.dat", 1, ISO2::LENGTH_I_2CH)');

    I1 = Invar(iip-1, ssp+1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spinup-isodowna.dat", 0, ISO2::LENGTH_I_1CH)',

           `RECALC_F_TAB("iso2/iso2-2ch-spinup-isodowna.dat", 0, ISO2::LENGTH_I_2CH);
	    RECALC_F_TAB("iso2/iso2-2ch-spinup-isodownb.dat", 1, ISO2::LENGTH_I_2CH)');

    I1 = Invar(iip-1, ssp-1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spindown-isodowna.dat", 0, ISO2::LENGTH_I_1CH)',
    
    	   `RECALC_F_TAB("iso2/iso2-2ch-spindown-isodowna.dat", 0, ISO2::LENGTH_I_2CH);
	    RECALC_F_TAB("iso2/iso2-2ch-spindown-isodownb.dat", 1, ISO2::LENGTH_I_2CH)');
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
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-triplets.dat", ISO2::LENGTH_T0_1CH, Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-triplets.dat", ISO2::LENGTH_T0_2CH, Invar(1, 3))');

    Ip = Invar(ii1, ss1+2);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-tripletp.dat", ISO2::LENGTH_Tpm_1CH, Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-tripletp.dat", ISO2::LENGTH_Tpm_2CH, Invar(1, 3))');

    Ip = Invar(ii1, ss1-2);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-tripletm.dat", ISO2::LENGTH_Tpm_1CH, Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-tripletm.dat", ISO2::LENGTH_Tpm_2CH, Invar(1, 3))');
  }
}
