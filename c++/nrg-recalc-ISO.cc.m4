// *** WARNING!!! Modify nrg-recalc-ISO.cc.m4, not nrg-recalc-ISO.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (I,S) subspaces

include(recalc-macros.m4)

namespace ISO1 {
#include "iso/iso-1ch-def.dat"
#include "iso/iso-2ch-def.dat"
#include "iso/iso-3ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
void SymmetryISO::recalc_doublet(DiagInfo &diag,
                    	 MatrixElements &cold,
                    	 MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1-1, ss1+1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletmp.dat", ISO1::LENGTH_D_1CH, Invar(2, 2))',
           `RECALC_TAB("iso/iso-2ch-doubletmp.dat", ISO1::LENGTH_D_2CH, Invar(2, 2))');
      
    Ip = Invar(ii1-1, ss1-1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletmm.dat", ISO1::LENGTH_D_1CH, Invar(2, 2))',
    	   `RECALC_TAB("iso/iso-2ch-doubletmm.dat", ISO1::LENGTH_D_2CH, Invar(2, 2))');

    Ip = Invar(ii1+1, ss1+1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletpp.dat", ISO1::LENGTH_D_1CH, Invar(2, 2))',
    	   `RECALC_TAB("iso/iso-2ch-doubletpp.dat", ISO1::LENGTH_D_2CH, Invar(2, 2))');

    Ip = Invar(ii1+1, ss1-1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletpm.dat", ISO1::LENGTH_D_1CH, Invar(2, 2))',
   	   `RECALC_TAB("iso/iso-2ch-doubletpm.dat", ISO1::LENGTH_D_2CH, Invar(2, 2))');
  }
}

// (ISO): Four calls of recalc_f() are necessary for each channel.

// Driver routine for recalc_f()
void SymmetryISO::recalc_irreduc(const DiagInfo &diag)
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
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spinup-isoupa.dat", 0, ISO1::LENGTH_I_1CH)',

          `RECALC_F_TAB("iso/iso-2ch-spinup-isoupa.dat", 0, ISO1::LENGTH_I_2CH);
	   RECALC_F_TAB("iso/iso-2ch-spinup-isoupb.dat", 1, ISO1::LENGTH_I_2CH)',

          `RECALC_F_TAB("iso/iso-3ch-spinup-isoupa.dat", 0, ISO1::LENGTH_I_3CH_0);
	   RECALC_F_TAB("iso/iso-3ch-spinup-isoupb.dat", 1, ISO1::LENGTH_I_3CH_1);
	   RECALC_F_TAB("iso/iso-3ch-spinup-isoupc.dat", 2, ISO1::LENGTH_I_3CH_2)');
    
    I1 = Invar(iip+1, ssp-1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spindown-isoupa.dat", 0, ISO1::LENGTH_I_1CH)',

          `RECALC_F_TAB("iso/iso-2ch-spindown-isoupa.dat", 0, ISO1::LENGTH_I_2CH);
	   RECALC_F_TAB("iso/iso-2ch-spindown-isoupb.dat", 1, ISO1::LENGTH_I_2CH)',
	   
          `RECALC_F_TAB("iso/iso-3ch-spindown-isoupa.dat", 0, ISO1::LENGTH_I_3CH_0);
	   RECALC_F_TAB("iso/iso-3ch-spindown-isoupb.dat", 1, ISO1::LENGTH_I_3CH_1);
	   RECALC_F_TAB("iso/iso-3ch-spindown-isoupc.dat", 2, ISO1::LENGTH_I_3CH_2)');

    I1 = Invar(iip-1, ssp+1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spinup-isodowna.dat", 0, ISO1::LENGTH_I_1CH)',

          `RECALC_F_TAB("iso/iso-2ch-spinup-isodowna.dat", 0, ISO1::LENGTH_I_2CH);
	   RECALC_F_TAB("iso/iso-2ch-spinup-isodownb.dat", 1, ISO1::LENGTH_I_2CH)',
	   
          `RECALC_F_TAB("iso/iso-3ch-spinup-isodowna.dat", 0, ISO1::LENGTH_I_3CH_0);
	   RECALC_F_TAB("iso/iso-3ch-spinup-isodownb.dat", 1, ISO1::LENGTH_I_3CH_1);
	   RECALC_F_TAB("iso/iso-3ch-spinup-isodownc.dat", 2, ISO1::LENGTH_I_3CH_2)');

    I1 = Invar(iip-1, ssp-1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spindown-isodowna.dat", 0, ISO1::LENGTH_I_1CH)',
    
    	  `RECALC_F_TAB("iso/iso-2ch-spindown-isodowna.dat", 0, ISO1::LENGTH_I_2CH);
	   RECALC_F_TAB("iso/iso-2ch-spindown-isodownb.dat", 1, ISO1::LENGTH_I_2CH)',

          `RECALC_F_TAB("iso/iso-3ch-spindown-isodowna.dat", 0, ISO1::LENGTH_I_3CH_0);
	   RECALC_F_TAB("iso/iso-3ch-spindown-isodownb.dat", 1, ISO1::LENGTH_I_3CH_1);
	   RECALC_F_TAB("iso/iso-3ch-spindown-isodownc.dat", 2, ISO1::LENGTH_I_3CH_2)');
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryISO::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1, ss1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-triplets.dat", ISO1::LENGTH_T0_1CH, Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-triplets.dat", ISO1::LENGTH_T0_2CH, Invar(1, 3))');

    Ip = Invar(ii1, ss1+2);
    ONETWO(`RECALC_TAB("iso/iso-1ch-tripletp.dat", ISO1::LENGTH_Tpm_1CH, Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-tripletp.dat", ISO1::LENGTH_Tpm_2CH, Invar(1, 3))');

    Ip = Invar(ii1, ss1-2);
    ONETWO(`RECALC_TAB("iso/iso-1ch-tripletm.dat", ISO1::LENGTH_Tpm_1CH, Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-tripletm.dat", ISO1::LENGTH_Tpm_2CH, Invar(1, 3))');
  }
}
