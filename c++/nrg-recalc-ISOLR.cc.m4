// *** WARNING!!! Modify nrg-recalc-ISOLR.cc.m4, not nrg-recalc-ISOLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Nov 2007
// This file pertains to (I,S,P) subspaces

include(recalc-macros.m4)

namespace ISOLR {
#include "isolr/isolr-2ch-def.dat"
}

// (ISOLR): 8 calls of recalc_f() are necessary: different parities are also possible!

// Driver routine for recalc_f()
void SymmetryISOLR::recalc_irreduc(const DiagInfo &diag)
{
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;
    
    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels. 

    Ispin iip = Ip.get("II");
    Sspin ssp = Ip.get("SS");

    // IMPORTANT NEW ELEMENT: the parity is important here!! 
    int pp = Ip.get("P");
       
    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    // Both parities yield non-zero <I+-1/2, S+-1/2, P| a^\mu_\nu
    // |I,S,P'>.  Coefficients *DO* depend on P,P', or more accurately,
    // on whether or not P and P' are the same.
    // 
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip+1, ssp+1, pp);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isoupa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isoupb.dat", 1, ISOLR::LENGTH_I_2CH);
    
    I1 = Invar(iip+1, ssp-1, pp);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isoupa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isoupb.dat", 1, ISOLR::LENGTH_I_2CH);

    I1 = Invar(iip-1, ssp+1, pp);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isodowna.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isodownb.dat", 1, ISOLR::LENGTH_I_2CH);

    I1 = Invar(iip-1, ssp-1, pp);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isodowna.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isodownb.dat", 1, ISOLR::LENGTH_I_2CH);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip+1, ssp+1, -pp);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isoupdiffa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isoupdiffb.dat", 1, ISOLR::LENGTH_I_2CH);
    
    I1 = Invar(iip+1, ssp-1, -pp);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isoupdiffa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isoupdiffb.dat", 1, ISOLR::LENGTH_I_2CH);
      
    I1 = Invar(iip-1, ssp+1, -pp);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isodowndiffa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spinup-isodowndiffb.dat", 1, ISOLR::LENGTH_I_2CH);

    I1 = Invar(iip-1, ssp-1, -pp);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isodowndiffa.dat", 0, ISOLR::LENGTH_I_2CH);
    RECALC_F_TAB("isolr/isolr-2ch-spindown-isodowndiffb.dat", 1, ISOLR::LENGTH_I_2CH);
  }
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
void SymmetryISOLR::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1-1, ss1+1, p1);
    RECALC_TAB("isolr/isolr-2ch-doubletmp.dat", ISOLR::LENGTH_D_2CH, 
    						Invar(2, 2, +1));

    Ip = Invar(ii1-1, ss1-1, p1);
    RECALC_TAB("isolr/isolr-2ch-doubletmm.dat", ISOLR::LENGTH_D_2CH, 
    						Invar(2, 2, +1));

    Ip = Invar(ii1+1, ss1+1, p1);
    RECALC_TAB("isolr/isolr-2ch-doubletpp.dat", ISOLR::LENGTH_D_2CH, 
    						Invar(2, 2, +1));

    Ip = Invar(ii1+1, ss1-1, p1);
    RECALC_TAB("isolr/isolr-2ch-doubletpm.dat", ISOLR::LENGTH_D_2CH, 
    						Invar(2, 2, +1));
  }
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
void SymmetryISOLR::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ss1, p1);
    RECALC_TAB("isolr/isolr-2ch-triplets.dat", ISOLR::LENGTH_T0_2CH, 
    					       Invar(1, 3, +1));

    Ip = Invar(ii1, ss1+2, p1);
    RECALC_TAB("isolr/isolr-2ch-tripletp.dat", ISOLR::LENGTH_Tpm_2CH, 
    					       Invar(1, 3, +1));

    Ip = Invar(ii1, ss1-2, p1);
    RECALC_TAB("isolr/isolr-2ch-tripletm.dat", ISOLR::LENGTH_Tpm_2CH, 
    					       Invar(1, 3, +1));
  }
}
