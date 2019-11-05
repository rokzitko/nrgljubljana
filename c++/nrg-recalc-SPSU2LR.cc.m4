// *** WARNING!!! Modify nrg-recalc-SPSU2LR.cc.m4, not nrg-recalc-SPSU2LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Sep 2015
// This file pertains to (S,P) subspaces

namespace SPSU2LR {
#include "spsu2lr/spsu2lr-2ch-def.dat"
}

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySPSU2LR::recalc_doublet(DiagInfo &diag,
                                     MatrixElements &cold,
				     MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1+1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-doubletp.dat", SPSU2LR::LENGTH_D_2CH, Invar(-1, 1));
      
    Ip = Invar(ss1-1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-doubletm.dat", SPSU2LR::LENGTH_D_2CH, Invar(+1, 1));
  }
}

// Driver routine for recalc_f()
void SymmetrySPSU2LR::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Sspin ssp = Ip.get("SS");
    int pp = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(ssp+1, pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupa.dat", 0, SPSU2LR::LENGTH_I_2CH);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupb.dat", 1, SPSU2LR::LENGTH_I_2CH);
    
    I1 = Invar(ssp-1, pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowna.dat", 0, SPSU2LR::LENGTH_I_2CH);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindownb.dat", 1, SPSU2LR::LENGTH_I_2CH);

    // CASE II: OPPOSITE PARITY

    I1 = Invar(ssp+1, -pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupdiffa.dat", 0, SPSU2LR::LENGTH_I_2CH);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spinupdiffb.dat", 1, SPSU2LR::LENGTH_I_2CH);
    
    I1 = Invar(ssp-1, -pp);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowndiffa.dat", 0, SPSU2LR::LENGTH_I_2CH);
    RECALC_F_TAB("spsu2lr/spsu2lr-2ch-spindowndiffb.dat", 1, SPSU2LR::LENGTH_I_2CH);
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetrySPSU2LR::recalc_triplet(DiagInfo &diag,
                                     MatrixElements &cold,
                                     MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-triplets.dat", SPSU2LR::LENGTH_T0_2CH, Invar(0, 1));
      
    Ip = Invar(ss1+2, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-tripletp.dat", SPSU2LR::LENGTH_Tpm_2CH, Invar(-2, 1));
      
    Ip = Invar(ss1-2, p1);
    RECALC_TAB("spsu2lr/spsu2lr-2ch-tripletm.dat", SPSU2LR::LENGTH_Tpm_2CH, Invar(+2, 1));
  }
}
