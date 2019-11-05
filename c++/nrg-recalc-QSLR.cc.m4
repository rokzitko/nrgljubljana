// *** WARNING!!! Modify nrg-recalc-QSLR.cc.m4, not nrg-recalc-QSLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

include(recalc-macros.m4)

namespace QSLR {
#include "qslr/qslr-2ch-def.dat"
}

// Driver routine for recalc_f()
void SymmetryQSLR::recalc_irreduc(const DiagInfo &diag)
{
  // CONVENTION: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels. 

    // IMPORTANT NEW ELEMENT: the parity is important here!! 
    int lrp = Ip.get("P");
       
    // Both parities yield non-zero <Q+1, S+-1/2, P| a^\dag_\nu
    // |Q,S,P'>.  Coefficients *DO* depend on P,P', or more
    // accurately, on whether or not P and P' are the same.
    // 
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(qp+1, ssp+1, lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spinupa.dat", 0, QSLR::LENGTH_I_2CH);
    RECALC_F_TAB("qslr/qslr-2ch-spinupb.dat", 1, QSLR::LENGTH_I_2CH);
    
    I1 = Invar(qp+1, ssp-1, lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spindowna.dat", 0, QSLR::LENGTH_I_2CH);
    RECALC_F_TAB("qslr/qslr-2ch-spindownb.dat", 1, QSLR::LENGTH_I_2CH);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp+1, ssp+1, -lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spinupdiffa.dat", 0, QSLR::LENGTH_I_2CH);
    RECALC_F_TAB("qslr/qslr-2ch-spinupdiffb.dat", 1, QSLR::LENGTH_I_2CH);
    
    I1 = Invar(qp+1, ssp-1, -lrp);
    RECALC_F_TAB("qslr/qslr-2ch-spindowndiffa.dat", 0, QSLR::LENGTH_I_2CH);
    RECALC_F_TAB("qslr/qslr-2ch-spindowndiffb.dat", 1, QSLR::LENGTH_I_2CH);
  }
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
void SymmetryQSLR::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1-1, ss1+1, p1);
    RECALC_TAB("qslr/qslr-2ch-doubletp.dat", QSLR::LENGTH_D_2CH, 
    					     Invar(1, 2, +1));

    Ip = Invar(q1-1, ss1-1, p1);
    RECALC_TAB("qslr/qslr-2ch-doubletm.dat", QSLR::LENGTH_D_2CH, 
    					     Invar(1, 2, +1));
  }
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
void SymmetryQSLR::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    RECALC_TAB("qslr/qslr-2ch-triplets.dat", QSLR::LENGTH_T0_2CH, 
    					     Invar(0, 3, +1));

    Ip = Invar(q1, ss1+2, p1);
    RECALC_TAB("qslr/qslr-2ch-tripletp.dat", QSLR::LENGTH_Tpm_2CH, 
    					     Invar(0, 3, +1));

    Ip = Invar(q1, ss1-2, p1);
    RECALC_TAB("qslr/qslr-2ch-tripletm.dat", QSLR::LENGTH_Tpm_2CH, 
    					     Invar(0, 3, +1));
  }
}
