// *** WARNING!!! Modify nrg-recalc-QSTZ.cc.m4, not nrg-recalc-QSTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2016
// This file pertains to (Q,S,Tz) subspaces

include(recalc-macros.m4)

namespace QSTZ {
#include "qstz/qstz-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSTZ::recalc_doublet(DiagInfo &diag,
                                  MatrixElements &cold,
                                  MatrixElements &cnew)
{
  nrglog('f', "QSTZ::recalc_doublet() called");
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;
    
    nrglog('f', "I1=" << I1);
    
    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet, 
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1-1, ss1+1, tz1-1);
    RECALC_TAB("qstz/qstz-doubletp-1.dat", QSTZ::LENGTH_D_3CH, Invar(1, 2, +1));
      
    Ip = Invar(q1-1, ss1-1, tz1-1);
    RECALC_TAB("qstz/qstz-doubletm-1.dat", QSTZ::LENGTH_D_3CH, Invar(1, 2, +1));

    Ip = Invar(q1-1, ss1+1, tz1);
    RECALC_TAB("qstz/qstz-doubletp0.dat",  QSTZ::LENGTH_D_3CH, Invar(1, 2, 0));
      
    Ip = Invar(q1-1, ss1-1, tz1);
    RECALC_TAB("qstz/qstz-doubletm0.dat",  QSTZ::LENGTH_D_3CH, Invar(1, 2, 0));

    Ip = Invar(q1-1, ss1+1, tz1+1);
    RECALC_TAB("qstz/qstz-doubletp+1.dat", QSTZ::LENGTH_D_3CH, Invar(1, 2, -1));
      
    Ip = Invar(q1-1, ss1-1, tz1+1);
    RECALC_TAB("qstz/qstz-doubletm+1.dat", QSTZ::LENGTH_D_3CH, Invar(1, 2, -1));
 }
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
void SymmetryQSTZ::recalc_irreduc(const DiagInfo &diag)
{
  nrglog('f', "QSTZ::recalc_irreduc() called");
  my_assert(!substeps);

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Tangmom tzp = Ip.get("TZ");
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces 
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    nrglog('f', "qp=" << qp << " ssp=" << ssp << " tzp=" << tzp);
    nrglog('f', "spinup+1");
    I1 = Invar(qp+1, ssp+1, tzp+1);
    RECALC_F_TAB("qstz/qstz-spinup+1.dat", 0, QSTZ::LENGTH_I_3CH);
    
    nrglog('f', "spinup0");
    I1 = Invar(qp+1, ssp+1, tzp);
    RECALC_F_TAB("qstz/qstz-spinup0.dat",  0, QSTZ::LENGTH_I_3CH);

    nrglog('f', "spinup-1");
    I1 = Invar(qp+1, ssp+1, tzp-1);
    RECALC_F_TAB("qstz/qstz-spinup-1.dat", 0, QSTZ::LENGTH_I_3CH);
    
    nrglog('f', "spindo+1");
    I1 = Invar(qp+1, ssp-1, tzp+1);
    RECALC_F_TAB("qstz/qstz-spindo+1.dat", 0, QSTZ::LENGTH_I_3CH);

    nrglog('f', "spindo0");
    I1 = Invar(qp+1, ssp-1, tzp);
    RECALC_F_TAB("qstz/qstz-spindo0.dat",  0, QSTZ::LENGTH_I_3CH);
    
    nrglog('f', "spindo-1");
    I1 = Invar(qp+1, ssp-1, tzp-1);
    RECALC_F_TAB("qstz/qstz-spindo-1.dat", 0, QSTZ::LENGTH_I_3CH);
  }
}


// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSTZ::recalc_triplet(DiagInfo &diag,
                                  MatrixElements &cold,
                                  MatrixElements &cnew)
{
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ss1, tz1);
    RECALC_TAB("qstz/qstz-triplets.dat", QSTZ::LENGTH_T0_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ss1+2, tz1);
    RECALC_TAB("qstz/qstz-tripletp.dat", QSTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ss1-2, tz1);
    RECALC_TAB("qstz/qstz-tripletm.dat", QSTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
   }
}


