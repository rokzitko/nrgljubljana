// *** WARNING!!! Modify nrg-recalc-QSZTZ.cc.m4, not nrg-recalc-QSZTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

include(recalc-macros.m4)

namespace QSZTZ {
#include "qsztz/qsztz-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSZTZ::recalc_doublet(DiagInfo &diag,
                                   MatrixElements &cold,
                                   MatrixElements &cnew)
{
  nrglog('f', "QSZTZ::recalc_doublet() called");
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ssz1 = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;
    
    nrglog('f', "I1=" << I1);
    
    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet, 
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1-1, ssz1+1, tz1-1);
    RECALC_TAB("qsztz/qsztz-doubletp-1.dat", QSZTZ::LENGTH_D_3CH, Invar(1, -1, +1));
      
    Ip = Invar(q1-1, ssz1-1, tz1-1);
    RECALC_TAB("qsztz/qsztz-doubletm-1.dat", QSZTZ::LENGTH_D_3CH, Invar(1, +1, +1));

    Ip = Invar(q1-1, ssz1+1, tz1);
    RECALC_TAB("qsztz/qsztz-doubletp0.dat",  QSZTZ::LENGTH_D_3CH, Invar(1, -1, 0));
      
    Ip = Invar(q1-1, ssz1-1, tz1);
    RECALC_TAB("qsztz/qsztz-doubletm0.dat",  QSZTZ::LENGTH_D_3CH, Invar(1, +1, 0));

    Ip = Invar(q1-1, ssz1+1, tz1+1);
    RECALC_TAB("qsztz/qsztz-doubletp+1.dat", QSZTZ::LENGTH_D_3CH, Invar(1, -1, -1));
      
    Ip = Invar(q1-1, ssz1-1, tz1+1);
    RECALC_TAB("qsztz/qsztz-doubletm+1.dat", QSZTZ::LENGTH_D_3CH, Invar(1, +1, -1));
 }
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
void SymmetryQSZTZ::recalc_irreduc(const DiagInfo &diag)
{
  nrglog('f', "QSZTZ::recalc_irreduc() called");
  my_assert(!substeps);

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin sszp = Ip.get("SZ");
    Tangmom tzp = Ip.get("TZ");
    Invar I1;

    nrglog('f', "qp=" << qp << " sszp=" << sszp << " tzp=" << tzp);

    nrglog('f', "spinup+1");
    I1 = Invar(qp+1, sszp+1, tzp+1);
    RECALC_F_TAB("qsztz/qsztz-spinup+1.dat", 0, QSZTZ::LENGTH_I_3CH);
    
    nrglog('f', "spinup0");
    I1 = Invar(qp+1, sszp+1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spinup0.dat",  0, QSZTZ::LENGTH_I_3CH);

    nrglog('f', "spinup-1");
    I1 = Invar(qp+1, sszp+1, tzp-1);
    RECALC_F_TAB("qsztz/qsztz-spinup-1.dat", 0, QSZTZ::LENGTH_I_3CH);
    
    nrglog('f', "spindo+1");
    I1 = Invar(qp+1, sszp-1, tzp+1);
    RECALC_F_TAB("qsztz/qsztz-spindo+1.dat", 0, QSZTZ::LENGTH_I_3CH);

    nrglog('f', "spindo0");
    I1 = Invar(qp+1, sszp-1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spindo0.dat",  0, QSZTZ::LENGTH_I_3CH);
    
    nrglog('f', "spindo-1");
    I1 = Invar(qp+1, sszp-1, tzp-1);
    RECALC_F_TAB("qsztz/qsztz-spindo-1.dat", 0, QSZTZ::LENGTH_I_3CH);
  }
}


// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSZTZ::recalc_triplet(DiagInfo &diag,
                                  MatrixElements &cold,
                                  MatrixElements &cnew)
{
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ssz1 = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ssz1, tz1);
    RECALC_TAB("qsztz/qsztz-triplets.dat", QSZTZ::LENGTH_T0_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ssz1+2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletp.dat", QSZTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ssz1-2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletm.dat", QSZTZ::LENGTH_Tpm_3CH, Invar(0, 3, 0));
   }
}
