// *** WARNING!!! Modify nrg-recalc-QST.cc.m4, not nrg-recalc-QST.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (Q,S,T) subspaces

include(recalc-macros.m4)

namespace QST {
#include "qst/qst-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQST::recalc_doublet(DiagInfo &diag,
       		                 MatrixElements &cold,
                                 MatrixElements &cnew)
{
  nrglog('f', "QST::recalc_doublet() called");
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double T = t1; // trick!
    double S = (ss1-1.)/2.;
    Invar Ip;
    
    nrglog('f', "I1=" << I1);
    
    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,1) is correct. 1 = add charge, 2 = doublet, 
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1-1, ss1+1, t1-1);
    RECALC_TAB("qst/qst-doubletp-1.dat", QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
      
    Ip = Invar(q1-1, ss1-1, t1-1);
    RECALC_TAB("qst/qst-doubletm-1.dat", QST::LENGTH_D_3CH_a, Invar(1, 2, 1));

    Ip = Invar(q1-1, ss1+1, t1);
    RECALC_TAB("qst/qst-doubletp0.dat",  QST::LENGTH_D_3CH_b, Invar(1, 2, 1));
      
    Ip = Invar(q1-1, ss1-1, t1);
    RECALC_TAB("qst/qst-doubletm0.dat",  QST::LENGTH_D_3CH_b, Invar(1, 2, 1));

    Ip = Invar(q1-1, ss1+1, t1+1);
    RECALC_TAB("qst/qst-doubletp+1.dat", QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
      
    Ip = Invar(q1-1, ss1-1, t1+1);
    RECALC_TAB("qst/qst-doubletm+1.dat", QST::LENGTH_D_3CH_a, Invar(1, 2, 1));
 }
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
void SymmetryQST::recalc_irreduc(const DiagInfo &diag)
{
  nrglog('f', "QST::recalc_irreduc() called");
  my_assert(!substeps);

  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Tangmom tp = Ip.get("T");
    double T = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces 
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    nrglog('f', "qp=" << qp << " ssp=" << ssp << " tp=" << tp);
    nrglog('f', "spinup+1");
    I1 = Invar(qp+1, ssp+1, tp+1);
    RECALC_F_TAB("qst/qst-spinup+1.dat", 0, QST::LENGTH_I_3CH_0);
    
    nrglog('f', "spinup0");
    I1 = Invar(qp+1, ssp+1, tp);
    RECALC_F_TAB("qst/qst-spinup0.dat",  0, QST::LENGTH_I_3CH_1);

    nrglog('f', "spinup-1");
    I1 = Invar(qp+1, ssp+1, tp-1);
    RECALC_F_TAB("qst/qst-spinup-1.dat", 0, QST::LENGTH_I_3CH_2);
    
    nrglog('f', "spindo+1");
    I1 = Invar(qp+1, ssp-1, tp+1);
    RECALC_F_TAB("qst/qst-spindo+1.dat",  0, QST::LENGTH_I_3CH_0);

    nrglog('f', "spindo0");
    I1 = Invar(qp+1, ssp-1, tp);
    RECALC_F_TAB("qst/qst-spindo0.dat",   0, QST::LENGTH_I_3CH_1);
    
    nrglog('f', "spindo-1");
    I1 = Invar(qp+1, ssp-1, tp-1);
    RECALC_F_TAB("qst/qst-spindo-1.dat",  0, QST::LENGTH_I_3CH_2);
  }
}


// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQST::recalc_triplet(DiagInfo &diag,
                                MatrixElements &cold,
                                MatrixElements &cnew)
{
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S = (ss1-1.)/2.;
    double T = t1; // trick!
    Invar Ip;

    Ip = Invar(q1, ss1, t1);
    RECALC_TAB("qst/qst-triplets.dat", QST::LENGTH_T0_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ss1+2, t1);
    RECALC_TAB("qst/qst-tripletp.dat", QST::LENGTH_Tpm_3CH, Invar(0, 3, 0));
      
    Ip = Invar(q1, ss1-2, t1);
    RECALC_TAB("qst/qst-tripletm.dat", QST::LENGTH_Tpm_3CH, Invar(0, 3, 0));
   }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQST::recalc_orb_triplet(DiagInfo &diag,
                                     MatrixElements &cold,
                                     MatrixElements &cnew)
{
   nrglog('r', "recalc_orb_triplet");

   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S = (ss1-1.)/2.;
    double T = t1; // trick!
    Invar Ip;

    // 0 = chargeless
    // 1 = spin singlet (deg=2S+1=1)
    // 1 = spin triplet (T=1)
    
    nrglog('r', "orb-triplets");
    Ip = Invar(q1, ss1, t1);
    RECALC_TAB("qst/qst-orb-triplets.dat", QST::LENGTH_OT0_3CH, Invar(0, 1, 1));
      
    nrglog('r', "orb-tripletp");
    Ip = Invar(q1, ss1, t1+1);
    RECALC_TAB("qst/qst-orb-tripletp.dat", QST::LENGTH_OTpm_3CH, Invar(0, 1, 1));
      
    nrglog('r', "orb-tripletm");
    Ip = Invar(q1, ss1, t1-1);
    RECALC_TAB("qst/qst-orb-tripletm.dat", QST::LENGTH_OTpm_3CH, Invar(0, 1, 1));
   }
}


