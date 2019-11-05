// *** WARNING!!! Modify nrg-recalc-QSZ.cc.m4, not nrg-recalc-QSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, 2006-2017
// This file pertains to (Q,SZ) subspaces

include(recalc-macros.m4)

namespace QSZ {
#include "qsz/qsz-1ch-def.dat"
#include "qsz/qsz-2ch-def.dat"
#include "qsz/qsz-3ch-def.dat"
}

// NOTE: p is ket (right side), 1 is bra (left side). OP is sandwiched
// in between. Thus Q[p] + Q[op] = Q[1].

// Recalculate matrix elements of a doublet tensor operator
void SymmetryQSZ::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
 if (!substeps) {
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    // In the case of (Q,S_z) basis, spin up and spin down are not
    // equivalent. The distinction appears in this recalculation code, but
    // also during the evaluation of the spectral densities, where two
    // spin-diagonal spectral densities (and also spin-flip spectral
    // density) can be defined.

    Ip = Invar(q1-1, ssz1+1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-doubletp.dat", QSZ::LENGTH_D_1CH, Invar(1, -1))',
          `RECALC_TAB("qsz/qsz-2ch-doubletp.dat", QSZ::LENGTH_D_2CH, Invar(1, -1))',
          `RECALC_TAB("qsz/qsz-3ch-doubletp.dat", QSZ::LENGTH_D_3CH, Invar(1, -1))');
      
    Ip = Invar(q1-1, ssz1-1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-doubletm.dat", QSZ::LENGTH_D_1CH, Invar(1, +1))',
          `RECALC_TAB("qsz/qsz-2ch-doubletm.dat", QSZ::LENGTH_D_2CH, Invar(1, +1))',
          `RECALC_TAB("qsz/qsz-3ch-doubletm.dat", QSZ::LENGTH_D_3CH, Invar(1, +1))');
  } // loop
 } else { // substeps
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q1-1, ssz1+1);
    RECALC_TAB("qsz/qsz-1ch-doubletp.dat", QSZ::LENGTH_D_1CH, Invar(1, -1));

    Ip = Invar(q1-1, ssz1-1);
    RECALC_TAB("qsz/qsz-1ch-doubletm.dat", QSZ::LENGTH_D_1CH, Invar(1, +1));
  } // loop
 } // if
}

// Driver routine for recalc_f()
void SymmetryQSZ::recalc_irreduc(const DiagInfo &diag)
{
  my_assert(!substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels. 

    I1 = Invar(qp+1, sszp+1);
    ONE23(`RECALC_F_TAB("qsz/qsz-1ch-spinupa.dat", 0, QSZ::LENGTH_I_1CH)',

          `RECALC_F_TAB("qsz/qsz-2ch-spinupa.dat", 0, QSZ::LENGTH_I_2CH);
      	   RECALC_F_TAB("qsz/qsz-2ch-spinupb.dat", 1, QSZ::LENGTH_I_2CH)',

          `RECALC_F_TAB("qsz/qsz-3ch-spinupa.dat", 0, QSZ::LENGTH_I_3CH);
           RECALC_F_TAB("qsz/qsz-3ch-spinupb.dat", 1, QSZ::LENGTH_I_3CH);
      	   RECALC_F_TAB("qsz/qsz-3ch-spinupc.dat", 2, QSZ::LENGTH_I_3CH)');

    I1 = Invar(qp+1, sszp-1);
    ONE23(`RECALC_F_TAB("qsz/qsz-1ch-spindowna.dat", 0, QSZ::LENGTH_I_1CH)',
    
          `RECALC_F_TAB("qsz/qsz-2ch-spindowna.dat", 0, QSZ::LENGTH_I_2CH);
	   RECALC_F_TAB("qsz/qsz-2ch-spindownb.dat", 1, QSZ::LENGTH_I_2CH)',

          `RECALC_F_TAB("qsz/qsz-3ch-spindowna.dat", 0, QSZ::LENGTH_I_3CH);
           RECALC_F_TAB("qsz/qsz-3ch-spindownb.dat", 1, QSZ::LENGTH_I_3CH);
      	   RECALC_F_TAB("qsz/qsz-3ch-spindownc.dat", 2, QSZ::LENGTH_I_3CH)');
  } // loop
}

void SymmetryQSZ::recalc_irreduc_substeps(const DiagInfo &diag, int M)
{
  my_assert(substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(qp+1, sszp+1);
    RECALC_F_TAB("qsz/qsz-1ch-spinupa.dat", M, QSZ::LENGTH_I_1CH);

    I1 = Invar(qp+1, sszp-1);
    RECALC_F_TAB("qsz/qsz-1ch-spindowna.dat", M, QSZ::LENGTH_I_1CH);
  } // loop
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetryQSZ::recalc_triplet(DiagInfo &diag,
                                 MatrixElements &cold,
				 MatrixElements &cnew)
{
 if (!substeps) {
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    SZspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q1, ssz1);
    ONE23(`RECALC_TAB("qsz/qsz-1ch-triplets.dat", QSZ::LENGTH_T0_1CH, Invar(0, 0))',
          `RECALC_TAB("qsz/qsz-2ch-triplets.dat", QSZ::LENGTH_T0_2CH, Invar(0, 0))',
          `RECALC_TAB("qsz/qsz-3ch-triplets.dat", QSZ::LENGTH_T0_3CH, Invar(0, 0))');

    Ip = Invar(q1, ssz1+2);
    ONETWO(`RECALC_TAB("qsz/qsz-1ch-tripletp.dat", QSZ::LENGTH_Tpm_1CH, Invar(0, -2))',
           `RECALC_TAB("qsz/qsz-2ch-tripletp.dat", QSZ::LENGTH_Tpm_2CH, Invar(0, -2))',
           `RECALC_TAB("qsz/qsz-3ch-tripletp.dat", QSZ::LENGTH_Tpm_3CH, Invar(0, -2))');

    Ip = Invar(q1, ssz1-2);
    ONETWO(`RECALC_TAB("qsz/qsz-1ch-tripletm.dat", QSZ::LENGTH_Tpm_1CH, Invar(0, +2))',
           `RECALC_TAB("qsz/qsz-2ch-tripletm.dat", QSZ::LENGTH_Tpm_2CH, Invar(0, +2))',
           `RECALC_TAB("qsz/qsz-3ch-tripletm.dat", QSZ::LENGTH_Tpm_3CH, Invar(0, +2))');
  } // loop 
 } else { // substeps
  my_error("Not implemented.");
 }
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetryQSZ::recalc_global(DiagInfo &diag,
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
#include "qsz/qsz-1ch-spinz.dat"
	 break;
	case 2:
#include "qsz/qsz-2ch-spinz.dat"
	 break;
	case 3:
#include "qsz/qsz-3ch-spinz.dat"
	 break;
        default:
	 my_assert_not_reached();
       }
     } // LOOP
   }
}
