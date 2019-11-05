// *** WARNING!!! Modify nrg-recalc-SPU1LR.cc.m4, not nrg-recalc-SPU1LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, May 2010
// This file pertains to (SZ,P) subspaces

namespace SPU1LR {
#include "spu1lr/spu1lr-1ch-def.dat"
#include "spu1lr/spu1lr-2ch-def.dat"
}

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySPU1LR::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    SZspin ssz1 = I1.get("SSZ");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1+1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-doubletp.dat", SPU1LR::LENGTH_D_1CH, Invar(-1, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-doubletp.dat", SPU1LR::LENGTH_D_2CH, Invar(-1, 1))');
      
    Ip = Invar(ssz1-1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-doubletm.dat", SPU1LR::LENGTH_D_1CH, Invar(+1, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-doubletm.dat", SPU1LR::LENGTH_D_2CH, Invar(+1, 1))');
  }
}

// Driver routine for recalc_f()
void SymmetrySPU1LR::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    SZspin sszp = Ip.get("SSZ");
    int pp = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(sszp+1, pp);
    ONETWO(`RECALC_F_TAB("spu1lr/spu1lr-1ch-spinupa.dat", 0, SPU1LR::LENGTH_I_1CH)',

           `RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupa.dat", 0, SPU1LR::LENGTH_I_2CH);
	    RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupb.dat", 1, SPU1LR::LENGTH_I_2CH)');
    
    I1 = Invar(sszp-1, pp);
    ONETWO(`RECALC_F_TAB("spu1lr/spu1lr-1ch-spindowna.dat", 0, SPU1LR::LENGTH_I_1CH)',

           `RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowna.dat", 0, SPU1LR::LENGTH_I_2CH);
            RECALC_F_TAB("spu1lr/spu1lr-2ch-spindownb.dat", 1, SPU1LR::LENGTH_I_2CH)');

   // CASE II: OPPOSITE PARITY

    if (channels == 2) {
      I1 = Invar(sszp+1, -pp);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupdiffa.dat", 0, SPU1LR::LENGTH_I_2CH);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupdiffb.dat", 1, SPU1LR::LENGTH_I_2CH);
    
      I1 = Invar(sszp-1, -pp);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowndiffa.dat", 0, SPU1LR::LENGTH_I_2CH);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowndiffb.dat", 1, SPU1LR::LENGTH_I_2CH);
    }
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetrySPU1LR::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    SZspin ssz1 = I1.get("SSZ");
    int p1 = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-triplets.dat", SPU1LR::LENGTH_T0_1CH, Invar(0, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-triplets.dat", SPU1LR::LENGTH_T0_2CH, Invar(0, 1))');
      
    Ip = Invar(ssz1+2);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-tripletp.dat", SPU1LR::LENGTH_Tpm_1CH, Invar(-2, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-tripletp.dat", SPU1LR::LENGTH_Tpm_2CH, Invar(-2, 1))');
      
    Ip = Invar(ssz1-2);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-tripletm.dat", SPU1LR::LENGTH_Tpm_1CH, Invar(+2, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-tripletm.dat", SPU1LR::LENGTH_Tpm_2CH, Invar(+2, 1))');
  }
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * psgn(getnn()+1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * psgn(getnn()+1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * psgn(getnn()+1))


void SymmetrySPU1LR::recalc_global(DiagInfo &diag,
				 string name,
				 MatrixElements &cnew)
{
  if (name == "Qtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-Qtot.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-Qtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Iztot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-Iztot.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-Iztot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

if (name == "Ixtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-Ixtot.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-Ixtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Iptot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-Iptot.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-Iptot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Imtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-Imtot.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-Imtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

}
