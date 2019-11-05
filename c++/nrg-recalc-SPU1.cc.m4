// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2008.
// This file pertains to (S) subspaces

namespace SPU1 {
#include "spu1/spu1-1ch-def.dat"
#include "spu1/spu1-2ch-def.dat"
}

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySPU1::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  if (!P::substeps) {
   LOOP(diag, is1) {
     Invar I1 = INVAR(is1);
     SZspin ssz1 = I1.get("SSZ");
     Invar Ip;

     Ip = Invar(ssz1+1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-doubletp.dat", SPU1::LENGTH_D_1CH, Invar(-1))',
           `RECALC_TAB("spu1/spu1-2ch-doubletp.dat", SPU1::LENGTH_D_2CH, Invar(-1))');
      
     Ip = Invar(ssz1-1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-doubletm.dat", SPU1::LENGTH_D_1CH, Invar(+1))',
            `RECALC_TAB("spu1/spu1-2ch-doubletm.dat", SPU1::LENGTH_D_2CH, Invar(+1))');
   } // loop
  } else {
    LOOP(diag, is1) {
     Invar I1 = INVAR(is1);
     SZspin ssz1 = I1.get("SSZ");
     Invar Ip;

     Ip = Invar(ssz1+1);
     RECALC_TAB("spu1/spu1-1ch-doubletp.dat", SPU1::LENGTH_D_1CH, Invar(-1));
      
     Ip = Invar(ssz1-1);
     RECALC_TAB("spu1/spu1-1ch-doubletm.dat", SPU1::LENGTH_D_1CH, Invar(+1));
    } // loop
  }
}

// Driver routine for recalc_f()
void SymmetrySPU1::recalc_irreduc(const DiagInfo &diag)
{
  my_assert(!P::substeps);

  LOOP_const(diag, isp) {
     Invar Ip = INVAR(isp);
     SZspin sszp = Ip.get("SSZ");
     Invar I1;

     I1 = Invar(sszp+1);
     ONETWO(`RECALC_F_TAB("spu1/spu1-1ch-spinupa.dat", 0, SPU1::LENGTH_I_1CH)',

            `RECALC_F_TAB("spu1/spu1-2ch-spinupa.dat", 0, SPU1::LENGTH_I_2CH);
 	     RECALC_F_TAB("spu1/spu1-2ch-spinupb.dat", 1, SPU1::LENGTH_I_2CH)');
    
     I1 = Invar(sszp-1);
     ONETWO(`RECALC_F_TAB("spu1/spu1-1ch-spindowna.dat", 0, SPU1::LENGTH_I_1CH)',

            `RECALC_F_TAB("spu1/spu1-2ch-spindowna.dat", 0, SPU1::LENGTH_I_2CH);
             RECALC_F_TAB("spu1/spu1-2ch-spindownb.dat", 1, SPU1::LENGTH_I_2CH)');
  } // loop
}

// Driver routine for recalc_f()
void SymmetrySPU1::recalc_irreduc_substeps(const DiagInfo &diag, int M)
{
  my_assert(P::substeps);
  
  LOOP_const(diag, isp) {
     Invar Ip = INVAR(isp);
     SZspin sszp = Ip.get("SSZ");
     Invar I1;

     I1 = Invar(sszp+1);
     RECALC_F_TAB("spu1/spu1-1ch-spinupa.dat", M, SPU1::LENGTH_I_1CH);

     I1 = Invar(sszp-1);
     RECALC_F_TAB("spu1/spu1-1ch-spindowna.dat", M, SPU1::LENGTH_I_1CH);
  } // loop
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetrySPU1::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  if (!P::substeps) {
   LOOP(diag, is1) {
     Invar I1 = INVAR(is1);
     SZspin ssz1 = I1.get("SSZ");
     Invar Ip;

     Ip = Invar(ssz1);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-triplets.dat", SPU1::LENGTH_T0_1CH, Invar(0))',
           `RECALC_TAB("spu1/spu1-2ch-triplets.dat", SPU1::LENGTH_T0_2CH, Invar(0))');
      
     Ip = Invar(ssz1+2);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-tripletp.dat", SPU1::LENGTH_Tpm_1CH, Invar(-2))',
           `RECALC_TAB("spu1/spu1-2ch-tripletp.dat", SPU1::LENGTH_Tpm_2CH, Invar(-2))');
      
     Ip = Invar(ssz1-2);
     ONETWO(`RECALC_TAB("spu1/spu1-1ch-tripletm.dat", SPU1::LENGTH_Tpm_1CH, Invar(+2))',
            `RECALC_TAB("spu1/spu1-2ch-tripletm.dat", SPU1::LENGTH_Tpm_2CH, Invar(+2))');
   } // loop
  } else {
   LOOP(diag, is1) {
     Invar I1 = INVAR(is1);
     SZspin ssz1 = I1.get("SSZ");
     Invar Ip;

     Ip = Invar(ssz1);
     RECALC_TAB("spu1/spu1-1ch-triplets.dat", SPU1::LENGTH_T0_1CH, Invar(0));
      
     Ip = Invar(ssz1+2);
     RECALC_TAB("spu1/spu1-1ch-tripletp.dat", SPU1::LENGTH_Tpm_1CH, Invar(-2));
      
     Ip = Invar(ssz1-2);
     RECALC_TAB("spu1/spu1-1ch-tripletm.dat", SPU1::LENGTH_Tpm_1CH, Invar(+2));
   } // loop
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

// 2-channel only

#undef QDIFF
#define QDIFF(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1UP
#define Q1UP(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2UP
#define Q2UP(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1DO
#define Q1DO(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2DO
#define Q2DO(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)


void SymmetrySPU1::recalc_global(DiagInfo &diag,
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
#include "spu1/spu1-1ch-Qtot.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-Qtot.dat"
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
#include "spu1/spu1-1ch-Iztot.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-Iztot.dat"
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
#include "spu1/spu1-1ch-Ixtot.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-Ixtot.dat"
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
#include "spu1/spu1-1ch-Iptot.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-Iptot.dat"
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
#include "spu1/spu1-1ch-Imtot.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-Imtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Qdiff") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-qdiff.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q1") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q1.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q2") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q2.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q1UP") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q1up.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q2UP") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q2up.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q1DO") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q1do.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Q2DO") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 2:
#include "spu1/spu1-2ch-q2do.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }


}

#undef Q2
