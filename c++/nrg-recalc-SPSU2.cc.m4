// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Dec 2007
// This file pertains to (S) subspaces

namespace SPSU2 {
#include "spsu2/spsu2-1ch-def.dat"
#include "spsu2/spsu2-2ch-def.dat"
#include "spsu2/spsu2-3ch-def.dat"
}

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySPSU2::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  if (!substeps) {
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ss1+1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-doubletp.dat", SPSU2::LENGTH_D_1CH, Invar(2))',
          `RECALC_TAB("spsu2/spsu2-2ch-doubletp.dat", SPSU2::LENGTH_D_2CH, Invar(2))',
          `RECALC_TAB("spsu2/spsu2-3ch-doubletp.dat", SPSU2::LENGTH_D_3CH, Invar(2))');
      
    Ip = Invar(ss1-1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-doubletm.dat", SPSU2::LENGTH_D_1CH, Invar(2))',
          `RECALC_TAB("spsu2/spsu2-2ch-doubletm.dat", SPSU2::LENGTH_D_2CH, Invar(2))',
          `RECALC_TAB("spsu2/spsu2-3ch-doubletm.dat", SPSU2::LENGTH_D_3CH, Invar(2))');
   }
  } else {
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    Invar Ip;
		  
    Ip = Invar(ss1+1);
    RECALC_TAB("spsu2/spsu2-1ch-doubletp.dat", SPSU2::LENGTH_D_1CH, Invar(2));

    Ip = Invar(ss1-1);
    RECALC_TAB("spsu2/spsu2-1ch-doubletm.dat", SPSU2::LENGTH_D_1CH, Invar(2));
   }
  }
}

// Driver routine for recalc_f()
void SymmetrySPSU2::recalc_irreduc(const DiagInfo &diag)
{
  my_assert(!substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp+1);
    ONE23(`RECALC_F_TAB("spsu2/spsu2-1ch-spinupa.dat", 0, SPSU2::LENGTH_I_1CH)',
          `RECALC_F_TAB("spsu2/spsu2-2ch-spinupa.dat", 0, SPSU2::LENGTH_I_2CH);
	   RECALC_F_TAB("spsu2/spsu2-2ch-spinupb.dat", 1, SPSU2::LENGTH_I_2CH)',
          `RECALC_F_TAB("spsu2/spsu2-3ch-spinupa.dat", 0, SPSU2::LENGTH_I_3CH_0);
	   RECALC_F_TAB("spsu2/spsu2-3ch-spinupb.dat", 1, SPSU2::LENGTH_I_3CH_1);
	   RECALC_F_TAB("spsu2/spsu2-3ch-spinupc.dat", 2, SPSU2::LENGTH_I_3CH_2)');
    
    I1 = Invar(ssp-1);
    ONE23(`RECALC_F_TAB("spsu2/spsu2-1ch-spindowna.dat", 0, SPSU2::LENGTH_I_1CH)',
          `RECALC_F_TAB("spsu2/spsu2-2ch-spindowna.dat", 0, SPSU2::LENGTH_I_2CH);
           RECALC_F_TAB("spsu2/spsu2-2ch-spindownb.dat", 1, SPSU2::LENGTH_I_2CH)',
          `RECALC_F_TAB("spsu2/spsu2-3ch-spindowna.dat", 0, SPSU2::LENGTH_I_3CH_0);
           RECALC_F_TAB("spsu2/spsu2-3ch-spindownb.dat", 1, SPSU2::LENGTH_I_3CH_1);
           RECALC_F_TAB("spsu2/spsu2-3ch-spindownc.dat", 2, SPSU2::LENGTH_I_3CH_2)');
	   
	   // Note: for 3ch cases, the lengths in the three channels are not the same!
	   // The same thing occurs for all SU(2)_spin cases, for instance for symtype=QS.
	   // RZ, oct 2015
   }
}

// Driver routine for recalc_f()
void SymmetrySPSU2::recalc_irreduc_substeps(const DiagInfo &diag, int M)
{
  my_assert(substeps);
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp+1);
    RECALC_F_TAB("spsu2/spsu2-1ch-spinupa.dat", M, SPSU2::LENGTH_I_1CH);
    
    I1 = Invar(ssp-1);
    RECALC_F_TAB("spsu2/spsu2-1ch-spindowna.dat", M, SPSU2::LENGTH_I_1CH);
  }
}

// Recalculate matrix elements of a triplet tenzor operator
void SymmetrySPSU2::recalc_triplet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  if (!substeps) {
   LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ss1);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-triplets.dat", SPSU2::LENGTH_T0_1CH, Invar(3))',
          `RECALC_TAB("spsu2/spsu2-2ch-triplets.dat", SPSU2::LENGTH_T0_2CH, Invar(3))',
          `RECALC_TAB("spsu2/spsu2-3ch-triplets.dat", SPSU2::LENGTH_T0_3CH, Invar(3))');
      
    Ip = Invar(ss1+2);
    ONE23(`RECALC_TAB("spsu2/spsu2-1ch-tripletp.dat", SPSU2::LENGTH_Tpm_1CH, Invar(3))',
          `RECALC_TAB("spsu2/spsu2-2ch-tripletp.dat", SPSU2::LENGTH_Tpm_2CH, Invar(3))',
          `RECALC_TAB("spsu2/spsu2-3ch-tripletp.dat", SPSU2::LENGTH_Tpm_3CH, Invar(3))');
      
    Ip = Invar(ss1-2);
    ONETWO(`RECALC_TAB("spsu2/spsu2-1ch-tripletm.dat", SPSU2::LENGTH_Tpm_1CH, Invar(3))',
           `RECALC_TAB("spsu2/spsu2-2ch-tripletm.dat", SPSU2::LENGTH_Tpm_2CH, Invar(3))',
           `RECALC_TAB("spsu2/spsu2-3ch-tripletm.dat", SPSU2::LENGTH_Tpm_3CH, Invar(3))');
   }
  } else {
   my_error("Not implemented.");
  }
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef QDIFF
#define QDIFF(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) \
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


void SymmetrySPSU2::recalc_global(DiagInfo &diag,
				  string name,
				  MatrixElements &cnew)
{
  // NOTE: none of these are implemented for substeps==true.

  if (name == "Qtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "spsu2/spsu2-1ch-Qtot.dat"
        break;
      case 2:
#include "spsu2/spsu2-2ch-Qtot.dat"
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
#include "spsu2/spsu2-2ch-qdiff.dat"
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
#include "spsu2/spsu2-2ch-q1.dat"
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
#include "spsu2/spsu2-2ch-q2.dat"
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
#include "spsu2/spsu2-1ch-Iztot.dat"
        break;
      case 2:
#include "spsu2/spsu2-2ch-Iztot.dat"
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
#include "spsu2/spsu2-1ch-Ixtot.dat"
        break;
      case 2:
#include "spsu2/spsu2-2ch-Ixtot.dat"
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
#include "spsu2/spsu2-1ch-Iptot.dat"
        break;
      case 2:
#include "spsu2/spsu2-2ch-Iptot.dat"
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
#include "spsu2/spsu2-1ch-Imtot.dat"
        break;
      case 2:
#include "spsu2/spsu2-2ch-Imtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

}
