// *** WARNING!!! Modify nrg-recalc-DBLISOSZ.cc.m4, not nrg-recalc-DBLISOSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2010
// This file pertains to (I1,I2,Sz) subspaces

include(recalc-macros.m4)

namespace DBLISOSZ {
#include "dblisosz/dblisosz-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
void SymmetryDBLISOSZ::recalc_doublet(DiagInfo &diag,
                    	 	    MatrixElements &cold,
                    	 	    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Ispin ii11 = I1.get("II1");
    Ispin ii21 = I1.get("II2");
    Sspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii11-1, ii21, ssz1-1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doubletm0m.dat", DBLISOSZ::LENGTH_D1_2CH, Invar(2, 0, -1));

    Ip = Invar(ii11-1, ii21, ssz1+1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doubletm0p.dat", DBLISOSZ::LENGTH_D1_2CH, Invar(2, 0, +1));

    Ip = Invar(ii11+1, ii21, ssz1-1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doubletp0m.dat", DBLISOSZ::LENGTH_D1_2CH, Invar(2, 0, -1));

    Ip = Invar(ii11+1, ii21, ssz1+1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doubletp0p.dat", DBLISOSZ::LENGTH_D1_2CH, Invar(2, 0, +1));

    Ip = Invar(ii11, ii21-1, ssz1-1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doublet0mm.dat", DBLISOSZ::LENGTH_D2_2CH, Invar(0, 2, -1));

    Ip = Invar(ii11, ii21-1, ssz1+1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doublet0mp.dat", DBLISOSZ::LENGTH_D2_2CH, Invar(0, 2, +1));

    Ip = Invar(ii11, ii21+1, ssz1-1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doublet0pm.dat", DBLISOSZ::LENGTH_D2_2CH, Invar(0, 2, -1));

    Ip = Invar(ii11, ii21+1, ssz1+1);
    RECALC_TAB("dblisosz/dblisosz-2ch-doublet0pp.dat", DBLISOSZ::LENGTH_D2_2CH, Invar(0, 2, +1));
  }
}


// Driver routine for recalc_f()
void SymmetryDBLISOSZ::recalc_irreduc(const DiagInfo &diag)
{
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;

    Ispin ii1p = Ip.get("II1");
    Ispin ii2p = Ip.get("II2");
    Sspin sszp = Ip.get("SSZ");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    // RECALC_F_TAB_... (filename, channel_number, array_length)

    I1 = Invar(ii1p+1, ii2p, sszp+1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type1-isoup-a.dat", 0, DBLISOSZ::LENGTH_I_2CH);
    
    I1 = Invar(ii1p+1, ii2p, sszp-1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type2-isoup-a.dat", 0, DBLISOSZ::LENGTH_I_2CH);

    I1 = Invar(ii1p, ii2p+1, sszp+1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type1-isoup-b.dat", 1, DBLISOSZ::LENGTH_I_2CH);

    I1 = Invar(ii1p, ii2p+1, sszp-1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type2-isoup-b.dat", 1, DBLISOSZ::LENGTH_I_2CH);
    
    I1 = Invar(ii1p-1, ii2p, sszp+1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type1-isodown-a.dat", 0, DBLISOSZ::LENGTH_I_2CH);

    I1 = Invar(ii1p-1, ii2p, sszp-1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type2-isodown-a.dat", 0, DBLISOSZ::LENGTH_I_2CH);

    I1 = Invar(ii1p, ii2p-1, sszp+1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type1-isodown-b.dat", 1, DBLISOSZ::LENGTH_I_2CH);

    I1 = Invar(ii1p, ii2p-1, sszp-1);
    RECALC_F_TAB("dblisosz/dblisosz-2ch-type2-isodown-b.dat", 1, DBLISOSZ::LENGTH_I_2CH);
  }
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetryDBLISOSZ::recalc_global(DiagInfo &diag,
                                string name,
                                MatrixElements &cnew)
{
   if (name == "SZtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 2:
#include "dblisosz/dblisosz-2ch-spinz.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

}
