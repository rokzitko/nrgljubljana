namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SPU1LR.cc.m4, not nrg-recalc-SPU1LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, May 2010
// This file pertains to (SZ,P) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetrySPU1LR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1 + 1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-doubletp.dat", Invar(-1, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-doubletp.dat", Invar(-1, 1))');

    Ip = Invar(ssz1-1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-doubletm.dat", Invar(+1, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-doubletm.dat", Invar(+1, 1))');
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySPU1LR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    int sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(sszp + 1, pp);
    ONETWO(`RECALC_F_TAB("spu1lr/spu1lr-1ch-spinupa.dat", 0)',

           `RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupa.dat", 0);
	          RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupb.dat", 1)');

    I1 = Invar(sszp-1, pp);
    ONETWO(`RECALC_F_TAB("spu1lr/spu1lr-1ch-spindowna.dat", 0)',

           `RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowna.dat", 0);
            RECALC_F_TAB("spu1lr/spu1lr-2ch-spindownb.dat", 1)');

   // CASE II: OPPOSITE PARITY

    if (P.channels == 2) {
      I1 = Invar(sszp + 1, -pp);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupdiffa.dat", 0);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spinupdiffb.dat", 1);

      I1 = Invar(sszp - 1, -pp);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowndiffa.dat", 0);
      RECALC_F_TAB("spu1lr/spu1lr-2ch-spindowndiffb.dat", 1);
    }
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetrySPU1LR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-triplets.dat", Invar(0, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-triplets.dat", Invar(0, 1))');

    Ip = Invar(ssz1+2);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-tripletp.dat", Invar(-2, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-tripletp.dat", Invar(-2, 1))');

    Ip = Invar(ssz1-2);
    ONETWO(`RECALC_TAB("spu1lr/spu1lr-1ch-tripletm.dat", Invar(+2, 1))',
           `RECALC_TAB("spu1lr/spu1lr-2ch-tripletm.dat", Invar(+2, 1))');
  }
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

template<typename SC>
void SymmetrySPU1LR<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Qtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Iztot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Ixtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Iptot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Imtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
