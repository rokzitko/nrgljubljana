#include <cmath>
#include <vector>
#include <string>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, May 2010
// This file pertains to (SZ,P) subspaces

template<typename SC>
MatrixElements<SC> SymmetrySPU1LR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1 + 1, p1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(-1, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-1ch-doubletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(-1, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-2ch-doubletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ssz1-1, p1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(+1, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-1ch-doubletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(+1, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-2ch-doubletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySPU1LR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(sszp + 1, pp);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-1ch-spinupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spinupb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(sszp-1, pp);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-1ch-spindowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spindownb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

   // CASE II: OPPOSITE PARITY

    if (P.channels == 2) {
      I1 = Invar(sszp + 1, -pp);
      this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spinupdiffa.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spinupdiffb.dat"
      }; });

      I1 = Invar(sszp - 1, -pp);
      this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spindowndiffa.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spu1lr/spu1lr-2ch-spindowndiffb.dat"
      }; });
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

    Ip = Invar(ssz1, p1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-1ch-triplets.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-2ch-triplets.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ssz1+2, p1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(-2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-1ch-tripletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(-2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-2ch-tripletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ssz1-2, p1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(+2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-1ch-tripletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(+2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spu1lr/spu1lr-2ch-tripletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

template<typename SC>
void SymmetrySPU1LR<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
