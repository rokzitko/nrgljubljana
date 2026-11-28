#include <cmath>
#include <vector>
#include <string>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Jan 2009
// This file pertains to (I,S) subspaces

template<typename SC>
MatrixElements<SC> SymmetryISOSZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1   = I1.get("II");
    int ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii1 - 1, ssz1 + 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-doubletmp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-doubletmp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1-1, ssz1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-doubletmm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-doubletmm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1+1, ssz1+1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-doubletpp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-doubletpp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1+1, ssz1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-doubletpm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-doubletpm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

// (ISOSZ): Four calls of recalc_f() are necessary for each channel.

template<typename SC>
Opch<SC> SymmetryISOSZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels.

    int iip   = Ip.get("II");
    int sszp = Ip.get("SSZ");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(iip + 1, sszp + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-1ch-spinup-isoupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spinup-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spinup-isoupb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip+1, sszp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-1ch-spindown-isoupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spindown-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spindown-isoupb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip-1, sszp+1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-1ch-spinup-isodowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spinup-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spinup-isodownb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip-1, sszp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-1ch-spindown-isodowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spindown-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isosz/isosz-2ch-spindown-isodownb.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryISOSZ<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1   = I1.get("II");
    int ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii1, ssz1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 0), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-triplets.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 0), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-triplets.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1, ssz1+2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -2), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-tripletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -2), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-tripletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1, ssz1-2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +2), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-1ch-tripletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +2), cnew, [&] { return RecalcTable<SC>{
#include "isosz/isosz-2ch-tripletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

template<typename SC>
void SymmetryISOSZ<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 1:
#include "isosz/isosz-1ch-spinz.dat"
          break;
        case 2:
#include "isosz/isosz-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
