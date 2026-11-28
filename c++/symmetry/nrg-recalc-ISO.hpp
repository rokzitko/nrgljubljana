#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (I,S) subspaces

template<typename SC>
MatrixElements<SC> SymmetryISO<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1 - 1, ss1 + 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-doubletmp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-doubletmp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1-1, ss1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-doubletmm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-doubletmm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1+1, ss1+1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-doubletpp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-doubletpp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1+1, ss1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-doubletpm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 2), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-doubletpm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryISO<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels.

    int iip = Ip.get("II");
    int ssp = Ip.get("SS");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(iip + 1, ssp + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-1ch-spinup-isoupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spinup-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spinup-isoupb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isoupb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isoupc.dat"
        }; });;
        break;
      default: my_assert_not_reached();
    }
    
    I1 = Invar(iip+1, ssp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-1ch-spindown-isoupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spindown-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spindown-isoupb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isoupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isoupb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isoupc.dat"
        }; });;
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip-1, ssp+1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-1ch-spinup-isodowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spinup-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spinup-isodownb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isodownb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spinup-isodownc.dat"
        }; });;
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip-1, ssp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-1ch-spindown-isodowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spindown-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-2ch-spindown-isodownb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isodowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isodownb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "iso/iso-3ch-spindown-isodownc.dat"
        }; });;
        break;
      default: my_assert_not_reached();
    }
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryISO<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1, ss1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-triplets.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-triplets.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1, ss1+2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-tripletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-tripletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1, ss1-2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-1ch-tripletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 3), cnew, [&] { return RecalcTable<SC>{
#include "iso/iso-2ch-tripletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

}
