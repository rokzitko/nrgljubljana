#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  return recalc_doublet_masked(diag, substruct, cold, qsztz_transfer::all);
}

// The local state is a spectator: both components of bra - ket are conserved together.
template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_doublet_masked(const DiagInfo<SC> &diag, const SubspaceStructure &substruct,
                                                                 const MatrixElements<SC> &cold, unsigned mask) const {
  MatrixElements<SC> cnew;
  if (mask == qsztz_transfer::none) return cnew;
  for(const auto &[I1, eig]: diag) {
    const auto q1 = I1.get("Q"), ssz1 = I1.get("SZ"), tz1 = I1.get("TZ");
    if (mask & qsztz_transfer::bit(-1, +1)) {
      const Invar Ip(q1 - 1, ssz1 - (-1), tz1 - (+1));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -1, +1), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletp-1.dat"
      }; });
    }
    if (mask & qsztz_transfer::bit(+1, +1)) {
      const Invar Ip(q1 - 1, ssz1 - (+1), tz1 - (+1));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +1, +1), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletm-1.dat"
      }; });
    }
    if (mask & qsztz_transfer::bit(-1, 0)) {
      const Invar Ip(q1 - 1, ssz1 - (-1), tz1 - (0));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -1, 0), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletp0.dat"
      }; });
    }
    if (mask & qsztz_transfer::bit(+1, 0)) {
      const Invar Ip(q1 - 1, ssz1 - (+1), tz1 - (0));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +1, 0), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletm0.dat"
      }; });
    }
    if (mask & qsztz_transfer::bit(-1, -1)) {
      const Invar Ip(q1 - 1, ssz1 - (-1), tz1 - (-1));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -1, -1), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletp+1.dat"
      }; });
    }
    if (mask & qsztz_transfer::bit(+1, -1)) {
      const Invar Ip(q1 - 1, ssz1 - (+1), tz1 - (-1));
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +1, -1), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-doubletm+1.dat"
      }; });
    }
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

template<typename SC>
Opch<SC> SymmetryQSZTZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int sszp  = Ip.get("SZ");
    int tzp = Ip.get("TZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1, tzp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spinup+1.dat"
    }; });

    I1 = Invar(qp + 1, sszp + 1, tzp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spinup0.dat"
    }; });

    I1 = Invar(qp + 1, sszp + 1, tzp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spinup-1.dat"
    }; });

    I1 = Invar(qp + 1, sszp - 1, tzp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spindo+1.dat"
    }; });

    I1 = Invar(qp + 1, sszp - 1, tzp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spindo0.dat"
    }; });

    I1 = Invar(qp + 1, sszp - 1, tzp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qsztz/qsztz-spindo-1.dat"
    }; });
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1   = I1.get("Q");
    int ssz1  = I1.get("SZ");
    int tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ssz1, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 0, 0), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-triplets.dat"
    }; });

    Ip = Invar(q1, ssz1 + 2, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, -2, 0), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-tripletp.dat"
    }; });

    Ip = Invar(q1, ssz1 - 2, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, +2, 0), cnew, [&] { return RecalcTable<SC>{
#include "qsztz/qsztz-tripletm.dat"
    }; });
  }
  return cnew;
}

}
