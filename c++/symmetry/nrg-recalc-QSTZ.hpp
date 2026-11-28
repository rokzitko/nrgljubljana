#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2016
// This file pertains to (Q,S,Tz) subspaces

template<typename SC>
MatrixElements<SC> SymmetryQSTZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  return recalc_doublet_masked(diag, substruct, cold, qstz_transfer::all);
}

template<typename SC>
MatrixElements<SC> SymmetryQSTZ<SC>::recalc_doublet_masked(const DiagInfo<SC> &diag, const SubspaceStructure &substruct,
                                                          const MatrixElements<SC> &cold, unsigned mask) const {
  MatrixElements<SC> cnew;
  if (mask == qstz_transfer::none) return cnew;
  for(const auto &[I1, eig]: diag) {
    int q1   = I1.get("Q");
    int ss1   = I1.get("SS");
    int tz1 = I1.get("TZ");
    Invar Ip;

    // The added local state is a spectator: bra TZ - ket TZ is unchanged from the seed.
    if (mask & qstz_transfer::plus) {
      Ip = Invar(q1 - 1, ss1 + 1, tz1 - 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, +1), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletp-1.dat"
      }; });
      Ip = Invar(q1 - 1, ss1 - 1, tz1 - 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, +1), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletm-1.dat"
      }; });
    }
    if (mask & qstz_transfer::zero) {
      Ip = Invar(q1 - 1, ss1 + 1, tz1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 0), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletp0.dat"
      }; });
      Ip = Invar(q1 - 1, ss1 - 1, tz1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 0), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletm0.dat"
      }; });
    }
    if (mask & qstz_transfer::minus) {
      Ip = Invar(q1 - 1, ss1 + 1, tz1 + 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, -1), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletp+1.dat"
      }; });
      Ip = Invar(q1 - 1, ss1 - 1, tz1 + 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, -1), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-doubletm+1.dat"
      }; });
    }
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

template<typename SC>
Opch<SC> SymmetryQSTZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int ssp   = Ip.get("SS");
    int tzp = Ip.get("TZ");
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    I1 = Invar(qp + 1, ssp + 1, tzp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spinup+1.dat"
    }; });

    I1 = Invar(qp + 1, ssp + 1, tzp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spinup0.dat"
    }; });

    I1 = Invar(qp + 1, ssp + 1, tzp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spinup-1.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tzp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spindo+1.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tzp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spindo0.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tzp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qstz/qstz-spindo-1.dat"
    }; });
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryQSTZ<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1   = I1.get("Q");
    int ss1   = I1.get("SS");
    int tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ss1, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-triplets.dat"
    }; });

    Ip = Invar(q1, ss1 + 2, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-tripletp.dat"
    }; });

    Ip = Invar(q1, ss1 - 2, tz1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qstz/qstz-tripletm.dat"
    }; });
  }
  return cnew;
}

}
