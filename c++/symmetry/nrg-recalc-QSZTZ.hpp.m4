#include <cmath>
#include <initializer_list>

namespace NRG {

using std::sqrt;

// Edit nrg-recalc-QSZTZ.hpp.m4, not the generated header.

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  return recalc_named_qsztz_doublet(diag, substruct, cold, qsztz_transfer::all);
}

// The local state is a spectator: both components of bra - ket are conserved together.
define(`QSZTZ_DOUBLET', `if (mask & qsztz_transfer::bit($1, $2)) {
      const Invar Ip(q1 - 1, ssz1 - ($1), tz1 - ($2));
      RECALC_TAB($3, Invar(1, $1, $2));
    }')dnl
template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_named_qsztz_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct,
                                                                 const MatrixElements<SC> &cold, qsztz_transfer::Mask mask) const {
  MatrixElements<SC> cnew;
  if (mask == qsztz_transfer::none) return cnew;
  for(const auto &[I1, eig]: diag) {
    const auto q1 = I1.get("Q"), ssz1 = I1.get("SZ"), tz1 = I1.get("TZ");
    QSZTZ_DOUBLET(-1, +1, "qsztz/qsztz-doubletp-1.dat")
    QSZTZ_DOUBLET(+1, +1, "qsztz/qsztz-doubletm-1.dat")
    QSZTZ_DOUBLET(-1,  0, "qsztz/qsztz-doubletp0.dat")
    QSZTZ_DOUBLET(+1,  0, "qsztz/qsztz-doubletm0.dat")
    QSZTZ_DOUBLET(-1, -1, "qsztz/qsztz-doubletp+1.dat")
    QSZTZ_DOUBLET(+1, -1, "qsztz/qsztz-doubletm+1.dat")
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
    RECALC_F_TAB("qsztz/qsztz-spinup+1.dat", 0);

    I1 = Invar(qp + 1, sszp + 1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spinup0.dat", 0);

    I1 = Invar(qp + 1, sszp + 1, tzp - 1);
    RECALC_F_TAB("qsztz/qsztz-spinup-1.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp + 1);
    RECALC_F_TAB("qsztz/qsztz-spindo+1.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spindo0.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp - 1);
    RECALC_F_TAB("qsztz/qsztz-spindo-1.dat", 0);
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
    RECALC_TAB("qsztz/qsztz-triplets.dat", Invar(0, 0, 0));

    Ip = Invar(q1, ssz1 + 2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletp.dat", Invar(0, -2, 0));

    Ip = Invar(q1, ssz1 - 2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletm.dat", Invar(0, +2, 0));
  }
  return cnew;
}

}
