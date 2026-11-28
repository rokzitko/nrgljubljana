#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (Q,S,T) subspaces

// Recalculate matrix elements of a doublet tensor operator
template<typename SC>
MatrixElements<SC> SymmetryQST<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1  = I1.get("Q");
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double T   = t1; // trick!
    double S   = (ss1 - 1.) / 2.;
    Invar Ip;

    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,1) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ss1 + 1, t1 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletp-1.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 - 1, t1 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletm-1.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 + 1, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletp0.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 - 1, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletm0.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 + 1, t1 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletp+1.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 - 1, t1 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-doubletm+1.dat"
    }; });
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

// Driver routine for recalc_f()
template<typename SC>
Opch<SC> SymmetryQST<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp  = Ip.get("Q");
    int ssp  = Ip.get("SS");
    int tp = Ip.get("T");
    double T   = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    I1 = Invar(qp + 1, ssp + 1, tp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spinup+1.dat"
    }; });

    I1 = Invar(qp + 1, ssp + 1, tp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spinup0.dat"
    }; });

    I1 = Invar(qp + 1, ssp + 1, tp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spinup-1.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spindo+1.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spindo0.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, tp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qst/qst-spindo-1.dat"
    }; });
  }
  return opch;
}

// Recalculate matrix elements of a triplet tenzor operator
template<typename SC>
MatrixElements<SC> SymmetryQST<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1  = I1.get("Q");
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    Ip = Invar(q1, ss1, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-triplets.dat"
    }; });

    Ip = Invar(q1, ss1 + 2, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-tripletp.dat"
    }; });

    Ip = Invar(q1, ss1 - 2, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, 0), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-tripletm.dat"
    }; });
  }
  return cnew;
}

// Recalculate matrix elements of an orbital triplet tenzor operator
template<typename SC>
MatrixElements<SC> SymmetryQST<SC>::recalc_orb_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1  = I1.get("Q");
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    // 0 = chargeless
    // 1 = spin singlet (deg=2S+1=1)
    // 1 = spin triplet (T=1)

    Ip = Invar(q1, ss1, t1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 1, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-orb-triplets.dat"
    }; });

    Ip = Invar(q1, ss1, t1 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 1, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-orb-tripletp.dat"
    }; });

    Ip = Invar(q1, ss1, t1 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 1, 1), cnew, [&] { return RecalcTable<SC>{
#include "qst/qst-orb-tripletm.dat"
    }; });
  }
  return cnew;
}

}
