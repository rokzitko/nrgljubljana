#include <cmath>
#include <vector>

namespace NRG {

using std::abs;
using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,J) subspaces

// Recalculate matrix elements of a doublet tensor operator
template<typename SC>
MatrixElements<SC> SymmetryQJ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int jj1 = I1.get("JJ");
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-doubletp.dat"
    }; });

    Ip = Invar(q1 - 1, jj1 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-doubletm.dat"
    }; });
  }
  return cnew;
}

#undef If
#define If(cond, a, b) (cond ? a : b)

// Recalculate matrix elements of a quadruplet tensor operator
template<typename SC>
MatrixElements<SC> SymmetryQJ<SC>::recalc_quadruplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int jj1 = I1.get("JJ");
    //    double J = (jj1-1.0)/2.0;
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 3);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 4), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-quad1.dat"
    }; });

    Ip = Invar(q1 - 1, jj1 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 4), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-quad2.dat"
    }; });

    Ip = Invar(q1 - 1, jj1 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 4), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-quad3.dat"
    }; });

    Ip = Invar(q1 - 1, jj1 - 3);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 4), cnew, [&] { return RecalcTable<SC>{
#include "qj/qj-quad4.dat"
    }; });
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryQJ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int jjp = Ip.get("JJ");
    double j  = J(jjp);
    Invar I1;

    I1 = Invar(qp + 1, jjp + 3);
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j3_2-jz3_2.dat"
    }; });

    I1 = Invar(qp + 1, jjp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j1_2-jz1_2.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j3_2-jz1_2.dat"
    }; });

    I1 = Invar(qp + 1, jjp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j1_2-jz-1_2.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j3_2-jz-1_2.dat"
    }; });

    I1 = Invar(qp + 1, jjp - 3);
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qj/qj-spin_j3_2-jz-3_2.dat"
    }; });
  }
  return opch;
}

}
