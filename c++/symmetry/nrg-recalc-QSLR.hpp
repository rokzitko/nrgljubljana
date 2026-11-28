#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

template<typename SC>
Opch<SC> SymmetryQSLR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    // IMPORTANT NEW ELEMENT: the parity is important here!!
    int lrp = Ip.get("P");

    // Both parities yield non-zero <Q+1, S+-1/2, P| a^\dag_\nu
    // |Q,S,P'>.  Coefficients *DO* depend on P,P', or more
    // accurately, on whether or not P and P' are the same.
    //
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(qp + 1, ssp + 1, lrp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spinupa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spinupb.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, lrp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spindowna.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spindownb.dat"
    }; });

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, ssp + 1, -lrp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spinupdiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spinupdiffb.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1, -lrp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spindowndiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qslr/qslr-2ch-spindowndiffb.dat"
    }; });
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryQSLR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1 - 1, ss1 + 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, +1), cnew, [&] { return RecalcTable<SC>{
#include "qslr/qslr-2ch-doubletp.dat"
    }; });

    Ip = Invar(q1 - 1, ss1 - 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2, +1), cnew, [&] { return RecalcTable<SC>{
#include "qslr/qslr-2ch-doubletm.dat"
    }; });
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryQSLR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, +1), cnew, [&] { return RecalcTable<SC>{
#include "qslr/qslr-2ch-triplets.dat"
    }; });

    Ip = Invar(q1, ss1 + 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, +1), cnew, [&] { return RecalcTable<SC>{
#include "qslr/qslr-2ch-tripletp.dat"
    }; });

    Ip = Invar(q1, ss1 - 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3, +1), cnew, [&] { return RecalcTable<SC>{
#include "qslr/qslr-2ch-tripletm.dat"
    }; });
  }
  return cnew;
}

}
