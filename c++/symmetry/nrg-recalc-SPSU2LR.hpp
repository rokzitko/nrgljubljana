#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Sep 2015
// This file pertains to (S,P) subspaces

template<typename SC>
MatrixElements<SC> SymmetrySPSU2LR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1 + 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spsu2lr/spsu2lr-2ch-doubletp.dat"
    }; });

    Ip = Invar(ss1 - 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 1), cnew, [&] { return RecalcTable<SC>{
#include "spsu2lr/spsu2lr-2ch-doubletm.dat"
    }; });
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySPSU2LR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    int pp    = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(ssp + 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spinupa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spinupb.dat"
    }; });

    I1 = Invar(ssp - 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spindowna.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spindownb.dat"
    }; });

    // CASE II: OPPOSITE PARITY

    I1 = Invar(ssp + 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spinupdiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spinupdiffb.dat"
    }; });

    I1 = Invar(ssp - 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spindowndiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2lr/spsu2lr-2ch-spindowndiffb.dat"
    }; });
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetrySPSU2LR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3, 1), cnew, [&] { return RecalcTable<SC>{
#include "spsu2lr/spsu2lr-2ch-triplets.dat"
    }; });

    Ip = Invar(ss1 + 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3, 1), cnew, [&] { return RecalcTable<SC>{
#include "spsu2lr/spsu2lr-2ch-tripletp.dat"
    }; });

    Ip = Invar(ss1 - 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3, 1), cnew, [&] { return RecalcTable<SC>{
#include "spsu2lr/spsu2lr-2ch-tripletm.dat"
    }; });
  }
  return cnew;
}

}
