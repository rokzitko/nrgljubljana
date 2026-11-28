#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2009
// This file pertains to (I,Sz,P) subspaces

template<typename SC>
Opch<SC> SymmetryISOSZLR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int iip   = Ip.get("II");
    int sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");

    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip + 1, sszp + 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isoupa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isoupb.dat"
    }; });

    I1 = Invar(iip + 1, sszp - 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isoupa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isoupb.dat"
    }; });

    I1 = Invar(iip - 1, sszp + 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isodowna.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isodownb.dat"
    }; });

    I1 = Invar(iip - 1, sszp - 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isodowna.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isodownb.dat"
    }; });

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip + 1, sszp + 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isoupdiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isoupdiffb.dat"
    }; });

    I1 = Invar(iip + 1, sszp - 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isoupdiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isoupdiffb.dat"
    }; });

    I1 = Invar(iip - 1, sszp + 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isodowndiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spinup-isodowndiffb.dat"
    }; });

    I1 = Invar(iip - 1, sszp - 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isodowndiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "isoszlr/isoszlr-2ch-spindown-isodowndiffb.dat"
    }; });
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISOSZLR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1   = I1.get("II");
    int ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1 - 1, ssz1 + 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-doubletmp.dat"
    }; });

    Ip = Invar(ii1 - 1, ssz1 - 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-doubletmm.dat"
    }; });

    Ip = Invar(ii1 + 1, ssz1 + 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, -1, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-doubletpp.dat"
    }; });

    Ip = Invar(ii1 + 1, ssz1 - 1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, +1, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-doubletpm.dat"
    }; });
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISOSZLR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1   = I1.get("II");
    int ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ssz1, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 0, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-triplets.dat"
    }; });

    Ip = Invar(ii1, ssz1 + 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, -2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-tripletp.dat"
    }; });

    Ip = Invar(ii1, ssz1 - 2, p1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, +2, +1), cnew, [&] { return RecalcTable<SC>{
#include "isoszlr/isoszlr-2ch-tripletm.dat"
    }; });
  }
  return cnew;
}

}
