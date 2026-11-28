#include <cmath>
#include <vector>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006
// This file pertains to (Q,Sz,P) subspaces

template<typename SC>
Opch<SC> SymmetryQSZLR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels.

    // ****** CASE I: SAME PARITY ******

    Invar I1;
    I1 = Invar(qp + 1, sszp + 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spinupa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spinupb.dat"
    }; });

    I1 = Invar(qp + 1, sszp - 1, pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spindowna.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spindownb.dat"
    }; });

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, sszp + 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spinupdiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spinupdiffb.dat"
    }; });

    I1 = Invar(qp + 1, sszp - 1, -pp);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spindowndiffa.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qszlr/qszlr-2ch-spindowndiffb.dat"
    }; });
  }
  return opch;
}

}
