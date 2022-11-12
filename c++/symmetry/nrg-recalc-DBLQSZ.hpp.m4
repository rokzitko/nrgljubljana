namespace NRG {

// *** WARNING!!! Modify nrg-recalc-DBLQSZ.cc.m4, not nrg-recalc-DBLQSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2022
// This file pertains to (Q1,Q2,Sz) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetryDBLQSZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q11 = I1.get("Q1");
    int q21 = I1.get("Q2");
    int ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q11 - 1, q21, ssz1 + 1);
    RECALC_TAB("dblqsz/dblqsz-2ch-doubletm0p.dat", Invar(1, 0, -1));

    Ip = Invar(q11 - 1, q21, ssz1 - 1);
    RECALC_TAB("dblqsz/dblqsz-2ch-doubletm0m.dat", Invar(1, 0, +1));

    Ip = Invar(q11, q21 - 1, ssz1 + 1);
    RECALC_TAB("dblqsz/dblqsz-2ch-doublet0mp.dat", Invar(0, 1, -1));

    Ip = Invar(q11, q21 - 1, ssz1 - 1);
    RECALC_TAB("dblqsz/dblqsz-2ch-doublet0mm.dat", Invar(0, 1, +1));
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryDBLQSZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int q1p = Ip.get("Q1");
    int q2p = Ip.get("Q2");
    int sszp = Ip.get("SSZ");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(q1p + 1, q2p, sszp + 1);
    RECALC_F_TAB("dblqsz/dblqsz-2ch-spinup-a.dat", 0);

    I1 = Invar(q1p + 1, q2p, sszp - 1);
    RECALC_F_TAB("dblqsz/dblqsz-2ch-spindown-a.dat", 0);

    I1 = Invar(q1p, q2p + 1, sszp + 1);
    RECALC_F_TAB("dblqsz/dblqsz-2ch-spinup-b.dat", 1);

    I1 = Invar(q1p, q2p + 1, sszp - 1);
    RECALC_F_TAB("dblqsz/dblqsz-2ch-spindown-b.dat", 1);
  }
  return opch;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryDBLQSZ<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblqsz/dblqsz-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
