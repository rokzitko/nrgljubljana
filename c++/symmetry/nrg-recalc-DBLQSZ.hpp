namespace NRG {

// *** WARNING!!! Modify nrg-recalc-DBLQSZ.cc.m4, not nrg-recalc-DBLQSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2022
// This file pertains to (Q1,Q2,Sz) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
MatrixElements<SC> SymmetryDBLQSZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q11 = I1.get("Q1");
    int q21 = I1.get("Q2");
    int ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(q11 - 1, q21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblqsz/dblqsz-2ch-doubletm0p.dat" << ", Iop=" << Invar(1, 0, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-doubletm0p.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 0, -1));
    }
  }
};

    Ip = Invar(q11 - 1, q21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblqsz/dblqsz-2ch-doubletm0m.dat" << ", Iop=" << Invar(1, 0, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-doubletm0m.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 0, +1));
    }
  }
};

    Ip = Invar(q11, q21 - 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblqsz/dblqsz-2ch-doublet0mp.dat" << ", Iop=" << Invar(0, 1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-doublet0mp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(0, 1, -1));
    }
  }
};

    Ip = Invar(q11, q21 - 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblqsz/dblqsz-2ch-doublet0mm.dat" << ", Iop=" << Invar(0, 1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-doublet0mm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(0, 1, +1));
    }
  }
};
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryDBLQSZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int q1p = Ip.get("Q1");
    int q2p = Ip.get("Q2");
    int sszp = Ip.get("SSZ");

    I1 = Invar(q1p + 1, q2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblqsz/dblqsz-2ch-spinup-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-spinup-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(q1p + 1, q2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblqsz/dblqsz-2ch-spindown-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-spindown-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(q1p, q2p + 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblqsz/dblqsz-2ch-spinup-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-spinup-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(q1p, q2p + 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblqsz/dblqsz-2ch-spindown-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblqsz/dblqsz-2ch-spindown-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryDBLQSZ<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
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
