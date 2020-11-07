namespace NRG {

// *** WARNING!!! Modify nrg-recalc-DBLISOSZ.cc.m4, not nrg-recalc-DBLISOSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2010
// This file pertains to (I1,I2,Sz) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
MatrixElements<SC> SymmetryDBLISOSZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii11 = I1.get("II1");
    int ii21 = I1.get("II2");
    int ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0m.dat" << ", Iop=" << Invar(2, 1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletm0m.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 1, +1));
    }
  }
};

    Ip = Invar(ii11 - 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0p.dat" << ", Iop=" << Invar(2, 1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletm0p.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 1, -1));
    }
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0m.dat" << ", Iop=" << Invar(2, 1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletp0m.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 1, +1));
    }
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0p.dat" << ", Iop=" << Invar(2, 1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletp0p.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 1, -1));
    }
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0mm.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
    }
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mp.dat" << ", Iop=" << Invar(1, 2, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0mp.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, -1));
    }
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0pm.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
    }
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pp.dat" << ", Iop=" << Invar(1, 2, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0pp.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, -1));
    }
  }
};
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryDBLISOSZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int ii1p = Ip.get("II1");
    int ii2p = Ip.get("II2");
    int sszp = Ip.get("SSZ");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(ii1p + 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isoup-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p + 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isoup-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isoup-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isoup-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isodown-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isodown-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isodown-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isodown-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryDBLISOSZ<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblisosz/dblisosz-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
