namespace NRG {

// *** WARNING!!! Modify nrg-recalc-DBLSU2.cc.m4, not nrg-recalc-DBLSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2009
// This file pertains to (I1,I2) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
MatrixElements<SC> SymmetryDBLSU2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii11 = I1.get("II1");
    int ii21 = I1.get("II2");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doubletm0.dat" << ", Iop=" << Invar(2, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-doubletm0.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 0));
    }
  }
};

    Ip = Invar(ii11 + 1, ii21);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doubletp0.dat" << ", Iop=" << Invar(2, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-doubletp0.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 0));
    }
  }
};

    Ip = Invar(ii11, ii21 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doublet0m.dat" << ", Iop=" << Invar(0, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-doublet0m.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(0, 2));
    }
  }
};

    Ip = Invar(ii11, ii21 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doublet0p.dat" << ", Iop=" << Invar(0, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-doublet0p.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(0, 2));
    }
  }
};
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryDBLSU2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int ii1p = Ip.get("II1");
    int ii2p = Ip.get("II2");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(ii1p + 1, ii2p);
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type1-isoup-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type1-isoup-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type2-isoup-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type2-isoup-a.dat"
      };
      opch[0][1][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type1-isoup-b.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type1-isoup-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type2-isoup-b.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type2-isoup-b.dat"
      };
      opch[1][1][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p - 1, ii2p);
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type1-isodown-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type1-isodown-a.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type2-isodown-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type2-isodown-a.dat"
      };
      opch[0][1][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ii1p, ii2p - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type1-isodown-b.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type1-isodown-b.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "dblsu2/dblsu2-2ch-type2-isodown-b.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "dblsu2/dblsu2-2ch-type2-isodown-b.dat"
      };
      opch[1][1][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

#undef SPINY
#define SPINY(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)
#undef Complex
#define Complex(x, y) cmpl(x, y)

template<typename SC>
void SymmetryDBLSU2<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "SYtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "SXtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
