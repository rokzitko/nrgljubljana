// *** WARNING!!! Modify nrg-recalc-DBLISOSZ.cc.m4, not nrg-recalc-DBLISOSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2010
// This file pertains to (I1,I2,Sz) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





template<typename SC>
MatrixElements_tmpl<SC> SymmetryDBLISOSZ_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii11 = I1.get("II1");
    Ispin ii21 = I1.get("II2");
    Sspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0m.dat" << ", Iop=" << Invar(2, 1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletm0m.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2, 1, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11 - 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0p.dat" << ", Iop=" << Invar(2, 1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletm0p.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2, 1, -1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0m.dat" << ", Iop=" << Invar(2, 1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletp0m.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2, 1, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0p.dat" << ", Iop=" << Invar(2, 1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doubletp0p.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2, 1, -1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0mm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mp.dat" << ", Iop=" << Invar(1, 2, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0mp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(1, 2, -1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0pm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pp.dat" << ", Iop=" << Invar(1, 2, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc> recalc_table = {
#include "dblisosz/dblisosz-2ch-doublet0pp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(1, 2, -1));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}

template<typename SC>
Opch_tmpl<SC> SymmetryDBLISOSZ_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    Ispin ii1p = Ip.get("II1");
    Ispin ii2p = Ip.get("II2");
    Sspin sszp = Ip.get("SSZ");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(ii1p + 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isoup-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p + 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isoup-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isoup-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isoup-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isodown-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isodown-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type1-isodown-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f> recalc_table = {
#include "dblisosz/dblisosz-2ch-type2-isodown-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
  }
  return opch;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

template<typename SC>
void SymmetryDBLISOSZ_tmpl<SC>::recalc_global(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const std::string name, MatrixElements_tmpl<SC> &cnew) {
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
