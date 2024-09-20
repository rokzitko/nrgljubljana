namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SPSU2T.cc.m4, not nrg-recalc-SPSU2T.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (S,T) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
MatrixElements<SC> SymmetrySPSU2T<SC>::recalc_doublet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double T   = t1; // trick!
    Invar Ip;

    // Two different lengths: D_3CH_a and D_3CH_b

    Ip = Invar(ss1 + 1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletp-1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletp-1.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};

    Ip = Invar(ss1 - 1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletm-1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletm-1.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};

    Ip = Invar(ss1 + 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletp0.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletp0.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};

    Ip = Invar(ss1 - 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletm0.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletm0.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};

    Ip = Invar(ss1 + 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletp+1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletp+1.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};

    Ip = Invar(ss1 - 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-doubletm+1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-doubletm+1.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(1, 2, 1));
    }
  }
};
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

template<typename SC>
Opch<SC> SymmetrySPSU2T<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp  = Ip.get("SS");
    int tp = Ip.get("T");
    double T   = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=SPSU2T.
    I1 = Invar(ssp + 1, tp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spinup+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spinup+1.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp + 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spinup0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spinup0.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp + 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spinup-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spinup-1.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp - 1, tp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spindo+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spindo+1.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp - 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spindo0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spindo0.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp - 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2t/spsu2t-spindo-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2t/spsu2t-spindo-1.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetrySPSU2T<SC>::recalc_triplet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double T   = t1; // trick!
    Invar Ip;

    Ip = Invar(ss1, t1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-triplets.dat" << ", Iop=" << Invar(3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-triplets.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3, 0));
    }
  }
};

    Ip = Invar(ss1 + 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-tripletp.dat" << ", Iop=" << Invar(3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-tripletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3, 0));
    }
  }
};

    Ip = Invar(ss1 - 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2t/spsu2t-tripletm.dat" << ", Iop=" << Invar(3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2t/spsu2t-tripletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3, 0));
    }
  }
};
  }
  return cnew;
}

}
