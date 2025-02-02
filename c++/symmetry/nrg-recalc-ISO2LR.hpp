namespace NRG {

// *** WARNING!!! Modify nrg-recalc-ISOLR.cc.m4, not nrg-recalc-ISOLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Nov 2007, May 2008
// This file pertains to (I,S,P) subspaces
// Version for EVEN number of impurities

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


double sign(double x) {
  if (x > 0.0) return +1.0;
  if (x < 0.0) return -1.0;
  my_assert_not_reached();
}

// (ISOLR): 8 calls of recalc_f() are necessary: different parities are also possible!

template<typename SC>
Opch<SC> SymmetryISO2LR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels.

    int iip = Ip.get("II");
    int ssp = Ip.get("SS");

    // IMPORTANT NEW ELEMENT: the parity is important here!!
    int pp = Ip.get("P");

    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // Both parities yield non-zero <I+-1/2, S+-1/2, P| a^\mu_\nu
    // |I,S,P'>.  Coefficients *DO* depend on P,P', or more accurately,
    // on whether or not P and P' are the same.
    //
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip + 1, ssp + 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isoupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isoupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isoupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isoupb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip + 1, ssp - 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isoupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isoupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isoupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isoupb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip - 1, ssp + 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isodowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isodowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isodownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isodownb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip - 1, ssp - 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isodowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isodowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isodownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isodownb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip + 1, ssp + 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isoupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isoupdiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isoupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isoupdiffb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip + 1, ssp - 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isoupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isoupdiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isoupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isoupdiffb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip - 1, ssp + 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isodowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isodowndiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spinup-isodowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spinup-isodowndiffb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(iip - 1, ssp - 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isodowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isodowndiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "iso2lr/iso2lr-2ch-spindown-isodowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-spindown-isodowndiffb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISO2LR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1 - 1, ss1 + 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-doubletmp.dat" << ", Iop=" << Invar(2, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-doubletmp.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(ii1 - 1, ss1 - 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-doubletmm.dat" << ", Iop=" << Invar(2, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-doubletmm.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(ii1 + 1, ss1 + 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-doubletpp.dat" << ", Iop=" << Invar(2, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-doubletpp.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(ii1 + 1, ss1 - 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-doubletpm.dat" << ", Iop=" << Invar(2, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-doubletpm.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(2, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISO2LR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ss1, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-triplets.dat" << ", Iop=" << Invar(1, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-triplets.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(ii1, ss1 + 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-tripletp.dat" << ", Iop=" << Invar(1, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-tripletp.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(ii1, ss1 - 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "iso2lr/iso2lr-2ch-tripletm.dat" << ", Iop=" << Invar(1, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "iso2lr/iso2lr-2ch-tripletm.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};
  }
  return cnew;
}

}
