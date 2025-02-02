namespace NRG {

// *** WARNING!!! Modify nrg-recalc-QSLR.cc.m4, not nrg-recalc-QSLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
Opch<SC> SymmetryQSLR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    // IMPORTANT NEW ELEMENT: the parity is important here!!
    int lrp = Ip.get("P");

    // Both parities yield non-zero <Q+1, S+-1/2, P| a^\dag_\nu
    // |Q,S,P'>.  Coefficients *DO* depend on P,P', or more
    // accurately, on whether or not P and P' are the same.
    //
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(qp + 1, ssp + 1, lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spinupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spinupb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spindowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spindownb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, ssp + 1, -lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spinupdiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spinupdiffb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, -lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spindowndiffa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qslr/qslr-2ch-spindowndiffb.dat"
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
MatrixElements<SC> SymmetryQSLR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1 - 1, ss1 + 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletp.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qslr/qslr-2ch-doubletp.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qslr/qslr-2ch-doubletm.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryQSLR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-triplets.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qslr/qslr-2ch-triplets.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(0, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(q1, ss1 + 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletp.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qslr/qslr-2ch-tripletp.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(0, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};

    Ip = Invar(q1, ss1 - 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletm.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qslr/qslr-2ch-tripletm.dat"
      };
      auto cn = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(0, 3, +1));
      if (cn) cnew[II] = *cn;
    }
  }
};
  }
  return cnew;
}

}
