namespace NRG {

// *** WARNING!!! Modify nrg-recalc-QJ.cc.m4, not nrg-recalc-QJ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,J) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


// Recalculate matrix elements of a doublet tensor operator
template<typename SC>
MatrixElements<SC> SymmetryQJ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int jj1 = I1.get("JJ");
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-doubletp.dat" << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-doubletp.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2));
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-doubletm.dat" << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-doubletm.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 2));
    }
  }
};
  }
  return cnew;
}

#undef If
#define If(cond, a, b) (cond ? a : b)

// Recalculate matrix elements of a quadruplet tensor operator
template<typename SC>
MatrixElements<SC> SymmetryQJ<SC>::recalc_quadruplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    int jj1 = I1.get("JJ");
    //    double J = (jj1-1.0)/2.0;
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 3);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad1.dat" << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-quad1.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 4));
    }
  }
};

    Ip = Invar(q1 - 1, jj1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad2.dat" << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-quad2.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 4));
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad3.dat" << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-quad3.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 4));
    }
  }
};

    Ip = Invar(q1 - 1, jj1 - 3);
    {
  nrglog('f', "RECALC(fn=" << "qj/qj-quad4.dat" << ", Iop=" << Invar(1, 4) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "qj/qj-quad4.dat"
      };
      cnew[II] = this->recalc_general(diag, substruct, cold, I1, Ip, recalc_table, Invar(1, 4));
    }
  }
};
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryQJ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int jjp = Ip.get("JJ");
    double j  = J(jjp);
    Invar I1;

    I1 = Invar(qp + 1, jjp + 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz3_2.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j3_2-jz3_2.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(qp + 1, jjp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz1_2.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j1_2-jz1_2.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz1_2.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j3_2-jz1_2.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(qp + 1, jjp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j1_2-jz-1_2.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j1_2-jz-1_2.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-1_2.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j3_2-jz-1_2.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(qp + 1, jjp - 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qj/qj-spin_j3_2-jz-3_2.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "qj/qj-spin_j3_2-jz-3_2.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, substruct, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch;
}

}
