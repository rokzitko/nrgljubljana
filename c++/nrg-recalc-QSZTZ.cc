// *** WARNING!!! Modify nrg-recalc-QSZTZ.cc.m4, not nrg-recalc-QSZTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





MatrixElements SymmetryQSZTZ::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1   = I1.get("Q");
    Sspin ssz1  = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp-1.dat" << ", Iop=" << Invar(1, -1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp-1.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm-1.dat" << ", Iop=" << Invar(1, +1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm-1.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1, +1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ssz1 + 1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp0.dat" << ", Iop=" << Invar(1, -1, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp0.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm0.dat" << ", Iop=" << Invar(1, +1, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm0.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletp+1.dat" << ", Iop=" << Invar(1, -1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletp+1.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1, -1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-doubletm+1.dat" << ", Iop=" << Invar(1, +1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-doubletm+1.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1, -1));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

Opch SymmetryQSZTZ::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp   = Ip.get("Q");
    Sspin sszp  = Ip.get("SZ");
    Tangmom tzp = Ip.get("TZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1, tzp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup+1.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, sszp + 1, tzp);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup0.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, sszp + 1, tzp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spinup-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spinup-1.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo+1.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo0.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, sszp - 1, tzp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsztz/qsztz-spindo-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qsztz/qsztz-spindo-1.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
  }
  return opch;
}

MatrixElements SymmetryQSZTZ::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1   = I1.get("Q");
    Sspin ssz1  = I1.get("SZ");
    Tangmom tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ssz1, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-triplets.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-triplets.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ssz1 + 2, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-tripletp.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-tripletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ssz1 - 2, tz1);
    {
  nrglog('f', "RECALC(fn=" << "qsztz/qsztz-tripletm.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qsztz/qsztz-tripletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}
