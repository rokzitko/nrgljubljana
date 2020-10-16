// *** WARNING!!! Modify nrg-recalc-SPSU2LR.cc.m4, not nrg-recalc-SPSU2LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Sep 2015
// This file pertains to (S,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





MatrixElements SymmetrySPSU2LR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1 + 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-doubletp.dat" << ", Iop=" << Invar(-1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-doubletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-1, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ss1 - 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-doubletm.dat" << ", Iop=" << Invar(+1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-doubletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+1, 1));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}

Opch SymmetrySPSU2LR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    int pp    = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(ssp + 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spinupa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spinupb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(ssp - 1, pp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spindowna.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spindownb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    // CASE II: OPPOSITE PARITY

    I1 = Invar(ssp + 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spinupdiffa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spinupdiffb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(ssp - 1, -pp);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spindowndiffa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2lr/spsu2lr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-spindowndiffb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
  }
  return opch;
}

MatrixElements SymmetrySPSU2LR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ss1, p1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-triplets.dat" << ", Iop=" << Invar(0, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-triplets.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ss1 + 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-tripletp.dat" << ", Iop=" << Invar(-2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-tripletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(ss1 - 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "spsu2lr/spsu2lr-2ch-tripletm.dat" << ", Iop=" << Invar(+2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spsu2lr/spsu2lr-2ch-tripletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+2, 1));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}
