// *** WARNING!!! Modify nrg-recalc-QST.cc.m4, not nrg-recalc-QST.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (Q,S,T) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryQST::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double T   = t1; // trick!
    double S   = (ss1 - 1.) / 2.;
    Invar Ip;

    // Two different lengths: D_3CH_a and D_3CH_b

    // Invar(1,2,1) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ss1 + 1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp-1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletp-1.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm-1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletm-1.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ss1 + 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp0.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletp0.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm0.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletm0.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ss1 + 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletp+1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletp+1.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-doubletm+1.dat" << ", Iop=" << Invar(1, 2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-doubletm+1.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, 1));
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

// Driver routine for recalc_f()
Opch SymmetryQST::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp  = Ip.get("Q");
    Sspin ssp  = Ip.get("SS");
    Tangmom tp = Ip.get("T");
    double T   = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=QST.
    I1 = Invar(qp + 1, ssp + 1, tp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spinup+1.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, ssp + 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spinup0.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, ssp + 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spinup-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spinup-1.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo+1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spindo+1.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo0.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spindo0.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, tp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qst/qst-spindo-1.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "qst/qst-spindo-1.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
  }
  return opch;
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetryQST::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    Ip = Invar(q1, ss1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-triplets.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-triplets.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ss1 + 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-tripletp.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-tripletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ss1 - 2, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-tripletm.dat" << ", Iop=" << Invar(0, 3, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-tripletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, 0));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}

// Recalculate matrix elements of an orbital triplet tenzor operator
MatrixElements SymmetryQST::recalc_orb_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1  = I1.get("Q");
    Sspin ss1  = I1.get("SS");
    Tangmom t1 = I1.get("T");
    double S   = (ss1 - 1.) / 2.;
    double T   = t1; // trick!
    Invar Ip;

    // 0 = chargeless
    // 1 = spin singlet (deg=2S+1=1)
    // 1 = spin triplet (T=1)

    Ip = Invar(q1, ss1, t1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-triplets.dat" << ", Iop=" << Invar(0, 1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-orb-triplets.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ss1, t1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-tripletp.dat" << ", Iop=" << Invar(0, 1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-orb-tripletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1, 1));
      }(); // immediately executed lambda
    }
  }
};

    Ip = Invar(q1, ss1, t1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "qst/qst-orb-tripletm.dat" << ", Iop=" << Invar(0, 1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "qst/qst-orb-tripletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1, 1));
      }(); // immediately executed lambda
    }
  }
};
  }
  return cnew;
}
