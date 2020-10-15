// *** WARNING!!! Modify nrg-recalc-QSLR.cc.m4, not nrg-recalc-QSLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Aug 2006
// This file pertains to (Q,S,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





// Driver routine for recalc_f()
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO
Opch SymmetryQSLR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
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
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spinupa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spinupb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spindowna.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spindownb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(qp + 1, ssp + 1, -lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spinupdiffa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spinupdiffb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1, -lrp);
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spindowndiffa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qslr/qslr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qslr/qslr-2ch-spindowndiffb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO
MatrixElements SymmetryQSLR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1 - 1, ss1 + 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletp.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qslr/qslr-2ch-doubletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1 - 1, ss1 - 1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-doubletm.dat" << ", Iop=" << Invar(1, 2, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qslr/qslr-2ch-doubletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, 2, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tenzor operator [EVEN PARITY]
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO
MatrixElements SymmetryQSLR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(q1, ss1, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-triplets.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qslr/qslr-2ch-triplets.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1, ss1 + 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletp.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qslr/qslr-2ch-tripletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(q1, ss1 - 2, p1);
    {
  nrglog('f', "RECALC(fn=" << "qslr/qslr-2ch-tripletm.dat" << ", Iop=" << Invar(0, 3, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qslr/qslr-2ch-tripletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 3, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
  }
  return cnew;
}
