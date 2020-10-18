// *** WARNING!!! Modify nrg-recalc-SL.cc.m4, not nrg-recalc-SL.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007
// This file pertains to the spinless-fermions code.

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





MatrixElements SymmetrySL::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Invar Ip  = Invar(q1 - 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "sl/sl-1ch-doublet.dat" << ", Iop=" << Invar(1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "sl/sl-1ch-doublet.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "sl/sl-2ch-doublet.dat" << ", Iop=" << Invar(1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "sl/sl-2ch-doublet.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "sl/sl-3ch-doublet.dat" << ", Iop=" << Invar(1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include "sl/sl-3ch-doublet.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

Opch SymmetrySL::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Invar I1  = Invar(qp + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-1ch-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-1ch-a.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-2ch-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-2ch-a.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
}; 
          {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-2ch-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-2ch-b.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-3ch-a.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
}; 
          {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-3ch-b.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
          {
  nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include "sl/sl-3ch-c.dat"
        };
        opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
}  } break;
  default: my_assert_not_reached();
  };
  }
  return opch;
}

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef QTOT
#define QTOT(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef N1
#define N1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef N2
#define N2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef N3
#define N3(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

void SymmetrySL::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "sl/sl-2ch-qdiff.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "sl/sl-2ch-qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N3") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N3.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
