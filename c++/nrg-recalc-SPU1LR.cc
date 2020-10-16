// *** WARNING!!! Modify nrg-recalc-SPU1LR.cc.m4, not nrg-recalc-SPU1LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, May 2010
// This file pertains to (SZ,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





MatrixElements SymmetrySPU1LR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1 + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-doubletp.dat" << ", Iop=" << Invar(-1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-doubletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-1, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-doubletp.dat" << ", Iop=" << Invar(-1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-doubletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-1, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-doubletm.dat" << ", Iop=" << Invar(+1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-doubletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+1, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-doubletm.dat" << ", Iop=" << Invar(+1, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-doubletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+1, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

Opch SymmetrySPU1LR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(sszp + 1, pp);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-1ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-1ch-spinupa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(sszp-1, pp);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-1ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-1ch-spindowna.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowna.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindownb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };

   // CASE II: OPPOSITE PARITY

    if (P.channels == 2) {
      I1 = Invar(sszp + 1, -pp);
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupdiffa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupdiffb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};

      I1 = Invar(sszp - 1, -pp);
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowndiffa.dat"
        };
        opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowndiffb.dat"
        };
        opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
};
    }
  }
  return opch;
}

MatrixElements SymmetrySPU1LR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-triplets.dat" << ", Iop=" << Invar(0, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-triplets.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-triplets.dat" << ", Iop=" << Invar(0, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-triplets.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1+2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-tripletp.dat" << ", Iop=" << Invar(-2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-tripletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-2, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-tripletp.dat" << ", Iop=" << Invar(-2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-tripletp.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(-2, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1-2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-tripletm.dat" << ", Iop=" << Invar(+2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-tripletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+2, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-tripletm.dat" << ", Iop=" << Invar(+2, 1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-tripletm.dat"
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(+2, 1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

void SymmetrySPU1LR::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Qtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Iztot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Ixtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Iptot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spu1lr/spu1lr-1ch-Imtot.dat"
          break;
        case 2:
#include "spu1lr/spu1lr-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
