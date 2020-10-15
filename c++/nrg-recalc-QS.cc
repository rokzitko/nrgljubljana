// *** WARNING!!! Modify nrg-recalc-QS.cc.m4, not nrg-recalc-QS.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (Q,S) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace QS {
#include "qs/qs-1ch-def.dat"
#include "qs/qs-2ch-def.dat"
#include "qs/qs-3ch-def.dat"
#include "qs/qs-4ch-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
MatrixElements SymmetryQS::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Number q1 = I1.get("Q");
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-doubletp.dat" << ", len=" << QS::LENGTH_D_1CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_1CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-2ch-doubletp.dat" << ", len=" << QS::LENGTH_D_2CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-2ch-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_2CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-3ch-doubletp.dat" << ", len=" << QS::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-3ch-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_3CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-4ch-doubletp.dat" << ", len=" << QS::LENGTH_D_4CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-4ch-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_4CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_4CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
};

    Ip = Invar(q1-1, ss1-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-doubletm.dat" << ", len=" << QS::LENGTH_D_1CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_1CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-2ch-doubletm.dat" << ", len=" << QS::LENGTH_D_2CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-2ch-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_2CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-3ch-doubletm.dat" << ", len=" << QS::LENGTH_D_3CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-3ch-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_3CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-4ch-doubletm.dat" << ", len=" << QS::LENGTH_D_4CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-4ch-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_4CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_4CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
};
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      Number q1 = I1.get("Q");
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
      {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-doubletp.dat" << ", len=" << QS::LENGTH_D_1CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-doubletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_1CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

      Ip = Invar(q1 - 1, ss1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-doubletm.dat" << ", len=" << QS::LENGTH_D_1CH << ", Iop=" << Invar(1, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-doubletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_D_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_D_1CH, Invar(1, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
    }
  }
  return cnew;
}

// (QS): Two calls of recalc_f() are necessary (for S+1/2 and S-1/2)
// for each channel.
// See Krishna-Murthy p. 1034, equation (B10).

// Driver routine for recalc_f()
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
Opch SymmetryQS::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    I1 = Invar(qp + 1, ssp + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-1ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_1CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-1ch-spinupa.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_1CH);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_1CH);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_2CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-2ch-spinupa.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_2CH);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_2CH);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_2CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-2ch-spinupb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_2CH);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_2CH);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_3CH_0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spinupa.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_0);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_0);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_3CH_1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spinupb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_1);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_1);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spinupc.dat" << ", ch=" << 2 << ", len=" << QS::LENGTH_I_3CH_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spinupc.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_2);
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_2);
    } else {
      opch[2][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spinupa.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_4CH_0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spinupa.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_0);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_0);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spinupb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_4CH_1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spinupb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_1);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_1);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spinupc.dat" << ", ch=" << 2 << ", len=" << QS::LENGTH_I_4CH_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spinupc.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_2);
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_2);
    } else {
      opch[2][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spinupd.dat" << ", ch=" << 3 << ", len=" << QS::LENGTH_I_4CH_3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spinupd.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_3);
      opch[3][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_3);
    } else {
      opch[3][0][II] = Matrix(0,0);
    }
  }
} } break;
  default: my_assert_not_reached();
};

    I1 = Invar(qp+1, ssp-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-1ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_1CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-1ch-spindowna.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_1CH);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_1CH);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_2CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-2ch-spindowna.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_2CH);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_2CH);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_2CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-2ch-spindownb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_2CH);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_2CH);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_3CH_0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spindowna.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_0);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_0);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_3CH_1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spindownb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_1);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_1);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-3ch-spindownc.dat" << ", ch=" << 2 << ", len=" << QS::LENGTH_I_3CH_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-3ch-spindownc.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_3CH_2);
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_3CH_2);
    } else {
      opch[2][0][II] = Matrix(0,0);
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spindowna.dat" << ", ch=" << 0 << ", len=" << QS::LENGTH_I_4CH_0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spindowna.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_0);
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_0);
    } else {
      opch[0][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spindownb.dat" << ", ch=" << 1 << ", len=" << QS::LENGTH_I_4CH_1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spindownb.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_1);
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_1);
    } else {
      opch[1][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spindownc.dat" << ", ch=" << 2 << ", len=" << QS::LENGTH_I_4CH_2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spindownc.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_2);
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_2);
    } else {
      opch[2][0][II] = Matrix(0,0);
    }
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-4ch-spindownd.dat" << ", ch=" << 3 << ", len=" << QS::LENGTH_I_4CH_3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-4ch-spindownd.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_4CH_3);
      opch[3][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_4CH_3);
    } else {
      opch[3][0][II] = Matrix(0,0);
    }
  }
} } break;
  default: my_assert_not_reached();
};
  }
  return opch;
}

// Driver routine for recalc_f()
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
OpchChannel SymmetryQS::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(qp + 1, ssp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-1ch-spinupa.dat" << ", ch=" << M << ", len=" << QS::LENGTH_I_1CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-1ch-spinupa.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_1CH);
      opch[M][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_1CH);
    } else {
      opch[M][0][II] = Matrix(0,0);
    }
  }
};

    I1 = Invar(qp + 1, ssp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qs/qs-1ch-spindowna.dat" << ", ch=" << M << ", len=" << QS::LENGTH_I_1CH << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
      struct Recalc_f recalc_table[] = {
#include "qs/qs-1ch-spindowna.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_I_1CH);
      opch[M][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QS::LENGTH_I_1CH);
    } else {
      opch[M][0][II] = Matrix(0,0);
    }
  }
};
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
MatrixElements SymmetryQS::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Number q1 = I1.get("Q");
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1, ss1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-triplets.dat" << ", len=" << QS::LENGTH_T0_1CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-triplets.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_T0_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_T0_1CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-2ch-triplets.dat" << ", len=" << QS::LENGTH_T0_2CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-2ch-triplets.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_T0_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_T0_2CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-3ch-triplets.dat" << ", len=" << QS::LENGTH_T0_3CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-3ch-triplets.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_T0_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_T0_3CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-4ch-triplets.dat" << ", len=" << QS::LENGTH_T0_4CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-4ch-triplets.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_T0_4CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_T0_4CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
};

    Ip = Invar(q1, ss1+2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-tripletp.dat" << ", len=" << QS::LENGTH_Tpm_1CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-tripletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_1CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-2ch-tripletp.dat" << ", len=" << QS::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-2ch-tripletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_2CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-3ch-tripletp.dat" << ", len=" << QS::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-3ch-tripletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_3CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-4ch-tripletp.dat" << ", len=" << QS::LENGTH_Tpm_4CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-4ch-tripletp.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_4CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_4CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
};

    Ip = Invar(q1, ss1-2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-1ch-tripletm.dat" << ", len=" << QS::LENGTH_Tpm_1CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-1ch-tripletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_1CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_1CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-2ch-tripletm.dat" << ", len=" << QS::LENGTH_Tpm_2CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-2ch-tripletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_2CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-3ch-tripletm.dat" << ", len=" << QS::LENGTH_Tpm_3CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-3ch-tripletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_3CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_3CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 4: { {
  nrglog('f', "RECALC(fn=" << "qs/qs-4ch-tripletm.dat" << ", len=" << QS::LENGTH_Tpm_4CH << ", Iop=" << Invar(0, 3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qs/qs-4ch-tripletm.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QS::LENGTH_Tpm_4CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, QS::LENGTH_Tpm_4CH, Invar(0, 3));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
};
    }
  } else my_assert_not_reached();
  return cnew;
}

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef QTOT
#define QTOT(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
void SymmetryQS::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-qdiff.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-q1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-q2.dat"
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
#include "qs/qs-2ch-qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

#undef Q2
