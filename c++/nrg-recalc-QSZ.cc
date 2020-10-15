// *** WARNING!!! Modify nrg-recalc-QSZ.cc.m4, not nrg-recalc-QSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, 2006-2020
// This file pertains to (Q,SZ) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





// NOTE: p is ket (right side), 1 is bra (left side). OP is sandwiched
// in between. Thus Q[p] + Q[op] = Q[1].

// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetryQSZ::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Number q1   = I1.get("Q");
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      // In the case of (Q,S_z) basis, spin up and spin down are not
      // equivalent. The distinction appears in this recalculation code, but
      // also during the evaluation of the spectral densities, where two
      // spin-diagonal spectral densities (and also spin-flip spectral
      // density) can be defined.

      Ip = Invar(q1 - 1, ssz1 + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletp.dat" << ", Iop=" << Invar(1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-doubletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-doubletp.dat" << ", Iop=" << Invar(1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-2ch-doubletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-doubletp.dat" << ", Iop=" << Invar(1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-3ch-doubletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(q1-1, ssz1-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletm.dat" << ", Iop=" << Invar(1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-doubletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-doubletm.dat" << ", Iop=" << Invar(1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-2ch-doubletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-doubletm.dat" << ", Iop=" << Invar(1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-3ch-doubletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };
    }      // loop
  } else { // substeps
    for(const auto &[I1, eig]: diag) {
      Number q1   = I1.get("Q");
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(q1 - 1, ssz1 + 1);
      {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletp.dat" << ", Iop=" << Invar(1, -1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-doubletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, -1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

      Ip = Invar(q1 - 1, ssz1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-doubletm.dat" << ", Iop=" << Invar(1, +1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-doubletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(1, +1));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
    } // loop
  }
  return cnew;
}

// Driver routine for recalc_f()
Opch SymmetryQSZ::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp   = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for
    // several channels.

    I1 = Invar(qp + 1, sszp + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-1ch-spinupa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-2ch-spinupa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
      	   {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-2ch-spinupb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spinupa.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spinupb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};
      	   {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spinupc.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spinupc.dat"
      };
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[2][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(qp+1, sszp-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-1ch-spindowna.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-2ch-spindowna.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
	         {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-2ch-spindownb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spindowna.dat"
      };
      opch[0][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[0][0][II] = Matrix(0,0); // ???
    }
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spindownb.dat"
      };
      opch[1][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[1][0][II] = Matrix(0,0); // ???
    }
  }
};
      	   {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-3ch-spindownc.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-3ch-spindownc.dat"
      };
      opch[2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[2][0][II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  } // loop
  return opch;
}

OpchChannel SymmetryQSZ::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp   = Ip.get("Q");
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spinupa.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-1ch-spinupa.dat"
      };
      opch[M][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[M][0][II] = Matrix(0,0); // ???
    }
  }
};

    I1 = Invar(qp + 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "qsz/qsz-1ch-spindowna.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc_f recalc_table[] = {
#include "qsz/qsz-1ch-spindowna.dat"
      };
      opch[M][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
    } else {
      opch[M][0][II] = Matrix(0,0); // ???
    }
  }
};
  } // loop
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetryQSZ::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      Number q1   = I1.get("Q");
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(q1, ssz1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-triplets.dat" << ", Iop=" << Invar(0, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-triplets.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 0));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-triplets.dat" << ", Iop=" << Invar(0, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-2ch-triplets.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 0));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-3ch-triplets.dat" << ", Iop=" << Invar(0, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-3ch-triplets.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, 0));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(q1, ssz1+2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-tripletp.dat" << ", Iop=" << Invar(0, -2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-tripletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, -2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-tripletp.dat" << ", Iop=" << Invar(0, -2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-2ch-tripletp.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, -2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(q1, ssz1-2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-1ch-tripletm.dat" << ", Iop=" << Invar(0, +2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-1ch-tripletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, +2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "qsz/qsz-2ch-tripletm.dat" << ", Iop=" << Invar(0, +2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "qsz/qsz-2ch-tripletm.dat"
      };
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), Invar(0, +2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
} } break;
  default: my_assert_not_reached();
  };
    }      // loop
  } else { // substeps
    my_assert_not_reached();
  }
  return cnew;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef Q1U
#define Q1U(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef Q1D
#define Q1D(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

void SymmetryQSZ::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-spinz.dat"
          break;
        case 2:
#include "qsz/qsz-2ch-spinz.dat"
          break;
        case 3:
#include "qsz/qsz-3ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1u") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-q1u.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1d") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "qsz/qsz-1ch-q1d.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
