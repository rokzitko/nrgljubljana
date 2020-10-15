// *** WARNING!!! Modify nrg-recalc-DBLSU2.cc.m4, not nrg-recalc-DBLSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2009
// This file pertains to (I1,I2) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace DBLSU2 {
#include "dblsu2/dblsu2-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
MatrixElements SymmetryDBLSU2::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii11 = I1.get("II1");
    Ispin ii21 = I1.get("II2");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doubletm0.dat" << ", len=" << DBLSU2::LENGTH_D1_2CH << ", Iop=" << Invar(2, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doubletm0.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D1_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, DBLSU2::LENGTH_D1_2CH, Invar(2, 0));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(ii11 + 1, ii21);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doubletp0.dat" << ", len=" << DBLSU2::LENGTH_D1_2CH << ", Iop=" << Invar(2, 0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doubletp0.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D1_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, DBLSU2::LENGTH_D1_2CH, Invar(2, 0));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(ii11, ii21 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doublet0m.dat" << ", len=" << DBLSU2::LENGTH_D2_2CH << ", Iop=" << Invar(0, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doublet0m.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D2_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, DBLSU2::LENGTH_D2_2CH, Invar(0, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};

    Ip = Invar(ii11, ii21 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblsu2/dblsu2-2ch-doublet0p.dat" << ", len=" << DBLSU2::LENGTH_D2_2CH << ", Iop=" << Invar(0, 2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doublet0p.dat"
      };
      BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D2_2CH);
      cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, DBLSU2::LENGTH_D2_2CH, Invar(0, 2));
    } else {
      cnew[II] = Matrix(0,0); // ???
    }
  }
};
  }
  return cnew;
}

// Override the recalc_f definition: we need to track the type (1 or 2) of
// the f-matrices.



// Driver routine for recalc_f()
Opch SymmetryDBLSU2::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    Ispin ii1p = Ip.get("II1");
    Ispin ii2p = Ip.get("II2");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(ii1p + 1, ii2p);
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isoup-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isoup-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[0][1][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p + 1);
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isoup-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isoup-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[1][1][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p - 1, ii2p);
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isodown-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isodown-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[0][1][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p - 1);
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isodown-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};
    {
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isodown-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
    opch[1][1][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, DBLSU2::LENGTH_I_2CH);
  }
};
  }
  return opch;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#ifdef NRG_COMPLEX
#undef SPINY
#define SPINY(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef Complex
#define Complex(x, y) cmpl(x, y)
#endif // NRG_COMPLEX

void SymmetryDBLSU2::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

#ifdef NRG_COMPLEX
  if (name == "SYtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif // NRG_COMPLEX

  if (name == "SXtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
