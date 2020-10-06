// *** WARNING!!! Modify nrg-recalc-QSC3.cc.m4, not nrg-recalc-QSC3.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (Q,S,P) subspaces, P modulo 3

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace QSC3 {
#include "qsc3/qsc3-def.dat"
}

#define xRECALC_F_TAB(a, b, c) 0;

// Driver routine for recalc_f()
Opch SymmetryQSC3::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
#ifdef NRG_COMPLEX

  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin ssp = Ip.get("SS");
    int p     = Ip.get("P");

    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(qp + 1, ssp + 1, (p + 0) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup0-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup0-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup0-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup0-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, ssp - 1, (p + 0) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown0-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown0-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown0-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown0-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, ssp + 1, (p + 1) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup1-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup1-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup1-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup1-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, ssp - 1, (p + 1) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown1-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown1-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown1-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown1-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, ssp + 1, (p + 2) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup2-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup2-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spinup2-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spinup2-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};

    I1 = Invar(qp + 1, ssp - 1, (p + 2) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-a.dat" << ", ch=" << 0 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown2-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[0][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-b.dat" << ", ch=" << 1 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown2-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[1][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "qsc3/qsc3-spindown2-c.dat" << ", ch=" << 2 << ", len=" << QSC3::LENGTH_I_3CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "qsc3/qsc3-spindown2-c.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == QSC3::LENGTH_I_3CH);
    opch[2][0][Twoinvar(I1, Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, QSC3::LENGTH_I_3CH);
  }
};
#undef Power
#undef sqrt
  }
#endif
  return opch;
}
