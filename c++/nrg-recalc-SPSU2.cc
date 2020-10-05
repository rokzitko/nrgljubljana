// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Dec 2007
// This file pertains to (S) subspaces

namespace SPSU2 {
#include "spsu2/spsu2-1ch-def.dat"
#include "spsu2/spsu2-2ch-def.dat"
#include "spsu2/spsu2-3ch-def.dat"
} // namespace SPSU2

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetrySPSU2::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!substeps) {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletp.dat" << ", len=" << SPSU2::LENGTH_D_1CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_1CH, Invar(2));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-doubletp.dat" << ", len=" << SPSU2::LENGTH_D_2CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-2ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_2CH, Invar(2));
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-doubletp.dat" << ", len=" << SPSU2::LENGTH_D_3CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-3ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_3CH, Invar(2));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1-1);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletm.dat" << ", len=" << SPSU2::LENGTH_D_1CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_1CH, Invar(2));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-doubletm.dat" << ", len=" << SPSU2::LENGTH_D_2CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-2ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_2CH, Invar(2));
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-doubletm.dat" << ", len=" << SPSU2::LENGTH_D_3CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-3ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_3CH, Invar(2));
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
      {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletp.dat" << ", len=" << SPSU2::LENGTH_D_1CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_1CH, Invar(2));
  }
};

      Ip = Invar(ss1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletm.dat" << ", len=" << SPSU2::LENGTH_D_1CH << ", Iop=" << Invar(2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_D_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_D_1CH, Invar(2));
  }
};
    }
  }
  return cnew;
}

// Driver routine for recalc_f()
Opch SymmetrySPSU2::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  my_assert(!substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-1ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-2ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_2CH);
  }
};
	   {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << SPSU2::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-2ch-spinupb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_2CH);
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_3CH_0 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_0);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_0);
  }
};
	   {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupb.dat" << ", ch=" << 1 << ", len=" << SPSU2::LENGTH_I_3CH_1 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spinupb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_1);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_1);
  }
};
	   {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupc.dat" << ", ch=" << 2 << ", len=" << SPSU2::LENGTH_I_3CH_2 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spinupc.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_2);
    recalc_f(diag, qsrmax, opch[2][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_2);
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(ssp-1);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-1ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-2ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_2CH);
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << SPSU2::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-2ch-spindownb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_2CH);
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPSU2::LENGTH_I_3CH_0 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_0);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_0);
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindownb.dat" << ", ch=" << 1 << ", len=" << SPSU2::LENGTH_I_3CH_1 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spindownb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_1);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_1);
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindownc.dat" << ", ch=" << 2 << ", len=" << SPSU2::LENGTH_I_3CH_2 << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-3ch-spindownc.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_3CH_2);
    recalc_f(diag, qsrmax, opch[2][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_3CH_2);
  }
} } break;
  default: my_assert_not_reached();
  };

    // Note: for 3ch cases, the lengths in the three channels are not the same!
    // The same thing occurs for all SU(2)_spin cases, for instance for symtype=QS.
    // RZ, oct 2015
  }
  return opch;
}

// Driver routine for recalc_f()
OpchChannel SymmetrySPSU2::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  my_assert(substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spinupa.dat" << ", ch=" << M << ", len=" << SPSU2::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-1ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[M][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_1CH);
  }
};

    I1 = Invar(ssp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spindowna.dat" << ", ch=" << M << ", len=" << SPSU2::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spsu2/spsu2-1ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[M][0], Ip, I1, recalc_table, SPSU2::LENGTH_I_1CH);
  }
};
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetrySPSU2::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!substeps) {
    for(const auto &[I1, eig]: diag) {
      Sspin ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-triplets.dat" << ", len=" << SPSU2::LENGTH_T0_1CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_T0_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_T0_1CH, Invar(3));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-triplets.dat" << ", len=" << SPSU2::LENGTH_T0_2CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-2ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_T0_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_T0_2CH, Invar(3));
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-triplets.dat" << ", len=" << SPSU2::LENGTH_T0_3CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-3ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_T0_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_T0_3CH, Invar(3));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1+2);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-tripletp.dat" << ", len=" << SPSU2::LENGTH_Tpm_1CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_Tpm_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_Tpm_1CH, Invar(3));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-tripletp.dat" << ", len=" << SPSU2::LENGTH_Tpm_2CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-2ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_Tpm_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_Tpm_2CH, Invar(3));
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-tripletp.dat" << ", len=" << SPSU2::LENGTH_Tpm_3CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-3ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_Tpm_3CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_Tpm_3CH, Invar(3));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1-2);
    switch (channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-tripletm.dat" << ", len=" << SPSU2::LENGTH_Tpm_1CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-1ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_Tpm_1CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_Tpm_1CH, Invar(3));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-tripletm.dat" << ", len=" << SPSU2::LENGTH_Tpm_2CH << ", Iop=" << Invar(3) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spsu2/spsu2-2ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2::LENGTH_Tpm_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, SPSU2::LENGTH_Tpm_2CH, Invar(3));
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else my_assert_not_reached();
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

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

void SymmetrySPSU2::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  // NOTE: none of these are implemented for substeps==true.

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spsu2/spsu2-1ch-Qtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spsu2/spsu2-2ch-qdiff.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spsu2/spsu2-2ch-q1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spsu2/spsu2-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iztot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spsu2/spsu2-1ch-Ixtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iptot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spsu2/spsu2-1ch-Imtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
