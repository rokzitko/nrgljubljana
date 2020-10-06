// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2008.
// This file pertains to (S) subspaces

namespace SPU1 {
#include "spu1/spu1-1ch-def.dat"
#include "spu1/spu1-2ch-def.dat"
} // namespace SPU1

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetrySPU1::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1 + 1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletp.dat" << ", len=" << SPU1::LENGTH_D_1CH << ", Iop=" << Invar(-1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_1CH, Invar(-1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-doubletp.dat" << ", len=" << SPU1::LENGTH_D_2CH << ", Iop=" << Invar(-1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-2ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_2CH, Invar(-1));
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1-1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletm.dat" << ", len=" << SPU1::LENGTH_D_1CH << ", Iop=" << Invar(+1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_1CH, Invar(+1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-doubletm.dat" << ", len=" << SPU1::LENGTH_D_2CH << ", Iop=" << Invar(+1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-2ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_2CH, Invar(+1));
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1 + 1);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletp.dat" << ", len=" << SPU1::LENGTH_D_1CH << ", Iop=" << Invar(-1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_1CH, Invar(-1));
  }
};

      Ip = Invar(ssz1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletm.dat" << ", len=" << SPU1::LENGTH_D_1CH << ", Iop=" << Invar(+1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_D_1CH, Invar(+1));
  }
};
    }
  }
  return cnew;
}

// Driver routine for recalc_f()
Opch SymmetrySPU1::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  my_assert(!P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPU1::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-1ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[0][0], I1, Ip, recalc_table, SPU1::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPU1::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-2ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], I1, Ip, recalc_table, SPU1::LENGTH_I_2CH);
  }
};
 	     {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << SPU1::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-2ch-spinupb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], I1, Ip, recalc_table, SPU1::LENGTH_I_2CH);
  }
} } break;
  default: my_assert_not_reached();
  };

     I1 = Invar(sszp-1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPU1::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-1ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[0][0], I1, Ip, recalc_table, SPU1::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPU1::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-2ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], I1, Ip, recalc_table, SPU1::LENGTH_I_2CH);
  }
};
             {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << SPU1::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-2ch-spindownb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], I1, Ip, recalc_table, SPU1::LENGTH_I_2CH);
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return opch;
}

// Driver routine for recalc_f()
OpchChannel SymmetrySPU1::recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) {
  my_assert(P.substeps);
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spinupa.dat" << ", ch=" << M << ", len=" << SPU1::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-1ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[M][0], I1, Ip, recalc_table, SPU1::LENGTH_I_1CH);
  }
};

    I1 = Invar(sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spindowna.dat" << ", ch=" << M << ", len=" << SPU1::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1/spu1-1ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_I_1CH);
    recalc_f(diag, qsrmax, opch[M][0], I1, Ip, recalc_table, SPU1::LENGTH_I_1CH);
  }
};
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetrySPU1::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-triplets.dat" << ", len=" << SPU1::LENGTH_T0_1CH << ", Iop=" << Invar(0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_T0_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_T0_1CH, Invar(0));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-triplets.dat" << ", len=" << SPU1::LENGTH_T0_2CH << ", Iop=" << Invar(0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-2ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_T0_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_T0_2CH, Invar(0));
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1+2);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletp.dat" << ", len=" << SPU1::LENGTH_Tpm_1CH << ", Iop=" << Invar(-2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_1CH, Invar(-2));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-tripletp.dat" << ", len=" << SPU1::LENGTH_Tpm_2CH << ", Iop=" << Invar(-2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-2ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_2CH, Invar(-2));
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1-2);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletm.dat" << ", len=" << SPU1::LENGTH_Tpm_1CH << ", Iop=" << Invar(+2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_1CH, Invar(+2));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-tripletm.dat" << ", len=" << SPU1::LENGTH_Tpm_2CH << ", Iop=" << Invar(+2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-2ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_2CH, Invar(+2));
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-triplets.dat" << ", len=" << SPU1::LENGTH_T0_1CH << ", Iop=" << Invar(0) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_T0_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_T0_1CH, Invar(0));
  }
};

      Ip = Invar(ssz1 + 2);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletp.dat" << ", len=" << SPU1::LENGTH_Tpm_1CH << ", Iop=" << Invar(-2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_1CH, Invar(-2));
  }
};

      Ip = Invar(ssz1 - 2);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletm.dat" << ", len=" << SPU1::LENGTH_Tpm_1CH << ", Iop=" << Invar(+2) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1/spu1-1ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1::LENGTH_Tpm_1CH, Invar(+2));
  }
};
    }
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

// 2-channel only

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1UP
#define Q1UP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2UP
#define Q2UP(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1DO
#define Q1DO(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2DO
#define Q2DO(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

void SymmetrySPU1::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 1:
#include "spu1/spu1-1ch-Qtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Qtot.dat"
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
#include "spu1/spu1-1ch-Iztot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Iztot.dat"
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
#include "spu1/spu1-1ch-Ixtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Ixtot.dat"
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
#include "spu1/spu1-1ch-Iptot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Iptot.dat"
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
#include "spu1/spu1-1ch-Imtot.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-Imtot.dat"
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
#include "spu1/spu1-2ch-qdiff.dat"
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
#include "spu1/spu1-2ch-q1.dat"
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
#include "spu1/spu1-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1UP") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spu1/spu1-2ch-q1up.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2UP") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spu1/spu1-2ch-q2up.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1DO") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spu1/spu1-2ch-q1do.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2DO") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "spu1/spu1-2ch-q2do.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

#undef Q2
