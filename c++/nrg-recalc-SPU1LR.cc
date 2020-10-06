// *** WARNING!!! Modify nrg-recalc-SPU1LR.cc.m4, not nrg-recalc-SPU1LR.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, May 2010
// This file pertains to (SZ,P) subspaces

namespace SPU1LR {
#include "spu1lr/spu1lr-1ch-def.dat"
#include "spu1lr/spu1lr-2ch-def.dat"
} // namespace SPU1LR

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





// Recalculate matrix elements of a doublet tensor operator
MatrixElements SymmetrySPU1LR::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1 + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-doubletp.dat" << ", len=" << SPU1LR::LENGTH_D_1CH << ", Iop=" << Invar(-1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_D_1CH, Invar(-1, 1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-doubletp.dat" << ", len=" << SPU1LR::LENGTH_D_2CH << ", Iop=" << Invar(-1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-doubletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_D_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_D_2CH, Invar(-1, 1));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-doubletm.dat" << ", len=" << SPU1LR::LENGTH_D_1CH << ", Iop=" << Invar(+1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_D_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_D_1CH, Invar(+1, 1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-doubletm.dat" << ", len=" << SPU1LR::LENGTH_D_2CH << ", Iop=" << Invar(+1, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-doubletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_D_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_D_2CH, Invar(+1, 1));
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

// Driver routine for recalc_f()
Opch SymmetrySPU1LR::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    int pp      = Ip.get("P");
    Invar I1;

    // CASE I: SAME PARITY

    I1 = Invar(sszp + 1, pp);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-1ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-1ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_1CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupa.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};
	    {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupb.dat" << ", ch=" << 1 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(sszp-1, pp);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-1ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_1CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-1ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_1CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_1CH);
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowna.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowna.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindownb.dat" << ", ch=" << 1 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindownb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
} } break;
  default: my_assert_not_reached();
  };

   // CASE II: OPPOSITE PARITY

    if (channels == 2) {
      I1 = Invar(sszp + 1, -pp);
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupdiffa.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupdiffa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spinupdiffb.dat" << ", ch=" << 1 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spinupdiffb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};

      I1 = Invar(sszp - 1, -pp);
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowndiffa.dat" << ", ch=" << 0 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowndiffa.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[0][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "spu1lr/spu1lr-2ch-spindowndiffb.dat" << ", ch=" << 1 << ", len=" << SPU1LR::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "spu1lr/spu1lr-2ch-spindowndiffb.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_I_2CH);
    opch[1][0][Twoinvar(I1,Ip)] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, SPU1LR::LENGTH_I_2CH);
  }
};
    }
  }
  return opch;
}

// Recalculate matrix elements of a triplet tenzor operator
MatrixElements SymmetrySPU1LR::recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    SZspin ssz1 = I1.get("SSZ");
    int p1      = I1.get("P");
    Invar Ip;

    Ip = Invar(ssz1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-triplets.dat" << ", len=" << SPU1LR::LENGTH_T0_1CH << ", Iop=" << Invar(0, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_T0_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_T0_1CH, Invar(0, 1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-triplets.dat" << ", len=" << SPU1LR::LENGTH_T0_2CH << ", Iop=" << Invar(0, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-triplets.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_T0_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_T0_2CH, Invar(0, 1));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1+2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-tripletp.dat" << ", len=" << SPU1LR::LENGTH_Tpm_1CH << ", Iop=" << Invar(-2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_Tpm_1CH, Invar(-2, 1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-tripletp.dat" << ", len=" << SPU1LR::LENGTH_Tpm_2CH << ", Iop=" << Invar(-2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-tripletp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_Tpm_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_Tpm_2CH, Invar(-2, 1));
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ssz1-2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-1ch-tripletm.dat" << ", len=" << SPU1LR::LENGTH_Tpm_1CH << ", Iop=" << Invar(+2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-1ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_Tpm_1CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_Tpm_1CH, Invar(+2, 1));
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1lr/spu1lr-2ch-tripletm.dat" << ", len=" << SPU1LR::LENGTH_Tpm_2CH << ", Iop=" << Invar(+2, 1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "spu1lr/spu1lr-2ch-tripletm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPU1LR::LENGTH_Tpm_2CH);
    cnew[Twoinvar(I1, Ip)] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, SPU1LR::LENGTH_Tpm_2CH, Invar(+2, 1));
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
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
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
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
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
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
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
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
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
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
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
