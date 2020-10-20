// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2008.
// This file pertains to (S) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





template<typename SC>
MatrixElements_tmpl<SC> SymmetrySPU1_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1 + 1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletp.dat" << ", Iop=" << Invar(-1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-doubletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-doubletp.dat" << ", Iop=" << Invar(-1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-doubletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-1));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1-1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletm.dat" << ", Iop=" << Invar(+1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-doubletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+1));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-doubletm.dat" << ", Iop=" << Invar(+1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-doubletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+1));
      }(); // immediately executed lambda
    }
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
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletp.dat" << ", Iop=" << Invar(-1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-doubletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-1));
      }(); // immediately executed lambda
    }
  }
};

      Ip = Invar(ssz1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-doubletm.dat" << ", Iop=" << Invar(+1) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-doubletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+1));
      }(); // immediately executed lambda
    }
  }
};
    }
  }
  return cnew;
}

template<typename SC>
Opch_tmpl<SC> SymmetrySPU1_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  my_assert(!P.substeps);
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-spinupa.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-spinupa.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
 	           {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-spinupb.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };

     I1 = Invar(sszp-1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-spindowna.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-spindowna.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
             {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-spindownb.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return opch;
}

template<typename SC>
OpchChannel_tmpl<SC> SymmetrySPU1_tmpl<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, int M) {
  my_assert(P.substeps);
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    SZspin sszp = Ip.get("SSZ");
    Invar I1;

    I1 = Invar(sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spinupa.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-spinupa.dat"
        };
        opch[M][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spu1/spu1-1ch-spindowna.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-spindowna.dat"
        };
        opch[M][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
  }
  return opch[M];
}

template<typename SC>
MatrixElements_tmpl<SC> SymmetrySPU1_tmpl<SC>::recalc_triplet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      SZspin ssz1 = I1.get("SSZ");
      Invar Ip;

      Ip = Invar(ssz1);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-triplets.dat" << ", Iop=" << Invar(0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-triplets.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(0));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-triplets.dat" << ", Iop=" << Invar(0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-triplets.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(0));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1+2);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletp.dat" << ", Iop=" << Invar(-2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-tripletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-2));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-tripletp.dat" << ", Iop=" << Invar(-2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-tripletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-2));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

     Ip = Invar(ssz1-2);
     switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletm.dat" << ", Iop=" << Invar(+2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-tripletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+2));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-2ch-tripletm.dat" << ", Iop=" << Invar(+2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-2ch-tripletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+2));
      }(); // immediately executed lambda
    }
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
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-triplets.dat" << ", Iop=" << Invar(0) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-triplets.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(0));
      }(); // immediately executed lambda
    }
  }
};

      Ip = Invar(ssz1 + 2);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletp.dat" << ", Iop=" << Invar(-2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-tripletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(-2));
      }(); // immediately executed lambda
    }
  }
};

      Ip = Invar(ssz1 - 2);
      {
  nrglog('f', "RECALC(fn=" << "spu1/spu1-1ch-tripletm.dat" << ", Iop=" << Invar(+2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "spu1/spu1-1ch-tripletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(+2));
      }(); // immediately executed lambda
    }
  }
};
    }
  }
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

// 2-channel only

#undef QDIFF
#define QDIFF(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1UP
#define Q1UP(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2UP
#define Q2UP(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q1DO
#define Q1DO(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef Q2DO
#define Q2DO(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

template<typename SC>
void SymmetrySPU1_tmpl<SC>::recalc_global(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const std::string name, MatrixElements_tmpl<SC> &cnew) {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
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

// XXX: use m4 macros rather than C++ macros
#undef Q1
#undef Q2
#undef Q3
