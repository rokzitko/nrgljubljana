// *** WARNING!!! Modify nrg-recalc-SU2.cc.m4, not nrg-recalc-SU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Sep 2009
// This file pertains to (I) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





template<typename SC>
MatrixElements_tmpl<SC> SymmetrySU2_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1 = I1.get("II");
    Invar Ip;

    Ip = Invar(ii1 - 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "su2/su2-1ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-doubletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "su2/su2-2ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-doubletm.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ii1+1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "su2/su2-1ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-doubletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2));
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "su2/su2-2ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-doubletp.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar(2));
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

template<typename SC>
Opch_tmpl<SC> SymmetrySU2_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    Ispin iip = Ip.get("II");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(iip + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-1ch-type1-isoup-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-type1-isoup-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-1ch-type2-isoup-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-type2-isoup-a.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type1-isoup-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type1-isoup-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type1-isoup-b.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type1-isoup-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type2-isoup-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type2-isoup-a.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type2-isoup-b.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type2-isoup-b.dat"
        };
        opch[1][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(iip-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-1ch-type1-isodown-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-type1-isodown-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-1ch-type2-isodown-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-1ch-type2-isodown-a.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type1-isodown-a.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type1-isodown-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type1-isodown-b.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type1-isodown-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type2-isodown-a.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type2-isodown-a.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "su2/su2-2ch-type2-isodown-b.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "su2/su2-2ch-type2-isodown-b.dat"
        };
        opch[1][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return opch;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef SPINY
#define SPINY(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

template<typename SC>
void SymmetrySU2_tmpl<SC>::recalc_global(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const std::string name, MatrixElements_tmpl<SC> &cnew) {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "su2/su2-1ch-spinz.dat"
          break;
        case 2:
#include "su2/su2-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

#ifdef NRG_COMPLEX
#undef Complex
#define Complex(x, y) cmpl(x, y)
  if (name == "SYtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "su2/su2-1ch-spiny.dat"
          break;
        case 2:
#include "su2/su2-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif // NRG_COMPLEX

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "su2/su2-1ch-spinx.dat"
          break;
        case 2:
#include "su2/su2-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
