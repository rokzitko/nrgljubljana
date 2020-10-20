// *** WARNING!!! Modify nrg-recalc-NONE.cc.m4, not nrg-recalc-NONE.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, April 2010
// This file pertains to the case with no symmetry.

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





template<typename SC>
Opch_tmpl<SC> SymmetryNONE_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1 = Invar();

    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "none/none-1ch-a-CR-DO.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-1ch-a-CR-DO.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "none/none-1ch-a-CR-UP.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-1ch-a-CR-UP.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
}; } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "none/none-2ch-a-CR-DO.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-2ch-a-CR-DO.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
	          {
  nrglog('f', "RECALC_F(fn=" << "none/none-2ch-b-CR-DO.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-2ch-b-CR-DO.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "none/none-2ch-a-CR-UP.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-2ch-a-CR-UP.dat"
        };
        opch[0][1][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
            {
  nrglog('f', "RECALC_F(fn=" << "none/none-2ch-b-CR-UP.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f_tmpl<SC>> recalc_table = {
#include "none/none-2ch-b-CR-UP.dat"
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

template<typename SC>
MatrixElements_tmpl<SC> SymmetryNONE_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Invar Ip = Invar();

    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "none/none-1ch-doublet.dat" << ", Iop=" << Invar() << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "none/none-1ch-doublet.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar());
      }(); // immediately executed lambda
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "none/none-2ch-doublet.dat" << ", Iop=" << Invar() << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_tmpl<SC>> recalc_table = {
#include "none/none-2ch-doublet.dat"
        };
        cnew[II] = this->recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, Invar());
      }(); // immediately executed lambda
    }
  }
} } break;
  default: my_assert_not_reached();
  };
  }
  return cnew;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

// Isospin operator need an appropriate phase factor (bipartite sublattice index) 
#define USEISOFACTOR

#if defined(USEISOFACTOR)
#define ISOFACTOR psgn(step.getnn() + 1)
#else
#define ISOFACTOR 1
#endif

#ifdef NRG_COMPLEX
#undef SPINY
#define SPINY(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINY
#define ISOSPINY(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *complex<double>(ISOFACTOR))

#undef Complex
#define Complex(x, y) cmpl(x, y)
#endif // NRG_COMPLEX

#undef CHARGE
#define CHARGE(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) this->recalc1_global(diag, qsrmax, I1, cn, i1, ip, value *ISOFACTOR)

template<typename SC>
void SymmetryNONE_tmpl<SC>::recalc_global(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const std::string name, MatrixElements_tmpl<SC> &cnew) {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-spinz.dat"
          break;
        case 2:
#include "none/none-2ch-spinz.dat"
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
        case 1:
#include "none/none-1ch-spiny.dat"
          break;
        case 2:
#include "none/none-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-spinx.dat"
          break;
        case 2:
#include "none/none-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Qtot.dat"
          break;
        case 2:
#include "none/none-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iztot.dat"
          break;
        case 2:
#include "none/none-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Ixtot.dat"
          break;
        case 2:
#include "none/none-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

#ifdef NRG_COMPLEX
  if (name == "Iytot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iytot.dat"
          break;
        case 2:
#include "none/none-2ch-Iytot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
#endif

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Iptot.dat"
          break;
        case 2:
#include "none/none-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "none/none-1ch-Imtot.dat"
          break;
        case 2:
#include "none/none-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  
  my_assert_not_reached();
}
