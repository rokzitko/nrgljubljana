namespace NRG {

// *** WARNING!!! Modify nrg-recalc-PP.cc.m4, not nrg-recalc-PP.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2017
// This file pertains to the case with only fermion number parities
// (one per channel).

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
Opch<SC> SymmetryPP<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int pa   = Ip.get("Pa");
    int pb   = Ip.get("Pb");

    {
      Invar I1 = Invar(-pa, pb);
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-a-CR-DO.dat" << ", ch=" << 0 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-a-CR-DO.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-a-CR-UP.dat" << ", ch=" << 0 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-a-CR-UP.dat"
      };
      opch[0][1][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-a-AN-DO.dat" << ", ch=" << 0 << ", n=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-a-AN-DO.dat"
      };
      opch[0][2][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-a-AN-UP.dat" << ", ch=" << 0 << ", n=" << 3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-a-AN-UP.dat"
      };
      opch[0][3][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    }

    {
      Invar I1 = Invar(pa, -pb);
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-b-CR-DO.dat" << ", ch=" << 1 << ", n=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-b-CR-DO.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-b-CR-UP.dat" << ", ch=" << 1 << ", n=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-b-CR-UP.dat"
      };
      opch[1][1][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-b-AN-DO.dat" << ", ch=" << 1 << ", n=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-b-AN-DO.dat"
      };
      opch[1][2][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
      {
  nrglog('f', "RECALC_F(fn=" << "pp/pp-2ch-b-AN-UP.dat" << ", ch=" << 1 << ", n=" << 3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "pp/pp-2ch-b-AN-UP.dat"
      };
      opch[1][3][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
    }
  }
  return opch;
}

#undef SPINX
#define SPINX(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

// Isospin operator need an appropriate phase factor (bipartite sublattice index) 
#define USEISOFACTOR

#if defined(USEISOFACTOR)
#define ISOFACTOR psgn(step.getnn() + 1)
#else
#define ISOFACTOR 1
#endif

#undef SPINY
#define SPINY(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINY
#define ISOSPINY(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value * std::complex<double>(ISOFACTOR))

#undef Complex
#define Complex(x, y) cmpl(x, y)

#undef CHARGE
#define CHARGE(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *ISOFACTOR)

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *ISOFACTOR)

template<typename SC>
void SymmetryPP<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "pp/pp-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "SYtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "pp/pp-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "pp/pp-2ch-spinx.dat"
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
#include "pp/pp-2ch-Qtot.dat"
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
        case 2:
#include "pp/pp-2ch-Iztot.dat"
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
        case 2:
#include "pp/pp-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "Iytot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "pp/pp-2ch-Iytot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "pp/pp-2ch-Iptot.dat"
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
        case 2:
#include "pp/pp-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  
  my_assert_not_reached();
}

}
