#include <complex>
#include <vector>
#include <string>
#include <type_traits>

namespace NRG {

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2017
// This file pertains to the case with only fermion number parities
// (one per channel).

template<typename SC>
Opch<SC> SymmetryPP<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int pa   = Ip.get("Pa");
    int pb   = Ip.get("Pb");

    {
      Invar I1 = Invar(-pa, pb);
      this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-a-CR-DO.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-a-CR-UP.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[0][2], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-a-AN-DO.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[0][3], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-a-AN-UP.dat"
      }; });
    }

    {
      Invar I1 = Invar(pa, -pb);
      this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-b-CR-DO.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[1][1], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-b-CR-UP.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[1][2], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-b-AN-DO.dat"
      }; });
      this->recalc_f_table(diag, I1, Ip, opch[1][3], [&] { return RecalcFTable<SC>{
#include "pp/pp-2ch-b-AN-UP.dat"
      }; });
    }
  }
  return opch;
}

// The transverse isospin components need a phase factor (bipartite sublattice index), see GlobalTerms.

#undef Complex
#define Complex(x, y) cmpl(x, y)

template<typename SC>
void SymmetryPP<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
