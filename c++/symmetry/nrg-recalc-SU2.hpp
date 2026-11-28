#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <type_traits>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Sep 2009
// This file pertains to (I) subspaces

template<typename SC>
MatrixElements<SC> SymmetrySU2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    Invar Ip;

    Ip = Invar(ii1 - 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "su2/su2-1ch-doubletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "su2/su2-2ch-doubletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ii1+1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "su2/su2-1ch-doubletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "su2/su2-2ch-doubletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySU2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int iip = Ip.get("II");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(iip + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-1ch-type1-isoup-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-1ch-type2-isoup-a.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type1-isoup-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type1-isoup-b.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type2-isoup-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type2-isoup-b.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(iip-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-1ch-type1-isodown-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-1ch-type2-isodown-a.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type1-isodown-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type1-isodown-b.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type2-isodown-a.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][1], [&] { return RecalcFTable<SC>{
#include "su2/su2-2ch-type2-isodown-b.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return opch;
}

template<typename SC>
void SymmetrySU2<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
#undef Complex
#define Complex(x, y) cmpl(x, y)
  if (name == "SYtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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
  }

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
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

}
