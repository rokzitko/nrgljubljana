#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <type_traits>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2009
// This file pertains to (I1,I2) subspaces

template<typename SC>
MatrixElements<SC> SymmetryDBLSU2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii11 = I1.get("II1");
    int ii21 = I1.get("II2");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 1), cnew, [&] { return RecalcTable<SC>{
#include "dblsu2/dblsu2-2ch-doubletm0.dat"
    }; });

    Ip = Invar(ii11 + 1, ii21);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2, 1), cnew, [&] { return RecalcTable<SC>{
#include "dblsu2/dblsu2-2ch-doubletp0.dat"
    }; });

    Ip = Invar(ii11, ii21 - 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "dblsu2/dblsu2-2ch-doublet0m.dat"
    }; });

    Ip = Invar(ii11, ii21 + 1);
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "dblsu2/dblsu2-2ch-doublet0p.dat"
    }; });
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryDBLSU2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    int ii1p = Ip.get("II1");
    int ii2p = Ip.get("II2");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(ii1p + 1, ii2p);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type1-isoup-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type2-isoup-a.dat"
    }; });

    I1 = Invar(ii1p, ii2p + 1);
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type1-isoup-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][1], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type2-isoup-b.dat"
    }; });

    I1 = Invar(ii1p - 1, ii2p);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type1-isodown-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[0][1], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type2-isodown-a.dat"
    }; });

    I1 = Invar(ii1p, ii2p - 1);
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type1-isodown-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][1], [&] { return RecalcFTable<SC>{
#include "dblsu2/dblsu2-2ch-type2-isodown-b.dat"
    }; });
  }
  return opch;
}

#undef Complex
#define Complex(x, y) cmpl(x, y)

template<typename SC>
void SymmetryDBLSU2<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "SYtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "SXtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
