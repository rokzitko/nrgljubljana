#include <vector>
#include <string>

namespace NRG {

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007, Oct 2010
// This file pertains to the spinless-fermions code.

template<typename SC>
MatrixElements<SC> SymmetrySL3<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q11 = I1.get("Q1");
    int q21 = I1.get("Q2");
    int q31 = I1.get("Q3");
    Invar Ip   = Invar(q11 - 1, q21, q31); // This is a channel 1 operator
    this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 0, 0), cnew, [&] { return RecalcTable<SC>{
#include "sl3/sl3-3ch-doublet.dat"
    }; });
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySL3<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int q1p = Ip.get("Q1");
    int q2p = Ip.get("Q2");
    int q3p = Ip.get("Q3");

    Invar I1;

    I1 = Invar(q1p + 1, q2p, q3p);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "sl3/sl3-3ch-a.dat"
    }; });

    I1 = Invar(q1p, q2p + 1, q3p);
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "sl3/sl3-3ch-b.dat"
    }; });

    I1 = Invar(q1p, q2p, q3p + 1);
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "sl3/sl3-3ch-c.dat"
    }; });
  }
  return opch;
}

template<typename SC>
void SymmetrySL3<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
#include "sl3/sl3-3ch-qtot.dat"
    }
    return;
  }

  if (name == "N1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
#include "sl3/sl3-3ch-N1.dat"
    }
    return;
  }

  if (name == "N2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
#include "sl3/sl3-3ch-N2.dat"
    }
    return;
  }

  if (name == "N3") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
#include "sl3/sl3-3ch-N3.dat"
    }
    return;
  }

  my_assert_not_reached();
}

}
