namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SL3.cc.m4, not nrg-recalc-SL3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007, Oct 2010
// This file pertains to the spinless-fermions code.

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetrySL3<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q11 = I1.get("Q1");
    int q21 = I1.get("Q2");
    int q31 = I1.get("Q3");
    Invar Ip   = Invar(q11 - 1, q21, q31); // This is a channel 1 operator
    RECALC_TAB("sl3/sl3-3ch-doublet.dat", Invar(1, 0, 0));
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySL3<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    int q1p = Ip.get("Q1");
    int q2p = Ip.get("Q2");
    int q3p = Ip.get("Q3");

    Invar I1;

    I1 = Invar(q1p + 1, q2p, q3p);
    RECALC_F_TAB("sl3/sl3-3ch-a.dat", 0);

    I1 = Invar(q1p, q2p + 1, q3p);
    RECALC_F_TAB("sl3/sl3-3ch-b.dat", 1);

    I1 = Invar(q1p, q2p, q3p + 1);
    RECALC_F_TAB("sl3/sl3-3ch-c.dat", 2);
  }
  return opch;
}

#undef QTOT
#define QTOT(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

#undef N1
#define N1(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

#undef N2
#define N2(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

#undef N3
#define N3(i1, ip, ch, value) this->recalc1_global(diag, substruct, I1, cn, i1, ip, value)

template<typename SC>
void SymmetrySL3<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const std::string name, MatrixElements<SC> &cnew) {
  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
#include "sl3/sl3-3ch-qtot.dat"
    }
    return;
  }

  if (name == "N1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
#include "sl3/sl3-3ch-N1.dat"
    }
    return;
  }

  if (name == "N2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
#include "sl3/sl3-3ch-N2.dat"
    }
    return;
  }

  if (name == "N3") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
#include "sl3/sl3-3ch-N3.dat"
    }
    return;
  }

  my_assert_not_reached();
}

}
