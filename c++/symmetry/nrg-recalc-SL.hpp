#include <vector>
#include <string>

namespace NRG {

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007
// This file pertains to the spinless-fermions code.

template<typename SC>
MatrixElements<SC> SymmetrySL<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold)const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1 = I1.get("Q");
    Invar Ip  = Invar(q1 - 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1), cnew, [&] { return RecalcTable<SC>{
#include "sl/sl-1ch-doublet.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1), cnew, [&] { return RecalcTable<SC>{
#include "sl/sl-2ch-doublet.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1), cnew, [&] { return RecalcTable<SC>{
#include "sl/sl-3ch-doublet.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySL<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    Invar I1  = Invar(qp + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-1ch-a.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-2ch-a.dat"
        }; }); 
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-2ch-b.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-3ch-a.dat"
        }; }); 
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-3ch-b.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "sl/sl-3ch-c.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return opch;
}

template<typename SC>
void SymmetrySL<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "sl/sl-2ch-qdiff.dat"
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
#include "sl/sl-2ch-qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "N3") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 3:
#include "sl/sl-3ch-N3.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
