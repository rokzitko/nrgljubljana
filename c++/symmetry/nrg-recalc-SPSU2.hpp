#include <cmath>
#include <vector>
#include <string>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Dec 2007
// This file pertains to (S) subspaces

template<typename SC>
MatrixElements<SC> SymmetrySPSU2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-doubletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-2ch-doubletp.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-3ch-doubletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ss1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-doubletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-2ch-doubletm.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-3ch-doubletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-doubletp.dat"
      }; });

      Ip = Invar(ss1 - 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(2), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-doubletm.dat"
      }; });
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySPSU2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  my_assert(!P.substeps);
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-1ch-spinupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-2ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-2ch-spinupb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spinupb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spinupc.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(ssp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-1ch-spindowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-2ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-2ch-spindownb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spindownb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-3ch-spindownc.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    // Note: for 3ch cases, the lengths in the three channels are not the same!
    // The same thing occurs for all SU(2)_spin cases, for instance for symtype=QS.
    // RZ, oct 2015
  }
  return opch;
}

template<typename SC>
OpchChannel<SC> SymmetrySPSU2<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const int M) const {
  my_assert(P.substeps);
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[M][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-1ch-spinupa.dat"
    }; });

    I1 = Invar(ssp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[M][0], [&] { return RecalcFTable<SC>{
#include "spsu2/spsu2-1ch-spindowna.dat"
    }; });
  }
  return opch[M];
}

template<typename SC>
MatrixElements<SC> SymmetrySPSU2<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-triplets.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-2ch-triplets.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-3ch-triplets.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ss1+2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-tripletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-2ch-tripletp.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-3ch-tripletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(ss1-2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-1ch-tripletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(3), cnew, [&] { return RecalcTable<SC>{
#include "spsu2/spsu2-2ch-tripletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
    }
  } else my_assert_not_reached();
  return cnew;
}

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

template<typename SC>
void SymmetrySPSU2<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  // NOTE: none of these are implemented for substeps==true.

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Qtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Qtot.dat"
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-qdiff.dat"
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q1.dat"
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q2.dat"
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
        case 1:
#include "spsu2/spsu2-1ch-Iztot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iztot.dat"
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
        case 1:
#include "spsu2/spsu2-1ch-Ixtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Ixtot.dat"
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
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iptot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iptot.dat"
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
        case 1:
#include "spsu2/spsu2-1ch-Imtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
