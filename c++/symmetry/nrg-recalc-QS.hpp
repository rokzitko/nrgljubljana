#include <cmath>
#include <vector>
#include <string>

namespace NRG {

using std::sqrt;

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (Q,S) subspaces

template<typename SC>
MatrixElements<SC> SymmetryQS<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-doubletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-2ch-doubletp.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-3ch-doubletp.dat"
        }; });
        break;
      case 4:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-4ch-doubletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(q1-1, ss1-1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-doubletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-2ch-doubletm.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-3ch-doubletm.dat"
        }; });
        break;
      case 4:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-4ch-doubletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1 - 1, ss1 + 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-doubletp.dat"
      }; });

      Ip = Invar(q1 - 1, ss1 - 1);
      this->recalc_table(diag, substruct, cold, I1, Ip, Invar(1, 2), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-doubletm.dat"
      }; });
    }
  }
  return cnew;
}

// (QS): Two calls of recalc_f() are necessary (for S+1/2 and S-1/2)
// for each channel.
// See Krishna-Murthy p. 1034, equation (B10).

// <I1 | f^dag | Ip>, hence q1=qp+1

// Driver routine for recalc_f()
template<typename SC>
Opch<SC> SymmetryQS<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    // NOTE: q,ss only couples to q+1,ss+-1 in general, even for
    // several channels.

    I1 = Invar(qp + 1, ssp + 1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-1ch-spinupa.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-2ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-2ch-spinupb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spinupb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spinupc.dat"
        }; });
        break;
      case 4:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spinupa.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spinupb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spinupc.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[3][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spinupd.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    I1 = Invar(qp+1, ssp-1);
    switch (P.channels) {
      case 1:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-1ch-spindowna.dat"
        }; });
        break;
      case 2:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-2ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-2ch-spindownb.dat"
        }; });
        break;
      case 3:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spindownb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-3ch-spindownc.dat"
        }; });
        break;
      case 4:
        this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spindowna.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spindownb.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spindownc.dat"
        }; });
        this->recalc_f_table(diag, I1, Ip, opch[3][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-4ch-spindownd.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
  }
  return opch;
}
 
// Driver routine for recalc_f() for substeps=true, i.e., chain by chain diagonalisations
template<typename SC>
OpchChannel<SC> SymmetryQS<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const int M) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int qp = Ip.get("Q");
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(qp + 1, ssp + 1);
    this->recalc_f_table(diag, I1, Ip, opch[M][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-1ch-spinupa.dat"
    }; });

    I1 = Invar(qp + 1, ssp - 1);
    this->recalc_f_table(diag, I1, Ip, opch[M][0], [&] { return RecalcFTable<SC>{
#include "qs/qs-1ch-spindowna.dat"
    }; });
  }
  return opch[M];
}

// Recalculate matrix elements of a triplet tenzor operator
template<typename SC>
MatrixElements<SC> SymmetryQS<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int q1 = I1.get("Q");
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(q1, ss1);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-triplets.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-2ch-triplets.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-3ch-triplets.dat"
        }; });
        break;
      case 4:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-4ch-triplets.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(q1, ss1+2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-tripletp.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-2ch-tripletp.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-3ch-tripletp.dat"
        }; });
        break;
      case 4:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-4ch-tripletp.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }

    Ip = Invar(q1, ss1-2);
    switch (P.channels) {
      case 1:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-1ch-tripletm.dat"
        }; });
        break;
      case 2:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-2ch-tripletm.dat"
        }; });
        break;
      case 3:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-3ch-tripletm.dat"
        }; });
        break;
      case 4:
        this->recalc_table(diag, substruct, cold, I1, Ip, Invar(0, 3), cnew, [&] { return RecalcTable<SC>{
#include "qs/qs-4ch-tripletm.dat"
        }; });
        break;
      default: my_assert_not_reached();
    }
    }
  } else my_assert_not_reached();
  return cnew;
}

template<typename SC>
void SymmetryQS<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, 
                                        std::string name, MatrixElements<SC> &cnew) const {
  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      const GlobalTerms<SC> terms(*this, diag, I1, cn, psgn(step.getnn() + 1));
      switch (P.channels) {
        case 2:
#include "qs/qs-2ch-qdiff.dat"
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
#include "qs/qs-2ch-q1.dat"
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
#include "qs/qs-2ch-q2.dat"
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
#include "qs/qs-2ch-qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
