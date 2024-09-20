namespace NRG {

// *** WARNING!!! Modify nrg-recalc-P.cc.m4, not nrg-recalc-P.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, July 2017
// This file pertains to the case with only fermion number parity.

include(recalc-macros.m4)

template<typename SC>
Opch<SC> SymmetryP<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int p    = Ip.get("P");

    Invar I1 = Invar(-p); // always the opposite fermion parity!

    ONETWO(`RECALC_F_TAB_N("p/p-1ch-a-CR-DO.dat", 0, 0);
            RECALC_F_TAB_N("p/p-1ch-a-CR-UP.dat", 0, 1);
            RECALC_F_TAB_N("p/p-1ch-a-AN-DO.dat", 0, 2);
            RECALC_F_TAB_N("p/p-1ch-a-AN-UP.dat", 0, 3);',

           `RECALC_F_TAB_N("p/p-2ch-a-CR-DO.dat", 0, 0);
	          RECALC_F_TAB_N("p/p-2ch-b-CR-DO.dat", 1, 0);
            RECALC_F_TAB_N("p/p-2ch-a-CR-UP.dat", 0, 1);
            RECALC_F_TAB_N("p/p-2ch-b-CR-UP.dat", 1, 1);
    	      RECALC_F_TAB_N("p/p-2ch-a-AN-DO.dat", 0, 2);
	          RECALC_F_TAB_N("p/p-2ch-b-AN-DO.dat", 1, 2);
            RECALC_F_TAB_N("p/p-2ch-a-AN-UP.dat", 0, 3);
            RECALC_F_TAB_N("p/p-2ch-b-AN-UP.dat", 1, 3);');
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryP<SC>::recalc_doublet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int p1   = I1.get("P");
    Invar Ip = Invar(-p1); // always the opposite fermion parity!

    ONETWO(`RECALC_TAB("p/p-1ch-doublet.dat", Invar(-1))',
           `RECALC_TAB("p/p-2ch-doublet.dat", Invar(-1))');
  }
  return cnew;
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
void SymmetryP<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  if (name == "SZtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-spinz.dat"
          break;
        case 2:
#include "p/p-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "SYtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-spiny.dat"
          break;
        case 2:
#include "p/p-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "SXtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-spinx.dat"
          break;
        case 2:
#include "p/p-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-Qtot.dat"
          break;
        case 2:
#include "p/p-2ch-Qtot.dat"
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
#include "p/p-1ch-Iztot.dat"
          break;
        case 2:
#include "p/p-2ch-Iztot.dat"
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
#include "p/p-1ch-Ixtot.dat"
          break;
        case 2:
#include "p/p-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
  if (name == "Iytot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-Iytot.dat"
          break;
        case 2:
#include "p/p-2ch-Iytot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II {I1, I1};
      Matrix &cn = cnew[II];
      switch (P.channels) {
        case 1:
#include "p/p-1ch-Iptot.dat"
          break;
        case 2:
#include "p/p-2ch-Iptot.dat"
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
#include "p/p-1ch-Imtot.dat"
          break;
        case 2:
#include "p/p-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }
  
  my_assert_not_reached();
}

}
