namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SPSU2C3.cc.m4, not nrg-recalc-SPSU2C3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (S,P) subspaces

include(recalc-macros.m4)

#define xRECALC_F_TAB(a, b, c) 0;

template<typename SC>
Opch<SC> SymmetrySPSU2C3<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax) {
  Opch<SC> opch = newopch<SC>(P);

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
 
  // CONVENTION: primed indeces are on the right side (ket)
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    int p     = Ip.get("P");

    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(ssp + 1, (p + 0) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup0-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup0-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup0-c.dat", 2);

    I1 = Invar(ssp - 1, (p + 0) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown0-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown0-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown0-c.dat", 2);

    I1 = Invar(ssp + 1, (p + 1) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup1-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup1-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup1-c.dat", 2);

    I1 = Invar(ssp - 1, (p + 1) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown1-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown1-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown1-c.dat", 2);

    I1 = Invar(ssp + 1, (p + 2) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup2-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup2-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spinup2-c.dat", 2);

    I1 = Invar(ssp - 1, (p + 2) % 3);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown2-a.dat", 0);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown2-b.dat", 1);
    RECALC_F_TAB("spsu2c3/spsu2c3-spindown2-c.dat", 2);
#undef sqrt
#undef Power
  }
  
  }
  return opch;
}

}
