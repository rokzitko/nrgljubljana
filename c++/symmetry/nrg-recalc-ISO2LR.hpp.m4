namespace NRG {

// *** WARNING!!! Modify nrg-recalc-ISOLR.cc.m4, not nrg-recalc-ISOLR.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, June 2006, Nov 2007, May 2008
// This file pertains to (I,S,P) subspaces
// Version for EVEN number of impurities

include(recalc-macros.m4)

double sign(double x) {
  if (x > 0.0) return +1.0;
  if (x < 0.0) return -1.0;
  my_assert_not_reached();
}

// (ISOLR): 8 calls of recalc_f() are necessary: different parities are also possible!

template<typename SC>
Opch<SC> SymmetryISO2LR<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels.

    int iip = Ip.get("II");
    int ssp = Ip.get("SS");

    // IMPORTANT NEW ELEMENT: the parity is important here!!
    int pp = Ip.get("P");

    // nn is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // Both parities yield non-zero <I+-1/2, S+-1/2, P| a^\mu_\nu
    // |I,S,P'>.  Coefficients *DO* depend on P,P', or more accurately,
    // on whether or not P and P' are the same.
    //
    // Observation: due to reflection symmetry, the coefficient for 'a' and
    // 'b' (2 channels) are either all the same or differ in sign.

    // ****** CASE I: SAME PARITY ******

    I1 = Invar(iip + 1, ssp + 1, pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isoupa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isoupb.dat", 1);

    I1 = Invar(iip + 1, ssp - 1, pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isoupa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isoupb.dat", 1);

    I1 = Invar(iip - 1, ssp + 1, pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isodowna.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isodownb.dat", 1);

    I1 = Invar(iip - 1, ssp - 1, pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isodowna.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isodownb.dat", 1);

    // ****** CASE II: DIFFERENT PARITY ******

    I1 = Invar(iip + 1, ssp + 1, -pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isoupdiffa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isoupdiffb.dat", 1);

    I1 = Invar(iip + 1, ssp - 1, -pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isoupdiffa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isoupdiffb.dat", 1);

    I1 = Invar(iip - 1, ssp + 1, -pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isodowndiffa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spinup-isodowndiffb.dat", 1);

    I1 = Invar(iip - 1, ssp - 1, -pp);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isodowndiffa.dat", 0);
    RECALC_F_TAB("iso2lr/iso2lr-2ch-spindown-isodowndiffb.dat", 1);
  }
  return opch;
}

// Recalculate matrix elements of a doublet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISO2LR<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1 - 1, ss1 + 1, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-doubletmp.dat", Invar(2, 2, +1));

    Ip = Invar(ii1 - 1, ss1 - 1, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-doubletmm.dat", Invar(2, 2, +1));

    Ip = Invar(ii1 + 1, ss1 + 1, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-doubletpp.dat", Invar(2, 2, +1));

    Ip = Invar(ii1 + 1, ss1 - 1, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-doubletpm.dat", Invar(2, 2, +1));
  }
  return cnew;
}

// Recalculate matrix elements of a triplet tensor operator [EVEN PARITY]
template<typename SC>
MatrixElements<SC> SymmetryISO2LR<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ii1 = I1.get("II");
    int ss1 = I1.get("SS");
    int p1    = I1.get("P");
    Invar Ip;

    Ip = Invar(ii1, ss1, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-triplets.dat", Invar(1, 3, +1));

    Ip = Invar(ii1, ss1 + 2, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-tripletp.dat", Invar(1, 3, +1));

    Ip = Invar(ii1, ss1 - 2, p1);
    RECALC_TAB("iso2lr/iso2lr-2ch-tripletm.dat", Invar(1, 3, +1));
  }
  return cnew;
}

}
