namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SPSU2T.cc.m4, not nrg-recalc-SPSU2T.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Aug 2015
// This file pertains to (S,T) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetrySPSU2T<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double T   = t1; // trick!
    Invar Ip;

    // Two different lengths: D_3CH_a and D_3CH_b

    Ip = Invar(ss1 + 1, t1 - 1);
    RECALC_TAB("spsu2t/spsu2t-doubletp-1.dat", Invar(1, 2, 1));

    Ip = Invar(ss1 - 1, t1 - 1);
    RECALC_TAB("spsu2t/spsu2t-doubletm-1.dat", Invar(1, 2, 1));

    Ip = Invar(ss1 + 1, t1);
    RECALC_TAB("spsu2t/spsu2t-doubletp0.dat", Invar(1, 2, 1));

    Ip = Invar(ss1 - 1, t1);
    RECALC_TAB("spsu2t/spsu2t-doubletm0.dat", Invar(1, 2, 1));

    Ip = Invar(ss1 + 1, t1 + 1);
    RECALC_TAB("spsu2t/spsu2t-doubletp+1.dat", Invar(1, 2, 1));

    Ip = Invar(ss1 - 1, t1 + 1);
    RECALC_TAB("spsu2t/spsu2t-doubletm+1.dat", Invar(1, 2, 1));
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

template<typename SC>
Opch<SC> SymmetrySPSU2T<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) const {
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp  = Ip.get("SS");
    int tp = Ip.get("T");
    double T   = tp; // trick!
    Invar I1;

    // The different files just correspond to contributions computed
    // for various d[CR,sz,tz] operators.
    // Check: there should not be any lines with equal subspaces
    // indexes in different files!! That's indeed the case for the
    // generated files for symtype=SPSU2T.
    I1 = Invar(ssp + 1, tp + 1);
    RECALC_F_TAB("spsu2t/spsu2t-spinup+1.dat", 0);

    I1 = Invar(ssp + 1, tp);
    RECALC_F_TAB("spsu2t/spsu2t-spinup0.dat", 0);

    I1 = Invar(ssp + 1, tp - 1);
    RECALC_F_TAB("spsu2t/spsu2t-spinup-1.dat", 0);

    I1 = Invar(ssp - 1, tp + 1);
    RECALC_F_TAB("spsu2t/spsu2t-spindo+1.dat", 0);

    I1 = Invar(ssp - 1, tp);
    RECALC_F_TAB("spsu2t/spsu2t-spindo0.dat", 0);

    I1 = Invar(ssp - 1, tp - 1);
    RECALC_F_TAB("spsu2t/spsu2t-spindo-1.dat", 0);
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetrySPSU2T<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int ss1  = I1.get("SS");
    int t1 = I1.get("T");
    double T   = t1; // trick!
    Invar Ip;

    Ip = Invar(ss1, t1);
    RECALC_TAB("spsu2t/spsu2t-triplets.dat", Invar(3, 0));

    Ip = Invar(ss1 + 2, t1);
    RECALC_TAB("spsu2t/spsu2t-tripletp.dat", Invar(3, 0));

    Ip = Invar(ss1 - 2, t1);
    RECALC_TAB("spsu2t/spsu2t-tripletm.dat", Invar(3, 0));
  }
  return cnew;
}

}
