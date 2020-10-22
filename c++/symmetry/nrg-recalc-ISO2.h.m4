// *** WARNING!!! Modify nrg-recalc-ISO2.cc.m4, not nrg-recalc-ISO2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, May 2008
// This file pertains to (I,S) subspaces
// Version for EVEN number of impurities.

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetryISO2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1 - 1, ss1 + 1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletmp.dat", Invar(2, 2))',
           `RECALC_TAB("iso2/iso2-2ch-doubletmp.dat", Invar(2, 2))');

    Ip = Invar(ii1-1, ss1-1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletmm.dat", Invar(2, 2))',
    	   `RECALC_TAB("iso2/iso2-2ch-doubletmm.dat", Invar(2, 2))');

    Ip = Invar(ii1+1, ss1+1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletpp.dat", Invar(2, 2))',
    	   `RECALC_TAB("iso2/iso2-2ch-doubletpp.dat", Invar(2, 2))');

    Ip = Invar(ii1+1, ss1-1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-doubletpm.dat", Invar(2, 2))',
   	   `RECALC_TAB("iso2/iso2-2ch-doubletpm.dat", Invar(2, 2))');
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetryISO2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax) {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    // NOTE: ii,ss only couples to ii+-1,ss+-1 in general, even for
    // several channels.

    Ispin iip = Ip.get("II");
    Sspin ssp = Ip.get("SS");
    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    I1 = Invar(iip + 1, ssp + 1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spinup-isoupa.dat", 0)',

           `RECALC_F_TAB("iso2/iso2-2ch-spinup-isoupa.dat", 0);
     	      RECALC_F_TAB("iso2/iso2-2ch-spinup-isoupb.dat", 1)');

    I1 = Invar(iip+1, ssp-1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spindown-isoupa.dat", 0)',

           `RECALC_F_TAB("iso2/iso2-2ch-spindown-isoupa.dat", 0);
	          RECALC_F_TAB("iso2/iso2-2ch-spindown-isoupb.dat", 1)');

    I1 = Invar(iip-1, ssp+1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spinup-isodowna.dat", 0)',

           `RECALC_F_TAB("iso2/iso2-2ch-spinup-isodowna.dat", 0);
	          RECALC_F_TAB("iso2/iso2-2ch-spinup-isodownb.dat", 1)');

    I1 = Invar(iip-1, ssp-1);
    ONETWO(`RECALC_F_TAB("iso2/iso2-1ch-spindown-isodowna.dat", 0)',

    	     `RECALC_F_TAB("iso2/iso2-2ch-spindown-isodowna.dat", 0);
	          RECALC_F_TAB("iso2/iso2-2ch-spindown-isodownb.dat", 1)');
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryISO2<SC>::recalc_triplet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1, ss1);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-triplets.dat", Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-triplets.dat", Invar(1, 3))');

    Ip = Invar(ii1, ss1+2);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-tripletp.dat", Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-tripletp.dat", Invar(1, 3))');

    Ip = Invar(ii1, ss1-2);
    ONETWO(`RECALC_TAB("iso2/iso2-1ch-tripletm.dat", Invar(1, 3))',
           `RECALC_TAB("iso2/iso2-2ch-tripletm.dat", Invar(1, 3))');
  }
  return cnew;
}
