// *** WARNING!!! Modify nrg-recalc-ISO.cc.m4, not nrg-recalc-ISO.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006
// This file pertains to (I,S) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements_tmpl<SC> SymmetryISO_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1 - 1, ss1 + 1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletmp.dat", Invar(2, 2))',
           `RECALC_TAB("iso/iso-2ch-doubletmp.dat", Invar(2, 2))');

    Ip = Invar(ii1-1, ss1-1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletmm.dat", Invar(2, 2))',
    	   `RECALC_TAB("iso/iso-2ch-doubletmm.dat", Invar(2, 2))');

    Ip = Invar(ii1+1, ss1+1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletpp.dat", Invar(2, 2))',
    	   `RECALC_TAB("iso/iso-2ch-doubletpp.dat", Invar(2, 2))');

    Ip = Invar(ii1+1, ss1-1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-doubletpm.dat", Invar(2, 2))',
   	   `RECALC_TAB("iso/iso-2ch-doubletpm.dat", Invar(2, 2))');
  }
  return cnew;
}

template<typename SC>
Opch_tmpl<SC> SymmetryISO_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
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
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spinup-isoupa.dat", 0)',

          `RECALC_F_TAB("iso/iso-2ch-spinup-isoupa.dat", 0);
      	   RECALC_F_TAB("iso/iso-2ch-spinup-isoupb.dat", 1)',

          `RECALC_F_TAB("iso/iso-3ch-spinup-isoupa.dat", 0);
	         RECALC_F_TAB("iso/iso-3ch-spinup-isoupb.dat", 1);
	         RECALC_F_TAB("iso/iso-3ch-spinup-isoupc.dat", 2);');
    
    I1 = Invar(iip+1, ssp-1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spindown-isoupa.dat", 0)',

          `RECALC_F_TAB("iso/iso-2ch-spindown-isoupa.dat", 0);
	         RECALC_F_TAB("iso/iso-2ch-spindown-isoupb.dat", 1)',
	
          `RECALC_F_TAB("iso/iso-3ch-spindown-isoupa.dat", 0);
	         RECALC_F_TAB("iso/iso-3ch-spindown-isoupb.dat", 1);
	         RECALC_F_TAB("iso/iso-3ch-spindown-isoupc.dat", 2);');

    I1 = Invar(iip-1, ssp+1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spinup-isodowna.dat", 0)',

          `RECALC_F_TAB("iso/iso-2ch-spinup-isodowna.dat", 0);
	         RECALC_F_TAB("iso/iso-2ch-spinup-isodownb.dat", 1)',

          `RECALC_F_TAB("iso/iso-3ch-spinup-isodowna.dat", 0);
	         RECALC_F_TAB("iso/iso-3ch-spinup-isodownb.dat", 1);
	         RECALC_F_TAB("iso/iso-3ch-spinup-isodownc.dat", 2);');

    I1 = Invar(iip-1, ssp-1);
    ONE23(`RECALC_F_TAB("iso/iso-1ch-spindown-isodowna.dat", 0)',

          `RECALC_F_TAB("iso/iso-2ch-spindown-isodowna.dat", 0);
      	   RECALC_F_TAB("iso/iso-2ch-spindown-isodownb.dat", 1)',

          `RECALC_F_TAB("iso/iso-3ch-spindown-isodowna.dat", 0);
	         RECALC_F_TAB("iso/iso-3ch-spindown-isodownb.dat", 1);
	         RECALC_F_TAB("iso/iso-3ch-spindown-isodownc.dat", 2);');
  }
  return opch;
}

template<typename SC>
MatrixElements_tmpl<SC> SymmetryISO_tmpl<SC>::recalc_triplet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii1 = I1.get("II");
    Sspin ss1 = I1.get("SS");
    Invar Ip;

    Ip = Invar(ii1, ss1);
    ONETWO(`RECALC_TAB("iso/iso-1ch-triplets.dat", Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-triplets.dat", Invar(1, 3))');

    Ip = Invar(ii1, ss1+2);
    ONETWO(`RECALC_TAB("iso/iso-1ch-tripletp.dat", Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-tripletp.dat", Invar(1, 3))');

    Ip = Invar(ii1, ss1-2);
    ONETWO(`RECALC_TAB("iso/iso-1ch-tripletm.dat", Invar(1, 3))',
           `RECALC_TAB("iso/iso-2ch-tripletm.dat", Invar(1, 3))');
  }
  return cnew;
}
