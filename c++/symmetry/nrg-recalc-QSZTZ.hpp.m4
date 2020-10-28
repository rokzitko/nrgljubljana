namespace NRG {

// *** WARNING!!! Modify nrg-recalc-QSZTZ.cc.m4, not nrg-recalc-QSZTZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,Sz,Tz) subspaces

include(recalc-macros.m4)

template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1   = I1.get("Q");
    int ssz1  = I1.get("SZ");
    int tz1 = I1.get("TZ");
    Invar Ip;

    // Invar(1,2,+-1,0) is correct. 1 = add charge, 2 = doublet,
    // 1 = triplet (because working with abs orbital momentum QNs)

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 - 1);
    RECALC_TAB("qsztz/qsztz-doubletp-1.dat", Invar(1, -1, +1));

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 - 1);
    RECALC_TAB("qsztz/qsztz-doubletm-1.dat", Invar(1, +1, +1));

    Ip = Invar(q1 - 1, ssz1 + 1, tz1);
    RECALC_TAB("qsztz/qsztz-doubletp0.dat", Invar(1, -1, 0));

    Ip = Invar(q1 - 1, ssz1 - 1, tz1);
    RECALC_TAB("qsztz/qsztz-doubletm0.dat", Invar(1, +1, 0));

    Ip = Invar(q1 - 1, ssz1 + 1, tz1 + 1);
    RECALC_TAB("qsztz/qsztz-doubletp+1.dat", Invar(1, -1, -1));

    Ip = Invar(q1 - 1, ssz1 - 1, tz1 + 1);
    RECALC_TAB("qsztz/qsztz-doubletm+1.dat", Invar(1, +1, -1));
  }
  return cnew;
}

// ch=1 <-> Tz=+1
// ch=2 <-> Tz=0
// ch=3 <-> Tz=-1

template<typename SC>
Opch<SC> SymmetryQSZTZ<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const SubspaceStructure &substruct) {
  Opch<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    int qp   = Ip.get("Q");
    int sszp  = Ip.get("SZ");
    int tzp = Ip.get("TZ");
    Invar I1;

    I1 = Invar(qp + 1, sszp + 1, tzp + 1);
    RECALC_F_TAB("qsztz/qsztz-spinup+1.dat", 0);

    I1 = Invar(qp + 1, sszp + 1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spinup0.dat", 0);

    I1 = Invar(qp + 1, sszp + 1, tzp - 1);
    RECALC_F_TAB("qsztz/qsztz-spinup-1.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp + 1);
    RECALC_F_TAB("qsztz/qsztz-spindo+1.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp);
    RECALC_F_TAB("qsztz/qsztz-spindo0.dat", 0);

    I1 = Invar(qp + 1, sszp - 1, tzp - 1);
    RECALC_F_TAB("qsztz/qsztz-spindo-1.dat", 0);
  }
  return opch;
}

template<typename SC>
MatrixElements<SC> SymmetryQSZTZ<SC>::recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) {
  MatrixElements<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    int q1   = I1.get("Q");
    int ssz1  = I1.get("SZ");
    int tz1 = I1.get("TZ");
    Invar Ip;

    Ip = Invar(q1, ssz1, tz1);
    RECALC_TAB("qsztz/qsztz-triplets.dat", Invar(0, 3, 0));

    Ip = Invar(q1, ssz1 + 2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletp.dat", Invar(0, 3, 0));

    Ip = Invar(q1, ssz1 - 2, tz1);
    RECALC_TAB("qsztz/qsztz-tripletm.dat", Invar(0, 3, 0));
  }
  return cnew;
}

}
