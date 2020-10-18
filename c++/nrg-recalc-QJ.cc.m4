// *** WARNING!!! Modify nrg-recalc-QJ.cc.m4, not nrg-recalc-QJ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2016
// This file pertains to (Q,J) subspaces

include(recalc-macros.m4)

// Recalculate matrix elements of a doublet tensor operator
template<typename SC>
MatrixElements_tmpl<SC> SymmetryQJ_tmpl<SC>::recalc_doublet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 1);
    RECALC_TAB("qj/qj-doubletp.dat", Invar(1, 2));

    Ip = Invar(q1 - 1, jj1 - 1);
    RECALC_TAB("qj/qj-doubletm.dat", Invar(1, 2));
  }
  return cnew;
}

#undef If
#define If(cond, a, b) (cond ? a : b)

// Recalculate matrix elements of a quadruplet tensor operator
template<typename SC>
MatrixElements_tmpl<SC> SymmetryQJ_tmpl<SC>::recalc_quadruplet(const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax, const MatrixElements_tmpl<SC> &cold) {
  MatrixElements_tmpl<SC> cnew;
  for(const auto &[I1, eig]: diag) {
    Number q1 = I1.get("Q");
    Sspin jj1 = I1.get("JJ");
    //    double J = (jj1-1.0)/2.0;
    Invar Ip;

    Ip = Invar(q1 - 1, jj1 + 3);
    RECALC_TAB("qj/qj-quad1.dat", Invar(1, 4));

    Ip = Invar(q1 - 1, jj1 + 1);
    RECALC_TAB("qj/qj-quad2.dat", Invar(1, 4));

    Ip = Invar(q1 - 1, jj1 - 1);
    RECALC_TAB("qj/qj-quad3.dat", Invar(1, 4));

    Ip = Invar(q1 - 1, jj1 - 3);
    RECALC_TAB("qj/qj-quad4.dat", Invar(1, 4));
  }
  return cnew;
}

template<typename SC>
Opch_tmpl<SC> SymmetryQJ_tmpl<SC>::recalc_irreduc(const Step &step, const DiagInfo_tmpl<SC> &diag, const QSrmax &qsrmax) {
  Opch_tmpl<SC> opch = newopch<SC>(P);
  for(const auto &[Ip, eig]: diag) {
    Number qp = Ip.get("Q");
    Sspin jjp = Ip.get("JJ");
    double j  = J(jjp);
    Invar I1;

    I1 = Invar(qp + 1, jjp + 3);
    RECALC_F_TAB("qj/qj-spin_j3_2-jz3_2.dat", 1);

    I1 = Invar(qp + 1, jjp + 1);
    RECALC_F_TAB("qj/qj-spin_j1_2-jz1_2.dat", 0);
    RECALC_F_TAB("qj/qj-spin_j3_2-jz1_2.dat", 1);

    I1 = Invar(qp + 1, jjp - 1);
    RECALC_F_TAB("qj/qj-spin_j1_2-jz-1_2.dat", 0);
    RECALC_F_TAB("qj/qj-spin_j3_2-jz-1_2.dat", 1);

    I1 = Invar(qp + 1, jjp - 3);
    RECALC_F_TAB("qj/qj-spin_j3_2-jz-3_2.dat", 1);
  }
  return opch;
}
