class SymmetrySPU1 : public SymField {
 private:
   outfield Sz2, Sz;
   
 public:
   template<typename ... Args> SymmetrySPU1(Args&& ... args) : SymField(std::forward<Args>(args)...),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2) {
       InvarStructure InvStruc[] = {
         {"SSZ", additive} // spin projection
       };
       initInvar(InvStruc, ARRAYLENGTH(InvStruc));
       InvarSinglet = Invar(0);
     }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1  = I1.get("SSZ");
    SZspin sszp  = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) override { return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")); }

  void load() override {
    if (!P.substeps) {
      switch (P.channels) {
        case 1:
#include "spu1/spu1-1ch-In2.dat"
#include "spu1/spu1-1ch-QN.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-In2.dat"
#include "spu1/spu1-2ch-QN.dat"
          break;
        default: my_assert_not_reached();
      } // switch
    } else {
#include "spu1/spu1-1ch-In2.dat"
#include "spu1/spu1-1ch-QN.dat"
    }
  }

  void make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef);
  void make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef);

  void calculate_TD(const Step &step, const DiagInfo &diag, const Stats &stats, double factor) override {
    bucket trSZ, trSZ2; // Tr[S_z], Tr[S_z^2]

    for (const auto &[I, eig]: diag) {
      const SZspin ssz  = I.get("SSZ");
      const double sumZ = calculate_Z(I, eig, factor);

      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * sqr(SZ(ssz));
    }

    Sz2 = trSZ2 / stats.Z;
    Sz  = trSZ / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_GLOBAL;
  HAS_SUBSTEPS;
};

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(step, i, j, ch, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

void SymmetrySPU1::make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  if (!P.substeps) {
    switch (P.channels) {
      case 1:
#include "spu1/spu1-1ch-offdiag.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag.dat"
#include "spu1/spu1-1ch-isospinx.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-diag.dat"
#include "spu1/spu1-2ch-offdiag.dat"
#include "spu1/spu1-2ch-anomalous.dat"
#include "spu1/spu1-2ch-isospinx.dat"
        break;
      default: my_assert_not_reached();
    }
  } else {
    my_assert(P.coeffactor == 1);
    const auto [Ntrue, M] = step.NM();

// XXX: note the (M == 1 ? -1 : 1) term.

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(step, i, j, M, t_matel(factor) * 2.0 * (M == 1 ? -1.0 : 1.0) * coef.delta(Ntrue + 1, M), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, M, 0, t_matel(factor) * coef.kappa(Ntrue, M), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, M, 0, t_matel(factor0) * coef.xi(Ntrue, M), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(step, i, M, number, coef.zeta(Ntrue + 1, M), h, qq)

#include "spu1/spu1-1ch-offdiag.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag.dat"
#include "spu1/spu1-1ch-isospinx.dat"
  }
}

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(step, i, j, ch, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG_UP
#undef OFFDIAG_DOWN
#undef DIAG_UP
#undef DIAG_DOWN

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DOWN(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiDOWN(step.N(), ch), h, qq, In, opch)

#define DIAG_UP(i, j, ch, number) diag_function_half(step, i, ch, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) diag_function_half(step, i, ch, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

void SymmetrySPU1::make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  my_assert(!P.substeps);
  switch (P.channels) {
  case 1:
#include "spu1/spu1-1ch-offdiag-UP.dat"
#include "spu1/spu1-1ch-offdiag-DOWN.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag-UP.dat"
#include "spu1/spu1-1ch-diag-DOWN.dat"
#include "spu1/spu1-1ch-isospinx.dat"
    break;
  case 2:
#include "spu1/spu1-2ch-diag-UP.dat"
#include "spu1/spu1-2ch-diag-DOWN.dat"
#include "spu1/spu1-2ch-offdiag-UP.dat"
#include "spu1/spu1-2ch-offdiag-DOWN.dat"
#include "spu1/spu1-2ch-anomalous.dat"
#include "spu1/spu1-2ch-isospinx.dat"
    break;
  default: my_assert_not_reached();
  }
}

void SymmetrySPU1::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  if (P.polarized) {
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  } else {
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
  }
}

#include "nrg-recalc-SPU1.cc"
