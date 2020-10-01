class SymmetrySPU1 : public SymField {
  private:
  outfield Sz2, Sz;

  public:
  SymmetrySPU1() : SymField() { all_syms["SPU1"] = this; }

  void init() override {
    Sz2.set("<Sz^2>", 1);
    Sz.set("<Sz>", 2);
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
      switch (channels) {
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

  void makematrix_polarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);
  void makematrix_nonpolarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);

  void calculate_TD(const DiagInfo &diag, double factor) override {
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

Symmetry *SymSPU1 = new SymmetrySPU1;

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(step.N() + 1, ch), h, qq)

void SymmetrySPU1::makematrix_nonpolarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (!P.substeps) {
    switch (channels) {
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
    int Ntrue, M;
    tie(Ntrue, M) = get_Ntrue_M(step.N());

// See above for the explanation for the (M == 1 ? -1 : 1) term.
#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, M, t_matel(factor) * 2.0 * (M == 1 ? -1.0 : 1.0) * delta(Ntrue + 1, M), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, M, 0, t_matel(factor) * kappa(Ntrue, M), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, M, 0, t_matel(factor0) * xi(Ntrue, M) / scale_fix(step.N()), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, M, number, zeta(Ntrue + 1, M), h, qq)

#include "spu1/spu1-1ch-offdiag.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag.dat"
#include "spu1/spu1-1ch-isospinx.dat"
  }
}

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG_UP
#undef OFFDIAG_DOWN
#undef DIAG_UP
#undef DIAG_DOWN

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DOWN(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xiDOWN(step.N(), ch), h, qq, In, opch)

#define DIAG_UP(i, j, ch, number) diag_function_half(i, ch, number, zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) diag_function_half(i, ch, number, zetaDOWN(step.N() + 1, ch), h, qq)

void SymmetrySPU1::makematrix_polarized(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (!P.substeps) {
    switch (channels) {
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
    } // switch
  } else {
    my_error("Not implemented.");
  }
}

void SymmetrySPU1::makematrix(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (P.polarized) {
    makematrix_polarized(h, qq, I, In, opch);
  } else {
    makematrix_nonpolarized(h, qq, I, In, opch);
  }
}

#include "nrg-recalc-SPU1.cc"
