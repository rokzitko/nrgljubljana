namespace NRG {

template<typename SC>
class SymmetrySPSU2T : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPSU2T(const Params &P) : Symmetry<SC>(P, std::vector{"<Sz^2>", "<Tz^2>"}, Invar(1,0), Invar(2,1)) {
       initInvar({
         {"SS", additive}, // spin
         {"T", additive}   // angular momentum
       });
     }

  // Multiplicity of the (SS,T) subspace is (2S+1 = SS) times (2T+1).
  size_t mult(const Invar &I) const override { return I.get("SS") * (2 * I.get("T") + 1); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && su2_triangle_inequality(2 * I1.get("T") + 1, 2 * I2.get("T") + 1, 2 * I3.get("T") + 1);
  }

  bool Invar_allowed(const Invar &I) const override {
    const bool spin_ok   = I.get("SS") > 0;
    const bool angmom_ok = I.get("T") >= 0;
    return spin_ok && angmom_ok;
  }

  void load() override {
    my_assert(!P.substeps);
    my_assert(P.channels == 3);
#include "spsu2t/spsu2t-In2.dat"
#include "spsu2t/spsu2t-QN.dat"
  } // load

  // Same as for SYMTYPE=QS, because spin operators are angular momentum singlets.
  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "T", 0);
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert((abs(ss1 - ssp) == 2 || ss1 == ssp));
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert(abs(ss1 - ssp) == 1);
    double spinfactor = (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
    const int tp = Ip.get("T");
    const int t1 = I1.get("T");
    const int ttp    = 2 * tp + 1;
    const int tt1    = 2 * t1 + 1;
    my_assert(abs(ttp - tt1) == 2 || ttp == tt1);
    double angmomfactor = switch3(tt1, ttp + 2, 1. + (ttp - 1) / 3., ttp, ttp / 3., ttp - 2, (-2. + ttp) / 3.);
    return spinfactor * angmomfactor;
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ2, trTZ2; // Tr[S_z^2], Tr[T_z^2]
    for (const auto &[I, eig]: diag) {
      const int ss    = I.get("SS");
      const int t   = I.get("T");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ2 += sumZ * (ss * ss - 1) / 12.; // [(2S+1)(2S+1)-1]/12=S(S+1)/3
      trTZ2 += sumZ * t * (t + 1) / 3.;
    }
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Tz^2>", trTZ2 / stats.Z);
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;

  bool recalc_f_coupled(const Invar &I1, const Invar &I2, const Invar &If) const override {
    return triangle_inequality(I1, I2, If);
  }
};

bool spsu2t_exception(const unsigned int i, const unsigned int j, const Invar &I) {
  // In these cases the subspace exists, but taking the T=2 or T=1 limit shows that the coefficient is actually zero,
  // so there is no contribution. (Directly computed factor is nan.) This exception handling is added io order to
  // avoid false positives in error detection assertions.

  // see spsu2t_exceptions.nb
  int T = I.get("T");
  if (i == 46 && j == 21 && T == 1) return true;
  if (i == 55 && j == 27 && T == 1) return true;
  if (i == 5 && j == 38 && T == 2) return true;

  return false;
}

#define offdiag_spsu2t(i, j, ch, fnr, factor0, h, qq, In, I, opch)                                                                                   \
  {                                                                                                                                                  \
    const bool contributes = qq.offdiag_contributes(i, j);                                                                                           \
    if (contributes) {                                                                                                                               \
      t_matel factor;                                                                                                                                \
      if (spsu2t_exception(i, j, I)) {                                                                                                               \
        factor = 0;                                                                                                                                  \
      } else {                                                                                                                                       \
        factor = factor0;                                                                                                                            \
      }                                                                                                                                              \
      this->offdiag_function_impl(step, i, j, ch, fnr, factor, h, qq, In, opch);                                                                           \
    }                                                                                                                                                \
  };

// We take the coefficients of the first channel (indexed as 0), because all three set are exactly the same due to
// orbital symmetry.
#undef OFFDIAG
#define OFFDIAG(i, j, factor0) offdiag_spsu2t(i, j, 0, 0, t_matel(factor0) * coef.xi(step.N(), 0), h, qq, In, I, opch)

#undef DIAG
#define DIAG(i, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, 0), h, qq)

#undef ISOSPINX
#define ISOSPINX(i, j, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, 0), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, factor) offdiag_function(step, i, j, 0, 0, t_matel(factor) * coef.kappa(step.N(), 0), h, qq, In, opch)

template<typename SC>
void SymmetrySPSU2T<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int ss  = I.get("SS");
  int t = I.get("T");
  double T  = t; // crucially important to use floating point!
  my_assert(!P.substeps);
  my_assert(P.channels == 3);
#include "spsu2t/spsu2t-offdiag.dat"
#include "spsu2t/spsu2t-diag.dat"
  if (!num_equal(coef.delta(step.N() + 1, 0), 0.0)) {
#include "spsu2t/spsu2t-isospinx.dat"
  }
  if (!num_equal(coef.kappa(step.N(), 0), 0.0)) {
#include "spsu2t/spsu2t-anomalous.dat"
  }
}

}

#include "nrg-recalc-SPSU2T.hpp"
