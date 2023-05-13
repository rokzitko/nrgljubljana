namespace NRG {

template<typename SC>
class SymmetrySPSU2 : public Symmetry<SC> {
  private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPSU2(const Params &P) : Symmetry<SC>(P, std::vector{"<Sz^2>"}, Invar(1)) {
     initInvar({
        {"SS", additive} // spin
     });
   }

  // Multiplicity of the I=(SS) subspace = 2S+1 = SS.
  size_t mult(const Invar &I) const override { return I.get("SS"); }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"));
  }

  void load() override {
    if (!P.substeps) {
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-In2.dat"
#include "spsu2/spsu2-1ch-QN.dat"
          break;

        case 2:
#include "spsu2/spsu2-2ch-In2.dat"
#include "spsu2/spsu2-2ch-QN.dat"
          break;

        case 3:
#include "spsu2/spsu2-3ch-In2.dat"
#include "spsu2/spsu2-3ch-QN.dat"
          break;

        default: my_assert_not_reached();
      }
    } else {
#include "spsu2/spsu2-1ch-In2.dat"
#include "spsu2/spsu2-1ch-QN.dat"
    }
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert(abs(ss1 - ssp) == 2 || ss1 == ssp);
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert(abs(ss1 - ssp) == 1);
    return (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ; // Tr[S_z^2]
    for (const auto &[I, eig]: diag) {
      const int ss    = I.get("SS");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * (ss * ss - 1) / 12.;
    }
    stats.td.set("<Sz^2>", trSZ / stats.Z);
  }

  bool project_subspace(const Invar &I, const std::string &p) const override {
    if (p == ""s || p == "trivial"s) {
      return true;
    } else if (p == "even"s) {  // even particule number is half-interger spin, SS=1,3,5,...
      return is_odd(I.get("SS"));
    } else if (p == "odd"s) {   // odd particle number is integer spin, SS=2,4,...
      return is_even(I.get("SS"));
    } else throw std::runtime_error(fmt::format("Unknown projection type {} for symmetry SPSU2.", p));
  }

  void show_coefficients(const Step &step, const Coef<SC> &coef) const override {
    Symmetry<SC>::show_coefficients(step, coef);
    if (!P.substeps) {
      for (size_t i = 0; i < P.coefchannels; i++) {
        const auto N = step.N();
        std::cout << "[" << i + 1 << "]"
          << " Delta(" << N+1 << ")=" << coef.delta(N+1, i) << " kappa(" << N << ")=" << coef.kappa(N, i) << std::endl;
      }
    } else {
      const auto [N, M] = step.NM();
      for (auto i = 0; i < P.coeffactor; i++) {
        const auto index = M + P.channels * i;
        std::cout << "[" << index << "]"
          << " Delta(" << N+1 << ")=" << coef.delta(N+1, index) << " kappa(" << N << ")=" << coef.kappa(N, index) << std::endl;
      }
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_GLOBAL;
  HAS_SUBSTEPS;
};

// NOTE: the additional factor of 2 for coef.delta is due to the fact that the isospin is in fact defined as 1/2 of
// the (d^\dag d^\dag + d d) pairing operator.

template<typename SC>
void SymmetrySPSU2<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int ss = I.get("SS");

  if (!P.substeps) {
// conjugate the coefficient, to bring it in line with the standard convention for superconducting phases (only relevant for complex-valued parameters)
#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * conj_me( coef.delta(step.N() + 1, ch) ), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

    switch (P.channels) {
      case 1:
#include "spsu2/spsu2-1ch-offdiag.dat"
#include "spsu2/spsu2-1ch-anomalous.dat"
#include "spsu2/spsu2-1ch-diag.dat"
#include "spsu2/spsu2-1ch-isospinx.dat"
        break;

      case 2:
#include "spsu2/spsu2-2ch-diag.dat"
#include "spsu2/spsu2-2ch-offdiag.dat"
#include "spsu2/spsu2-2ch-anomalous.dat"
#include "spsu2/spsu2-2ch-isospinx.dat"
        break;

      case 3:
#include "spsu2/spsu2-3ch-diag.dat"
#include "spsu2/spsu2-3ch-offdiag.dat"
#include "spsu2/spsu2-3ch-anomalous.dat"
#include "spsu2/spsu2-3ch-isospinx.dat"
        break;

      default: my_assert_not_reached();
    }
  } else {
    my_assert(P.coeffactor == 1);
    const auto [Ntrue, M] = step.NM();

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * conj_me( coef.delta(Ntrue + 1, M) ), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, M, 0, t_matel(factor) * coef.kappa(Ntrue, M), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, M, 0, t_matel(factor0) * coef.xi(Ntrue, M), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(Ntrue + 1, M), h, qq)

#include "spsu2/spsu2-1ch-offdiag.dat"
#include "spsu2/spsu2-1ch-anomalous.dat"
#include "spsu2/spsu2-1ch-diag.dat"
#include "spsu2/spsu2-1ch-isospinx.dat"
  }
}

}

#include "nrg-recalc-SPSU2.hpp"

