class SymmetrySPSU2C3 : public SymC3 {
 private:
   outfield Sz2;

 public:
   template<typename ... Args> SymmetrySPSU2C3(Args&& ... args) : SymC3(std::forward<Args>(args)...),
     Sz2(P, allfields, "<Sz^2>", 1) {
       InvarStructure InvStruc[] = {
         {"SS", additive}, // spin
         {"P", mod3}       // C_3 rep
       };
       initInvar(InvStruc, ARRAYLENGTH(InvStruc));
       InvarSinglet = Invar(1, 0); // spin-singlet, C_3 P=0
     }

  size_t mult(const Invar &I) const override { return I.get("SS"); }

  bool Invar_allowed(const Invar &I) override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) && c3_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 3);
#include "spsu2c3/spsu2c3-In2.dat"
#include "spsu2c3/spsu2c3-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert((abs(ss1 - ssp) == 2 || ss1 == ssp));
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    return (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
  }

  void calculate_TD(const Step &step, const DiagInfo &diag, const Stats &stats, double factor) override {
    bucket trSZ2; // Tr[S_z^2]
    for (const auto &[I, eig]: diag) {
      const Sspin ss    = I.get("SS");
      const double sumZ = calculate_Z(I, eig, factor);
      trSZ2 += sumZ * (ss * ss - 1) / 12.;
    }
    Sz2 = trSZ2 / stats.Z;
  }

  DECL;
  // HAS_DOUBLET; HAS_TRIPLET;
};

#undef ISOSPINX
#define ISOSPINX(i, j, factor) diag_offdiag_function(step, i, j, 0, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, 0), h, qq)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, number) diag_function(step, i, 0, number, coef.zeta(step.N() + 1, 0), h, qq)

void SymmetrySPSU2C3::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) {
  if constexpr (std::is_same_v(t_matel, std::complex<double>)) {
    my_assert(P.channels == 3);
    Sspin ss = I.get("SS");
#undef Complex
#define Complex(x, y) cmpl(x, y)
#define sqrt(x) csqrt(x)
#include "spsu2c3/spsu2c3-diag.dat"
#include "spsu2c3/spsu2c3-offdiag.dat"
#include "spsu2c3/spsu2c3-isospinx.dat"
#undef sqrt
  } else my_assert_not_reached();
}

#include "nrg-recalc-SPSU2C3.cc"
