namespace NRG {

template<typename SC>
class SymmetryQSLR : public SymLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSLR(const Params &P) : SymLR<SC>(P, std::vector{"<Sz^2>", "<Q>", "<Q^2>"}, Invar(0,1,1)) {
     initInvar({
        {"Q", additive},      // charge
        {"SS", additive},     // spin
        {"P", multiplicative} // parity
     });
   }

  // Multiplicity of the I=(Q,SS,P) subspace = 2S+1 = SS.
  size_t mult(const Invar &I) const override {
    return I.get("SS"); // spin multiplicity
  }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "qslr/qslr-2ch-In2.dat"
#include "qslr/qslr-2ch-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert((abs(ss1 - ssp) == 2 || ss1 == ssp));
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    return (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trQ, trQ2; // Tr[S_z^2], Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const int ss    = I.get("SS");
      const int q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
      trSZ += sumZ * (ss * ss - 1) / 12.;
    }
    stats.td.set("<Sz^2>", trSZ / stats.Z);
    stats.td.set("<Q>",    trQ / stats.Z);
    stats.td.set("<Q^2>",  trQ2 / stats.Z);
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)
#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryQSLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int ss = I.get("SS");
#include "qslr/qslr-2ch-diag.dat"
#include "qslr/qslr-2ch-offdiag.dat"
}

}

#include "nrg-recalc-QSLR.hpp"
