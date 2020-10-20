template<typename SC>
class SymmetrySU2_tmpl : public Symmetry_tmpl<SC> {
 private:
   outfield Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySU2_tmpl(const Params &P, Allfields &allfields) : Symmetry_tmpl<SC>(P),
     Q2(P, allfields, "<Q^2>", 1) {
       initInvar({
         {"II", additive} // isospin
       });
       this->InvarSinglet = Invar(1);
     }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"));
  }

  // Multiplicity of the I=(II) subspace = (2I+1) = II.
  size_t mult(const Invar &I) const override {
    return I.get("II"); // isospin multiplicity
  }

  // We always must have I >= 0.
  bool Invar_allowed(const Invar &I) const override { return I.get("II") > 0; }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");
    my_assert(abs(ii1 - iip) == 1);
    const double isofactor = (ii1 == iip + 1 ? ISO(iip) + 1.0 : ISO(iip));
    return isofactor;
  }

  void load() override {
    switch (P.channels) {
      case 1:
#include "su2/su2-1ch-In2.dat"
#include "su2/su2-1ch-QN.dat"
        break;
      case 2:
#include "su2/su2-2ch-In2.dat"
#include "su2/su2-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats_tmpl<SC> &stats, const double factor) override {
    bucket trIZ2; // Tr[I_z^2]
    for (const auto &[I, eig]: diag) {
      const Number ii   = I.get("II");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trIZ2 += sumZ * (ii * ii - 1) / 12.;
    }
    Q2 = (4 * trIZ2) / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

// For SU2, the <||f||> depend on the "type"!
// OFFDIAG_1 corresponds to the first combination, OFFDIAG_2 to the
// second combination of operators. Each of them has contributions
// in each channel.

#undef OFFDIAG_1
#define OFFDIAG_1(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)
#undef OFFDIAG_2
#define OFFDIAG_2(i, j, ch, factor) offdiag_function(step, i, j, ch, 1, t_matel(factor) * coef.xi(step.N(), ch), h, qq, In, opch)

template<typename SC>
void SymmetrySU2_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  Ispin ii = I.get("II");
  int NN   = step.getnn();
  switch (P.channels) {
    case 1:
#include "su2/su2-1ch-offdiag-1.dat"
#include "su2/su2-1ch-offdiag-2.dat"
      break;
    case 2:
#include "su2/su2-2ch-offdiag-1.dat"
#include "su2/su2-2ch-offdiag-2.dat"
      break;
    default: my_assert_not_reached();
  }
}

#include "nrg-recalc-SU2.cc"
