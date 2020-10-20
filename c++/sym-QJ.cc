template<typename SC>
class SymmetryQJ_tmpl : public Symmetry_tmpl<SC> {
 private:
   outfield Jz2, Q, Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQJ_tmpl(const Params &P, Allfields &allfields) : Symmetry_tmpl<SC>(P),
     Jz2(P, allfields, "<Jz^2>", 1), Q(P, allfields, "<Q>", 2), Q2(P, allfields, "<Q^2>", 3) {
       initInvar({
         {"Q", additive}, // charge
         {"JJ", additive} // total angular momentum
       });
       this->InvarSinglet = Invar(0, 1);
       this->Invar_f      = Invar(1, 4);
     }

  // Multiplicity of the (Q,JJ) subspace is 2J+1 = JJ.
  size_t mult(const Invar &I) const override { return I.get("JJ"); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("JJ"), I2.get("JJ"), I3.get("JJ"));
  }

  bool Invar_allowed(const Invar &I) const override {
    const bool spin_ok = I.get("JJ") > 0;
    return spin_ok;
  }

  void load() override {
#include "qj/qj-In2.dat"
#include "qj/qj-QN.dat"
  }

  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats &stats, const double factor) override {
    bucket trJZ2, trQ, trQ2; // Tr[J_z^2], Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const Sspin jj    = I.get("JJ");
      const Number q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
      trJZ2 += sumZ * (jj * jj - 1) / 12.;
    }
    Jz2 = trJZ2 / stats.Z;
    Q   = trQ / stats.Z;
    Q2  = trQ2 / stats.Z;
  }

  // ClebschGordan[ket (p), op, bra (1)]
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const Sspin jjp = Ip.get("JJ");
    const Sspin jj1 = I1.get("JJ");
    my_assert(abs(jj1 - jjp) == 1);
    return (jj1 == jjp + 1 ? S(jjp) + 1.0 : S(jjp));
  }

  // See cg_factors_doublet_triplet_quadruplet.nb
  double specdensquad_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const Sspin jjp = Ip.get("JJ");
    const Sspin jj1 = I1.get("JJ");
    my_assert(abs(jj1 - jjp) == 1 || abs(jj1 - jjp) == 3);
    if (jj1 == jjp + 3) return S(jjp) / 2.0 + 1.0;
    if (jj1 == jjp + 1) {
      my_assert(jjp >= 2); // singlet in kvadruplet se ne moreta sklopiti v dublet
      return S(jjp) / 2.0 + 0.5;
    }
    if (jj1 == jjp - 1) {
      my_assert(jjp >= 2); // trikotniska
      return S(jjp) / 2.0;
    }
    if (jj1 == jjp - 3) {
      my_assert(jjp >= 4); // trikotniska
      return S(jjp) / 2.0 - 0.5;
    }
    my_assert_not_reached();
    return 0;
  }

  // AAA: const
  void offdiag_function_QJ(const Step &step, unsigned int i, unsigned int j, unsigned int ch, unsigned int fnr, t_matel factor, Matrix &h, const Rmaxvals &qq,
                           const InvarVec &In, const Opch_tmpl<SC> &opch);

  HAS_DOUBLET;
  HAS_QUADRUPLET;
  DECL;
};

#include <boost/math/special_functions/factorials.hpp>

double Factorial(const double x) { return boost::math::factorial<double>(round(x)); }

template<typename SC>
void SymmetryQJ_tmpl<SC>::offdiag_function_QJ(const Step &step, unsigned int i, unsigned int j,
                                     unsigned int ch,  // channel number
                                     unsigned int fnr, // extra index for <||f||>, usually 0
                                     t_matel factor,   // may be complex (in principle)
                                     Matrix &h, const Rmaxvals &qq, const InvarVec &In, const Opch_tmpl<SC> &opch) {
  const Invar Iop     = (ch == 0 ? Invar(1, 2) : Invar(1, 4));
  const Invar I1      = In[i];
  const Invar I2      = In[j];
  const bool triangle = triangle_inequality(I1, I2, Iop); // I1 = I2+Iop
  if (triangle) { offdiag_function(step, i, j, ch, fnr, factor, h, qq, In, opch); }
}

// *** Helper macros for make_matrix() members in matrix.cc
// Jndx = 0 for doublet, Jndx = 1 for quadruplet
#undef OFFDIAG
#define OFFDIAG(i, j, Jndx, factor0) offdiag_function_QJ(step, i, j, Jndx, 0, t_matel(factor0) * coef.xi(step.N(), 0), h, qq, In, opch)

#undef DIAG
#define DIAG(i, number) this->diag_function(step, i, 0, number, coef.zeta(step.N() + 1, 0), h, qq)

inline double J(int JJ) {
  return (JJ - 1.0) / 2.0; // JJ=2J+1
}

template<typename SC>
void SymmetryQJ_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  Sspin jj = I.get("JJ");
#include "qj/qj-offdiag.dat"
#include "qj/qj-diag.dat"
}

#include "nrg-recalc-QJ.cc"
