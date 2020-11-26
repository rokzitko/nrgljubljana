#include <boost/math/special_functions/factorials.hpp>

namespace NRG {

template<typename SC>
class SymmetryQJ : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQJ(const Params &P) : Symmetry<SC>(P, std::vector{"<Jz^2>", "<Q>", "<Q^2>"}, Invar(0,1), Invar(1,4)) {
     initInvar({
        {"Q", additive}, // charge
        {"JJ", additive} // total angular momentum
     });
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

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trJZ2, trQ, trQ2; // Tr[J_z^2], Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const int jj    = I.get("JJ");
      const int q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
      trJZ2 += sumZ * (jj * jj - 1) / 12.;
    }
    stats.td.set("<Jz^2>", trJZ2 / stats.Z);
    stats.td.set("<Q>",    trQ / stats.Z);
    stats.td.set("<Q^2>",  trQ2 / stats.Z);
  }

  // ClebschGordan[ket (p), op, bra (1)]
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const int jjp = Ip.get("JJ");
    const int jj1 = I1.get("JJ");
    my_assert(abs(jj1 - jjp) == 1);
    return (jj1 == jjp + 1 ? S(jjp) + 1.0 : S(jjp));
  }

  // See cg_factors_doublet_triplet_quadruplet.nb
  double specdensquad_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const int jjp = Ip.get("JJ");
    const int jj1 = I1.get("JJ");
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

   void offdiag_function_QJ(const Step &step, const unsigned int i, const unsigned int j, const unsigned int ch, const unsigned int fnr, const t_matel factor, Matrix &h, const SubspaceDimensions &qq,
                            const InvarVec &In, const Opch<SC> &opch) const
   {
     const Invar Iop     = ch == 0 ? Invar(1, 2) : Invar(1, 4);
     const Invar I1      = In[i];
     const Invar I2      = In[j];
     const bool triangle = triangle_inequality(I1, I2, Iop); // I1 = I2+Iop
     if (triangle) { offdiag_function(step, i, j, ch, fnr, factor, h, qq, In, opch); }
   }

  HAS_DOUBLET;
  HAS_QUADRUPLET;
  DECL;
};

double Factorial(const double x) { return boost::math::factorial<double>(round(x)); }

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
void SymmetryQJ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int jj = I.get("JJ");
#include "qj/qj-offdiag.dat"
#include "qj/qj-diag.dat"
}

}

#include "nrg-recalc-QJ.hpp"
