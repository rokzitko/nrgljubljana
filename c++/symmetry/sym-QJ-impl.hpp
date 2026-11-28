#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include <boost/math/special_functions/factorials.hpp>

#include "table_terms.hpp"

namespace NRG {

using std::abs;
using std::round;
using std::sqrt;

template<typename SC>
class SymmetryQJ : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQJ(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Jz^2>", "JJ"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1), Invar(1,4)) {
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
    my_assert(!P.substeps);
    my_assert(P.channels == 3); // three p-orbitals (64 combinations)
#include "qj/qj-In2.dat"
#include "qj/qj-QN.dat"
  }

  // ClebschGordan[ket (p), op, bra (1)]
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    return su2_doublet_factor(Ip.get("JJ"), I1.get("JJ"));
  }

  // See cg_factors_doublet_triplet_quadruplet.nb
  double specdensquad_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    const int jjp = Ip.get("JJ");
    const int jj1 = I1.get("JJ");
    my_assert(std::abs(jj1 - jjp) == 1 || std::abs(jj1 - jjp) == 3);
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

  HAS_DOUBLET;
  HAS_QUADRUPLET;
  DECL;
};

double Factorial(const double x) { return boost::math::factorial<double>(round(x)); }

// The table entries OFFDIAG(i, j, Jndx, factor) refer to the doublet (Jndx = 0) or quadruplet (Jndx = 1) <||f||>,
// stored as channels 0 and 1. The hopping only contributes if the triangle inequality is satisfied for the ancestor
// subspaces In[i-1], In[j-1] of the combinations i, j (1-based), whose <||f||> block is used.
template<typename SC>
class QJTerms : public ShellTerms<SC> {
   const SymmetryQJ<SC> &qj;
 public:
   template<typename ... Args> QJTerms(const SymmetryQJ<SC> &qj, Args && ... args) : ShellTerms<SC>(qj, std::forward<Args>(args)...), qj(qj) {}
   void offdiag(const size_t i, const size_t j, const size_t jndx, const coef_traits<SC> factor) const {
     const Invar Iop = jndx == 0 ? Invar(1, 2) : Invar(1, 4);
     if (qj.triangle_inequality(this->In[i-1], this->In[j-1], Iop)) this->hop(i, j, jndx, 0, factor * this->coef.xi(this->N, 0)); // I1 = I2+Iop
   }
};

inline double J(int JJ) {
  return (JJ - 1.0) / 2.0; // JJ=2J+1
}

template<typename SC>
void SymmetryQJ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const QJTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int jj = I.get("JJ");
#include "qj/qj-offdiag.dat"
#include "qj/qj-diag.dat"
}

}

#include "nrg-recalc-QJ.hpp"
