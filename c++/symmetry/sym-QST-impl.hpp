#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetryQST : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQST(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), TDTrace{"<Tz^2>", [](const Invar &I) { const auto t = I.get("T"); return t * (t + 1) / 3.; }}, td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1,0), Invar(1,2,1)) {
     initInvar({
        {"Q", additive},  // charge
        {"SS", additive}, // spin
        {"T", additive}   // angular momentum
     });
   }

  // Multiplicity of the (Q,SS,T) subspace is (2S+1 = SS) times (2T+1).
  size_t mult(const Invar &I) const override { return I.get("SS") * (2 * I.get("T") + 1); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
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
#include "qst/qst-In2.dat"
#include "qst/qst-QN.dat"
  } // load

  // QS spin trace times the multiplicity of the spectator orbital multiplet.
  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "T", 0);
    return (2 * Ip.get("T") + 1) * su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  double dynamic_orb_susceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "SS", 0);
    return su2_triplet_factor(2 * Ip.get("T") + 1, 2 * I1.get("T") + 1);
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  // See clebsch_gordan_qst.nb
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS")) * orbital_triplet_factor(2 * Ip.get("T") + 1, 2 * I1.get("T") + 1);
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_ORB_TRIPLET;

  bool recalc_f_coupled(const Invar &I1, const Invar &I2, const Invar &If) const override {
    return triangle_inequality(I1, I2, If);
  }
};

inline bool qst_exception(const unsigned int i, const unsigned int j, const Invar &I) {
  // In these cases the subspace exists, but taking the T=2 or T=1
  // limit shows that the coefficient is actually zero, so there is
  // no contribution. (Directly computed factor is nan.)
  // This exception handling is added to avoid false positives
  // in error detection assertions.
  int T = I.get("T");
  if (i == 9 && j == 34 && T == 2) return true;
  if (i == 17 && j == 23 && T == 1) return true;
  if (i == 52 && j == 58 && T == 1) return true;
  return false;
}

// The table entries OFFDIAG(i, j, factor) skip the exceptional blocks and the blocks with zero factor. We take the
// coefficients of the first channel (indexed as 0), because all three sets are exactly the same due to orbital
// symmetry.
template<typename SC>
class QSTTerms : public ShellTerms<SC> {
   const Invar &I;
 public:
   template<typename ... Args> QSTTerms(const Invar &I, Args && ... args) : ShellTerms<SC>(std::forward<Args>(args)...), I(I) {}
   [[nodiscard]] bool contributes_offdiag(const size_t i, const size_t j) const { return this->contributes(i, j) && !qst_exception(i, j, I); }
   void offdiag(const size_t i, const size_t j, const coef_traits<SC> factor) const {
     const auto f = factor * this->coef.xi(this->N, 0);
     if (std::abs(f) > 0.0) this->hop(i, j, 0, 0, f);
   }
};

template<typename SC>
void SymmetryQST<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const QSTTerms<SC> terms(I, *this, P, step, h, qq, In, opch, coef);
  int ss  = I.get("SS");
  int t = I.get("T");
  double T  = t; // crucially important to use floating point!
  my_assert(!P.substeps);
  my_assert(P.channels == 3);
#include "qst/qst-offdiag.dat"
#include "qst/qst-diag.dat"
}

}

#include "nrg-recalc-QST.hpp"
