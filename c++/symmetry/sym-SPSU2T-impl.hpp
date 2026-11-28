#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetrySPSU2T : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPSU2T(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), TDTrace{"<Tz^2>", [](const Invar &I) { const auto t = I.get("T"); return t * (t + 1) / 3.; }}}, Invar(1,0), Invar(2,1)) {
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

  // QS spin trace times the multiplicity of the spectator orbital multiplet.
  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "T", 0);
    return (2 * Ip.get("T") + 1) * su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS")) * orbital_triplet_factor(2 * Ip.get("T") + 1, 2 * I1.get("T") + 1);
  }

  DECL;
  bool uses_pairing() const override { return true; }
  bool uses_anomalous_hopping() const override { return true; }
  bool complex_pairing_supported() const override { return false; }
  HAS_DOUBLET;
  HAS_TRIPLET;

  bool recalc_f_coupled(const Invar &I1, const Invar &I2, const Invar &If) const override {
    return triangle_inequality(I1, I2, If);
  }
};

inline bool spsu2t_exception(const unsigned int i, const unsigned int j, const Invar &I) {
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

// The table entries OFFDIAG(i, j, factor) skip the exceptional blocks. We take the coefficients of the first channel
// (indexed as 0), because all three sets are exactly the same due to orbital symmetry.
//
// In the SPSU2T isospinx table the lower-indexed state of a pair does not always have fewer electrons (unlike in the
// other superconducting symmetry types), so a single phase convention for complex delta cannot be applied to all
// entries. Only real delta is supported (complex_pairing_supported()), and delta is not conjugated.
template<typename SC>
class SPSU2TTerms : public ShellTerms<SC> {
   const Invar &I;
 public:
   template<typename ... Args> SPSU2TTerms(const Invar &I, Args && ... args) : ShellTerms<SC>(std::forward<Args>(args)..., false), I(I) {}
   [[nodiscard]] bool contributes_offdiag(const size_t i, const size_t j) const { return this->contributes(i, j) && !spsu2t_exception(i, j, I); }
};

template<typename SC>
void SymmetrySPSU2T<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const SPSU2TTerms<SC> terms(I, *this, P, step, h, qq, In, opch, coef);
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
