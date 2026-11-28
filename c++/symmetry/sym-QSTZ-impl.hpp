#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "doublet_transfer.hpp"

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetryQSTZ : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSTZ(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), td::square("<Tz^2>", "TZ"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1,0), Invar(1,2,1)) {
     initInvar({
        {"Q", additive},  // charge
        {"SS", additive}, // spin
        {"TZ", additive}  // angular momentum
     });
   }

  // Multiplicity of the (Q,SS,TZ) subspace is (2S+1 = SS).
  size_t mult(const Invar &I) const override { return I.get("SS"); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && u1_equality(I1.get("TZ"), I2.get("TZ"), I3.get("TZ"));
  }

  bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

  void load() override {
    my_assert(!P.substeps);
    my_assert(P.channels == 3);
#include "qstz/qstz-In2.dat"
#include "qstz/qstz-QN.dat"
  } // load

  // Same as for SYMTYPE=QS, because spin operators are angular momentum singlets.
  double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "TZ", 0);
    return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_diff(Ip, I1, "Q", 1);
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
  }

  DECL;
  HAS_DOUBLET;
  MatrixElements<SC> recalc_doublet_masked(const DiagInfo<SC> &diag, const SubspaceStructure &substruct,
                                          const MatrixElements<SC> &cold, unsigned mask) const override;
  // The transfers of TZ (QSTZ) present in the seed of a doublet, see doublet_transfer.hpp
  std::optional<unsigned> doublet_transfer_mask(const std::string &name, const MatrixElements<SC> &seed) const override {
    return qstz_transfer::infer(name, seed);
  }
  HAS_TRIPLET;
};

// We take the coefficients of the first channel (indexed as 0),
// because all three set are exactly the same due to orbital
// symmetry. [XXXX True for QST, for QSTZ we could relax this approx.]

template<typename SC>
void SymmetryQSTZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ss = I.get("SS");
  my_assert(!P.substeps);
  my_assert(P.channels == 3);
#include "qstz/qstz-offdiag.dat"
#include "qstz/qstz-diag.dat"
}

}

#include "nrg-recalc-QSTZ.hpp"
