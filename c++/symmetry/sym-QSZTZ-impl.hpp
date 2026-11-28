#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
#include "doublet_transfer.hpp"

#include "table_terms.hpp"

namespace NRG {

using std::pow;
using std::sqrt;

template<typename SC>
class SymmetryQSZTZ : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSZTZ(const Params &P) : Symmetry<SC>(P, {td::sz("<Sz>", "SZ"), td::sz_square("<Sz^2>", "SZ"), td::mean("<Tz>", "TZ"), td::square("<Tz^2>", "TZ"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,0,0), Invar(1,2,1)) {
     initInvar({
        {"Q", additive},  // charge
        {"SZ", additive}, // spin
        {"TZ", additive}  // angular momentum
     });
   }
  size_t mult(const Invar &I) const override { return 1; }
  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SZ"), I2.get("SZ"), I3.get("SZ"))
       && u1_equality(I1.get("TZ"), I2.get("TZ"), I3.get("TZ"));
  }
  bool Invar_allowed(const Invar &I) const override { return true; }
  void load() override {
    my_assert(!P.substeps);
    my_assert(P.channels == 3);
#include "qsztz/qsztz-In2.dat"
#include "qsztz/qsztz-QN.dat"
  } // load
  DECL;
  HAS_DOUBLET;
  MatrixElements<SC> recalc_doublet_masked(const DiagInfo<SC> &diag, const SubspaceStructure &substruct,
                                                const MatrixElements<SC> &cold, unsigned mask) const override;
  // The transfers of TZ and SZ (QSZTZ) present in the seed of a doublet, see doublet_transfer.hpp
  std::optional<unsigned> doublet_transfer_mask(const std::string &name, const MatrixElements<SC> &seed) const override {
    return qsztz_transfer::infer(name, seed);
  }
  HAS_TRIPLET;
};

// We take the coefficients of the first channel (indexed as 0), because all three set are exactly the same due to
// orbital symmetry. [Can be generalized for symmetry-broken cases if the need arises.]

template<typename SC>
void SymmetryQSZTZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  my_assert(!P.substeps);
  my_assert(P.channels == 3);
#include "qsztz/qsztz-offdiag.dat"
#include "qsztz/qsztz-diag.dat"
}

}

#include "nrg-recalc-QSZTZ.hpp"
