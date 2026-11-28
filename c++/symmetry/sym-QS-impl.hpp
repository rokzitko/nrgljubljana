#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;
using std::string_literals::operator""s;

template<typename SC>
class SymmetryQS : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQS(const Params &P) : Symmetry<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}, Invar(0,1), Invar(1,2)) {
     initInvar({
       {"Q", additive}, // charge
       {"SS", additive} // spin
     });
   }

   // Multiplicity of the (Q,SS) subspace is 2S+1 = SS.
   size_t mult(const Invar &I) const override { return I.get("SS"); }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"));
   }

   bool Invar_allowed(const Invar &I) const override { return I.get("SS") > 0; }

   void load() override {
     if (!P.substeps) {
      switch (P.channels) {
        case 1:
#include "qs/qs-1ch-In2.dat"
#include "qs/qs-1ch-QN.dat"
          break;
        case 2:
#include "qs/qs-2ch-In2.dat"
#include "qs/qs-2ch-QN.dat"
          break;
        case 3:
#include "qs/qs-3ch-In2.dat"
#include "qs/qs-3ch-QN.dat"
          break;
        case 4:
#include "qs/qs-4ch-In2.dat"
#include "qs/qs-4ch-QN.dat"
          break;
        default: my_assert_not_reached();
      }
    } else {
#include "qs/qs-1ch-In2.dat"
#include "qs/qs-1ch-QN.dat"
    }
   }

   double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const override {
     check_diff(Ip, I1, "Q", 0);
     return su2_triplet_factor(Ip.get("SS"), I1.get("SS"));
   }

   double specdens_factor(const Invar &Ip, const Invar &I1) const override {
     check_diff(Ip, I1, "Q", 1);
     return su2_doublet_factor(Ip.get("SS"), I1.get("SS"));
   }

   bool project_subspace(const Invar &I, const std::string &p) const override {
     if (p == ""s || p == "trivial"s) {
       return true;
     } else if (p == "evenQ"s) {
       return is_even(I.get("Q"));
     } else if (p == "oddQ"s) {
       return is_odd(I.get("Q"));
     } else throw std::runtime_error(fmt::format("Unknown projection type {} for symmetry QS.", p));
   }

   DECL;
   HAS_DOUBLET;
   HAS_TRIPLET;
   HAS_GLOBAL;
   HAS_SUBSTEPS;
   void show_coefficients(const Step &, const Coef<SC> &) const override;
};


/* i - subspace index
   ch - channel (0 or 1)
   number - number of electrons added in channel 'ch' in subspace 'i' */

template<typename SC>
void SymmetryQS<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, 
                                     const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  auto ss = I.get("SS");

  if (!P.substeps) {
    switch (P.channels) {
      case 1:
#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"
        break;
      case 2:
#include "qs/qs-2ch-diag.dat"
#include "qs/qs-2ch-offdiag.dat"
        if (P.rungs) {
#include "qs/qs-2ch-offdiag-mix.dat"
#include "qs/qs-2ch-runghop.dat"
        }
        break;
      case 3:
#include "qs/qs-3ch-diag.dat"
#include "qs/qs-3ch-offdiag.dat"
        break;
      case 4:
#include "qs/qs-4ch-diag.dat"
#include "qs/qs-4ch-offdiag.dat"
        break;
      default: my_assert_not_reached();
    }
  } else { // one site of channel M, see ShellTerms
    my_assert(P.coeffactor == 1);
    if (P.rungs) my_assert_not_reached();
#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"
  }
}

template<typename SC>
void SymmetryQS<SC>::show_coefficients(const Step &step, const Coef<SC> &coef) const {
  Symmetry<SC>::show_coefficients(step, coef);
  if (P.silent) return;
  if (P.rungs)
    for (auto i = 0; i < P.channels; i++)
      std::cout << "[" << i + 1 << "]"
           << " xi_rung(" << step.N() << ")=" << coef.xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" 
           << coef.zetaR(step.N() + 1, i) << std::endl;
}

}

#include "nrg-recalc-QS.hpp"
