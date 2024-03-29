namespace NRG {

template<typename SC>
class SymmetryQS : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQS(const Params &P) : Symmetry<SC>(P, std::vector{"<Sz^2>", "<Q>", "<Q^2>"}, Invar(0,1), Invar(1,2)) {
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
     const auto ssp = Ip.get("SS");
     const auto ss1 = I1.get("SS");
     my_assert((abs(ss1 - ssp) == 2 || ss1 == ssp));
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
   }

   double specdens_factor(const Invar &Ip, const Invar &I1) const override {
     check_diff(Ip, I1, "Q", 1);
     const auto ssp = Ip.get("SS");
     const auto ss1 = I1.get("SS");
     my_assert(abs(ss1 - ssp) == 1);
     return (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
   }

   void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
     auto trSZ = 0.0, trQ = 0.0, trQ2 = 0.0; // Tr[S_z^2], Tr[Q], Tr[Q^2]
     for (const auto &[I, eig]: diag) {
       const auto ss   = I.get("SS");
       const auto q    = I.get("Q");
       const auto sumZ = this->calculate_Z(I, eig, factor);
       trQ  += sumZ * q;
       trQ2 += sumZ * q * q;
       trSZ += sumZ * (ss * ss - 1) / 12.;
     }
     stats.td.set("<Sz^2>", trSZ / stats.Z);
     stats.td.set("<Q>",    trQ  / stats.Z);
     stats.td.set("<Q^2>",  trQ2 / stats.Z);
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

// *** Helper macros for make_matrix() members in matrix.cc
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

/* i - subspace index
   ch - channel (0 or 1)
   number - number of electrons added in channel 'ch' in subspace 'i' */
#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

#undef OFFDIAG_MIX
#define OFFDIAG_MIX(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xiR(step.N(), ch), h, qq, In, opch)

#undef RUNGHOP
#define RUNGHOP(i, j, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * coef.zetaR(step.N() + 1, 0), h, qq)

template<typename SC>
void SymmetryQS<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, 
                                     const Opch<SC> &opch, const Coef<SC> &coef) const {
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
  } else {
    my_assert(P.coeffactor == 1);
    const auto [N, M] = step.NM();

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, M, 0, t_matel(factor0) * coef.xi(N, M), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(N + 1, M), h, qq)

#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"

    if (P.rungs) my_assert_not_reached();
  }
}

template<typename SC>
void SymmetryQS<SC>::show_coefficients(const Step &step, const Coef<SC> &coef) const {
  Symmetry<SC>::show_coefficients(step, coef);
  if (P.rungs)
    for (auto i = 0; i < P.channels; i++)
      std::cout << "[" << i + 1 << "]"
           << " xi_rung(" << step.N() << ")=" << coef.xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" 
           << coef.zetaR(step.N() + 1, i) << std::endl;
}

}

#include "nrg-recalc-QS.hpp"
