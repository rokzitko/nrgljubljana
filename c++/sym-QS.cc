template<typename SC>
class SymmetryQS_tmpl : public Symmetry_tmpl<SC> {
 private:
   outfield Sz2, Q, Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename Symmetry_tmpl<SC>::Matrix;
   explicit SymmetryQS_tmpl(const Params &P, Allfields &allfields) : Symmetry_tmpl<SC>(P),
     Sz2(P, allfields, "<Sz^2>", 1), Q(P, allfields, "<Q>", 2), Q2(P, allfields, "<Q^2>", 3) {
       initInvar({
         {"Q", additive}, // charge
         {"SS", additive} // spin
       });
       this->InvarSinglet = Invar(0, 1);
       this->Invar_f      = Invar(1, 2);
     }

   // Multiplicity of the (Q,SS) subspace is 2S+1 = SS.
   size_t mult(const Invar &I) const override { return I.get("SS"); }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"));
   }

   bool Invar_allowed(const Invar &I) const override {
     const bool spin_ok = I.get("SS") > 0;
     return spin_ok;
   }

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
     const Sspin ssp = Ip.get("SS");
     const Sspin ss1 = I1.get("SS");
     my_assert((abs(ss1 - ssp) == 2 || ss1 == ssp));
    return switch3(ss1, ssp + 2, 1. + (ssp - 1) / 3., ssp, ssp / 3., ssp - 2, (-2. + ssp) / 3.);
   }

   double specdens_factor(const Invar &Ip, const Invar &I1) const override {
     check_diff(Ip, I1, "Q", 1);
     const Sspin ssp = Ip.get("SS");
     const Sspin ss1 = I1.get("SS");
     my_assert(abs(ss1 - ssp) == 1);
     return (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
   }

   void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats_tmpl<SC> &stats, const double factor) override {
     bucket trSZ, trQ, trQ2; // Tr[S_z^2], Tr[Q], Tr[Q^2]
     for (const auto &[I, eig]: diag) {
       const Sspin ss    = I.get("SS");
       const Number q    = I.get("Q");
       const double sumZ = this->calculate_Z(I, eig, factor);
       trQ += sumZ * q;
       trQ2 += sumZ * q * q;
       trSZ += sumZ * (ss * ss - 1) / 12.;
     }
     Sz2 = trSZ / stats.Z;
     Q   = trQ / stats.Z;
     Q2  = trQ2 / stats.Z;
   }

   DECL;
   HAS_DOUBLET;
   HAS_TRIPLET;
   HAS_GLOBAL;
   HAS_SUBSTEPS;
   void show_coefficients(const Step &, const Coef_tmpl<SC> &) override;
};

// *** Helper macros for make_matrix() members in matrix.cc
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

/* i - subspace index
   ch - channel (0 or 1)
   number - number of electrons added in channel 'ch' in subspace 'i' */
#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

#undef OFFDIAG_MIX
#define OFFDIAG_MIX(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xiR(step.N(), ch), h, qq, In, opch)

#undef RUNGHOP
#define RUNGHOP(i, j, factor) this->diag_offdiag_function(step, i, j, 0, t_matel(factor) * coef.zetaR(step.N() + 1, 0), h, qq)

template<typename SC>
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives; must appear after template
void SymmetryQS_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                                     const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  Sspin ss = I.get("SS");

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
#define DIAG(i, ch, number) this->diag_function(step, i, M, number, coef.zeta(N + 1, M), h, qq)

#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"

    if (P.rungs) my_assert_not_reached();
  }
}

template<typename SC>
void SymmetryQS_tmpl<SC>::show_coefficients(const Step &step, const Coef_tmpl<SC> &coef) {
  Symmetry_tmpl<SC>::show_coefficients(step, coef);
  if (P.rungs)
    for (unsigned int i = 0; i < P.channels; i++)
      std::cout << "[" << i + 1 << "]"
           << " xi_rung(" << step.N() << ")=" << coef.xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" 
           << coef.zetaR(step.N() + 1, i) << std::endl;
}

#include "nrg-recalc-QS.cc"
