template<typename SC>
class SymmetryQSZ_tmpl : public SymField_tmpl<SC> {
 private:
   outfield Sz2, Sz, Q, Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   SymmetryQSZ_tmpl(const Params &P, Allfields &allfields) : SymField_tmpl<SC>(P),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2), Q(P, allfields, "<Q>", 3), Q2(P, allfields, "<Q^2>", 4) {
       initInvar({
         {"Q", additive},  // charge
         {"SSZ", additive} // spin projection
       });
       this->InvarSinglet = Invar(0, 0);
       this->Invar_f      = Invar(1, 2);
     }

   bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
     // The spin projection of the operator is defined by the difference in Sz of both the invariant subspaces.
     SZspin ssz1  = I1.get("SSZ");
     SZspin sszp  = Ip.get("SSZ");
     SZspin sszop = ssz1 - sszp;
     return sszop == SPIN;
   }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"));
   }

   void load() override {
     if (!P.substeps) {
       switch (P.channels) {
       case 1:
#include "qsz/qsz-1ch-In2.dat"
#include "qsz/qsz-1ch-QN.dat"
         break;
       case 2:
#include "qsz/qsz-2ch-In2.dat"
#include "qsz/qsz-2ch-QN.dat"
         break;
       case 3:
#include "qsz/qsz-3ch-In2.dat"
#include "qsz/qsz-3ch-QN.dat"
         break;
       default: my_assert_not_reached();
       } // switch
     } else {
#include "qsz/qsz-1ch-In2.dat"
#include "qsz/qsz-1ch-QN.dat"
     } // if
   }

   void make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                              const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef);
   void make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                                 const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef);

   void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats_tmpl<SC> &stats, const double factor) override {
     bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[(S_z)^2], etc.
     for (const auto &[I, eig]: diag) {
       const SZspin ssz  = I.get("SSZ");
       const Number q    = I.get("Q");
       const double sumZ = this->calculate_Z(I, eig, factor);
       trSZ += sumZ * SZ(ssz);
       trSZ2 += sumZ * sqr(SZ(ssz));
       trQ += sumZ * q;
       trQ2 += sumZ * sqr(q);
     }
     Sz2 = trSZ2 / stats.Z;
     Sz  = trSZ / stats.Z;
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

// Note ch indexes the <||f||> matrix which is used to construct the Hamiltonian matrix in the new step, i.e., the
// f_{N} from the f^\dag_{N_1} f_{N} hopping term.
#undef OFFDIAG_MIX
#define OFFDIAG_MIX(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.xiR(step.N(), ch), h, qq, In, opch)

#undef RUNGHOP
#define RUNGHOP(i, j, factor) this->diag_offdiag_function(step, i, j, 0, t_matel(factor) * coef.zetaR(step.N() + 1, 0), h, qq)

// "non-polarized" here means that the coefficients xi do not depend on spin. Note, however, that there is support
// for a global magnetic field, cf. P.globalB.
template<typename SC>
void SymmetryQSZ_tmpl<SC>::make_matrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                                                    const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  if (!P.substeps) {
    switch (P.channels) {
      case 1:
#include "qsz/qsz-1ch-offdiag.dat"
#include "qsz/qsz-1ch-diag.dat"
        break;
      case 2:
#include "qsz/qsz-2ch-offdiag.dat"
#include "qsz/qsz-2ch-diag.dat"
        if (P.rungs) {
#include "qsz/qsz-2ch-offdiag-mix.dat"
#include "qsz/qsz-2ch-runghop.dat"
        }
        break;
      case 3:
#include "qsz/qsz-3ch-offdiag.dat"
#include "qsz/qsz-3ch-diag.dat"
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

#include "qsz/qsz-1ch-offdiag.dat"
#include "qsz/qsz-1ch-diag.dat"

    if (P.rungs) my_assert_not_reached();
  }
}

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DOWN(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiDOWN(step.N(), ch), h, qq, In, opch)

#define DIAG_UP(i, j, ch, number) this->diag_function_half(step, i, ch, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) this->diag_function_half(step, i, ch, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryQSZ_tmpl<SC>::make_matrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                                                 const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  my_assert(!P.substeps); // not implemented!
  switch (P.channels) {
    case 1:
#include "qsz/qsz-1ch-offdiag-UP.dat"
#include "qsz/qsz-1ch-offdiag-DOWN.dat"
#include "qsz/qsz-1ch-diag-UP.dat"
#include "qsz/qsz-1ch-diag-DOWN.dat"
      break;
    case 2:
#include "qsz/qsz-2ch-offdiag-UP.dat"
#include "qsz/qsz-2ch-offdiag-DOWN.dat"
#include "qsz/qsz-2ch-diag-UP.dat"
#include "qsz/qsz-2ch-diag-DOWN.dat"
      if (P.rungs) {
        //#include "qsz/qsz-2ch-offdiag-mix-UP.dat"
        //#include "qsz/qsz-2ch-offdiag-mix-DOWN.dat"
        //#include "qsz/qsz-2ch-runghop-UP.dat"
        //#include "qsz/qsz-2ch-runghop-DOWN.dat"
        my_assert_not_reached();
      }
      break;
    case 3:
#include "qsz/qsz-3ch-offdiag-UP.dat"
#include "qsz/qsz-3ch-offdiag-DOWN.dat"
#include "qsz/qsz-3ch-diag-UP.dat"
#include "qsz/qsz-3ch-diag-DOWN.dat"
      break;
    default: my_assert_not_reached();
  }
}

template<typename SC>
void SymmetryQSZ_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In,
                                       const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  if (P.polarized) 
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  else
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
}

template<typename SC>
void SymmetryQSZ_tmpl<SC>::show_coefficients(const Step &step, const Coef_tmpl<SC> &coef) {
  Symmetry_tmpl<SC>::show_coefficients(step, coef);
  if (P.rungs) 
    for (unsigned int i = 0; i < P.channels; i++) {
      std::cout << "[" << i + 1 << "]"
        << " xi_rung(" << step.N() << ")=" << coef.xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" << coef.zetaR(step.N() + 1, i) << std::endl;
    }
}

#include "nrg-recalc-QSZ.cc"
