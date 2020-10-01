class SymmetryQSZ : public SymField {
  private:
  outfield Sz2, Sz, Q, Q2;

  public:
  SymmetryQSZ() : SymField() { all_syms["QSZ"] = this; }

  void init() override {
    Sz2.set("<Sz^2>", 1);
    Sz.set("<Sz>", 2);
    Q.set("<Q>", 3);
    Q2.set("<Q^2>", 4);
    InvarStructure InvStruc[] = {
       {"Q", additive},  // charge
       {"SSZ", additive} // spin projection
    };
    initInvar(InvStruc, ARRAYLENGTH(InvStruc));
    InvarSinglet = Invar(0, 0);
    Invar_f      = Invar(1, 2);
  }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1  = I1.get("SSZ");
    SZspin sszp  = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"));
  }

  void load() override {
    if (!substeps) {
      switch (channels) {
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

  void makematrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);
  void makematrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch);

  void calculate_TD(const DiagInfo &diag, double factor) override {
    bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[(S_z)^2], etc.

    for (const auto &[I, eig]: diag) {
      const SZspin ssz  = I.get("SSZ");
      const Number q    = I.get("Q");
      const double sumZ = calculate_Z(I, eig, factor);

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

  void show_coefficients(const Step &, const Params &) override;
};

Symmetry *SymQSZ = new SymmetryQSZ;

// *** Helper macros for makematrix() members in matrix.cc
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xi(step.N(), ch), h, qq, In, opch)

/* i - subspace index
   ch - channel (0 or 1)
   number - number of electrons added in channel 'ch' in subspace 'i' */

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(step.N() + 1, ch), h, qq)

#undef SPINZ
#define SPINZ(i, j, ch, factor) spinz_function(i, j, ch, t_matel(factor), h, qq)

// Note ch indexes the <||f||> matrix which is used to construct the Hamiltonian
// matrix in the new step, i.e., the f_{N} from the f^\dag_{N_1} f_{N} hopping
// term.
#undef OFFDIAG_MIX
#define OFFDIAG_MIX(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * xiR(step.N(), ch), h, qq, In, opch)

#undef RUNGHOP
#define RUNGHOP(i, j, factor) diag_offdiag_function(i, j, 0, t_matel(factor) * zetaR(step.N() + 1, 0), h, qq)

// "non-polarized" here means that the coefficients xi do not depend on
// spin. Note, however, that there is support for a global magnetic field,
// cf. P.globalB.
void SymmetryQSZ::makematrix_nonpolarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (!substeps) {
    switch (channels) {
      case 1:
#include "qsz/qsz-1ch-offdiag.dat"
#include "qsz/qsz-1ch-diag.dat"
#include "qsz/qsz-1ch-spinz.dat" // for P.globalB
        break;

      case 2:
#include "qsz/qsz-2ch-offdiag.dat"
#include "qsz/qsz-2ch-diag.dat"
#include "qsz/qsz-2ch-spinz.dat" // for P.globalB
        if (P.rungs) {
#include "qsz/qsz-2ch-offdiag-mix.dat"
#include "qsz/qsz-2ch-runghop.dat"
        }
        break;

      case 3:
#include "qsz/qsz-3ch-offdiag.dat"
#include "qsz/qsz-3ch-diag.dat"
#include "qsz/qsz-3ch-spinz.dat" // for P.globalB
        break;

      default: my_assert_not_reached();
    }
  } else { // substeps
    my_assert(P.coeffactor == 1);
    const auto [N, M] = step.NM();

    // Overrides. See sym-QS.cc for explanations!
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, M, 0, t_matel(factor0) * xi(N, M) / step.scale_fix(), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, M, number, zeta(N + 1, M), h, qq)

#undef SPINZ
#define SPINZ(i, j, ch, factor) spinz_function(i, j, M, t_matel(factor), h, qq)

#include "qsz/qsz-1ch-offdiag.dat"
#include "qsz/qsz-1ch-diag.dat"
#include "qsz/qsz-1ch-spinz.dat" // for P.globalB

    if (P.rungs) my_error("Not implemented.");
  }
}

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DOWN(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xiDOWN(step.N(), ch), h, qq, In, opch)

#define DIAG_UP(i, j, ch, number) diag_function_half(i, ch, number, zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) diag_function_half(i, ch, number, zetaDOWN(step.N() + 1, ch), h, qq)

#undef SPINZ
#define SPINZ(i, j, ch, factor) spinz_function(i, j, ch, t_matel(factor), h, qq)

void SymmetryQSZ::makematrix_polarized(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  my_assert(!substeps); // not implemented!

  switch (channels) {
    case 1:
#include "qsz/qsz-1ch-offdiag-UP.dat"
#include "qsz/qsz-1ch-offdiag-DOWN.dat"
#include "qsz/qsz-1ch-diag-UP.dat"
#include "qsz/qsz-1ch-diag-DOWN.dat"
#include "qsz/qsz-1ch-spinz.dat"
      break;

    case 2:
#include "qsz/qsz-2ch-offdiag-UP.dat"
#include "qsz/qsz-2ch-offdiag-DOWN.dat"
#include "qsz/qsz-2ch-diag-UP.dat"
#include "qsz/qsz-2ch-diag-DOWN.dat"
#include "qsz/qsz-2ch-spinz.dat"
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
#include "qsz/qsz-3ch-spinz.dat"
      break;

    default: my_assert_not_reached();
  }
}

void SymmetryQSZ::makematrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  if (P.polarized) 
    makematrix_polarized(h, step, qq, I, In, opch);
  else
    makematrix_nonpolarized(h, step, qq, I, In, opch);
}

void SymmetryQSZ::show_coefficients(const Step &step, const Params &P) {
  Symmetry::show_coefficients(step, P);
  if (P.rungs) 
    for (unsigned int i = 0; i < P.channels; i++) 
      cout << "[" << i + 1 << "]"
           << " xi_rung(" << step.N() << ")=" << xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" << zetaR(step.N() + 1, i) << endl;
}

#include "nrg-recalc-QSZ.cc"
