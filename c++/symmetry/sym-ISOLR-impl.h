template<typename SC>
class SymmetryISOLRcommon : public SymLR<SC> {
 private:
   outfield Sz2, Q2;

 protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   SymmetryISOLRcommon(const Params &P, Allfields &allfields) : SymLR<SC>(P, Invar(1,1,1)),
     Sz2(P, allfields, "<Sz^2>", 1), Q2(P, allfields, "<Q^2>", 2) {
       initInvar({
         {"II", additive},     // isospin
         {"SS", additive},     // spin
         {"P", multiplicative} // parity
       });
     }

  // Multiplicity of the I=(II,SS,P) subspace = (2I+1)(2S+1) = II SS.
  size_t mult(const Invar &I) const override {
    int mi = I.get("II"); // isospin multiplicity
    int ms = I.get("SS"); // spin multiplicity
    return mi * ms;
  }

  // We always must have S >= 0 and I >= 0.
  bool Invar_allowed(const Invar &I) const override {
    const bool isospin_ok = I.get("II") > 0;
    const bool spin_ok    = I.get("SS") > 0;
    return isospin_ok && spin_ok;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II")) && su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert(abs(ss1 - ssp) == 1);
    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");
    my_assert(abs(ii1 - iip) == 1);
    const double spinfactor = (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
    const double isofactor  = (ii1 == iip + 1 ? ISO(iip) + 1.0 : ISO(iip));
    return spinfactor * isofactor;
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, const Stats<SC> &stats, const double factor) override {
    bucket trSZ, trIZ; // Tr[S_z^2], Tr[I_z^2]
    for (const auto &[I, eig]: diag) {
      const Ispin ii    = I.get("II");
      const Sspin ss    = I.get("SS");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * (ss * ss - 1) / 12.; // isospin multiplicity contained in sumZ
      trIZ += sumZ * (ii * ii - 1) / 12.; // spin multiplicity contained in sumZ
    }
    Sz2 = trSZ / stats.Z;
    Q2  = (4 * trIZ) / stats.Z;
  }
};

template<typename SC>
class SymmetryISOLR : public SymmetryISOLRcommon<SC> {
 private:
   using SymmetryISOLRcommon<SC>::P;
   using SymmetryISOLRcommon<SC>::In;
   using SymmetryISOLRcommon<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISOLR(const Params &P, Allfields &allfields) : SymmetryISOLRcommon<SC>(P, allfields) {}

  void load() override {
    my_assert(P.channels == 2);
#include "isolr/isolr-2ch-In2.dat"
#include "isolr/isolr-2ch-QN.dat"
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

template<typename SC>
class SymmetryISO2LR : public SymmetryISOLRcommon<SC> {
 private:
   using SymmetryISOLRcommon<SC>::P;
   using SymmetryISOLRcommon<SC>::In;
   using SymmetryISOLRcommon<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISO2LR(const Params &P, Allfields &allfields) : SymmetryISOLRcommon<SC>(P, allfields) {}

  void load() override {
    my_assert(P.channels == 2);
#include "iso2lr/iso2lr-2ch-In2.dat"
#include "iso2lr/iso2lr-2ch-QN.dat"
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

// *** Helper macros for make_matrix() members in matrix.cc
#undef OFFIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

template<typename SC>
void SymmetryISOLR<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) {
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");
#include "isolr/isolr-2ch-offdiag.dat"
}

template<typename SC>
void SymmetryISO2LR<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) {
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");
#include "iso2lr/iso2lr-2ch-offdiag.dat"
}

#include "nrg-recalc-ISOLR.h"
#include "nrg-recalc-ISO2LR.h"
