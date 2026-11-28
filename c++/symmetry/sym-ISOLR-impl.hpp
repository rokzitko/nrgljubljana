#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "table_terms.hpp"

namespace NRG {

using std::sqrt;

template<typename SC>
class SymmetryISOLRcommon : public SymLR<SC> {
 protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
 public:
   SymmetryISOLRcommon(const Params &P) : SymLR<SC>(P, {td::su2_z_square("<Sz^2>", "SS"), td::su2_z_square("<Q^2>", "II", 4)}, Invar(1,1,1)) {
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
    return su2_doublet_factor(Ip.get("SS"), I1.get("SS")) * su2_doublet_factor(Ip.get("II"), I1.get("II"));
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
   SymmetryISOLR(const Params &P) : SymmetryISOLRcommon<SC>(P) {}

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
   SymmetryISO2LR(const Params &P) : SymmetryISOLRcommon<SC>(P) {}

  void load() override {
    my_assert(P.channels == 2);
#include "iso2lr/iso2lr-2ch-In2.dat"
#include "iso2lr/iso2lr-2ch-QN.dat"
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};


template<typename SC>
void SymmetryISOLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ss = I.get("SS");
  int ii = I.get("II");
#include "isolr/isolr-2ch-offdiag.dat"
}

template<typename SC>
void SymmetryISO2LR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  const ShellTerms<SC> terms(*this, P, step, h, qq, In, opch, coef);
  int ss = I.get("SS");
  int ii = I.get("II");
#include "iso2lr/iso2lr-2ch-offdiag.dat"
}

}

#include "nrg-recalc-ISOLR.hpp"
#include "nrg-recalc-ISO2LR.hpp"
