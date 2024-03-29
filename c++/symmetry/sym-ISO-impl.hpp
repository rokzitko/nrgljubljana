namespace NRG {

template<typename SC>
class SymmetryISOcommon : public Symmetry<SC> {
  protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

  public:
  SymmetryISOcommon(const Params &P) : Symmetry<SC>(P, std::vector{"<Sz^2>", "<Q^2>"}, Invar(1, 1), Invar(2, 2)) {
    initInvar({
       {"II", additive}, // isospin
       {"SS", additive}  // spin
    });
  }

  // Multiplicity of the I=(II,SS) subspace = (2I+1)(2S+1) = II SS.
  size_t mult(const Invar &I) const override {
    int mi = I.get("II"); // isospin multiplicity
    int ms = I.get("SS"); // spin multiplicity
    return mi * ms;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) && su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"));
  }

  // We always must have S >= 0 and I >= 0.
  bool Invar_allowed(const Invar &I) const override {
    const bool isospin_ok = I.get("II") > 0;
    const bool spin_ok    = I.get("SS") > 0;
    return isospin_ok && spin_ok;
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    const int ssp = Ip.get("SS");
    const int ss1 = I1.get("SS");
    my_assert(abs(ss1 - ssp) == 1);

    const int iip = Ip.get("II");
    const int ii1 = I1.get("II");
    my_assert(abs(ii1 - iip) == 1);

    const double spinfactor = (ss1 == ssp + 1 ? S(ssp) + 1.0 : S(ssp));
    const double isofactor  = (ii1 == iip + 1 ? ISO(iip) + 1.0 : ISO(iip));

    return spinfactor * isofactor;
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trIZ; // Tr[S_z^2], Tr[I_z^2]
    for (const auto &[I, eig]: diag) {
      const int ii    = I.get("II");
      const int ss    = I.get("SS");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * (ss * ss - 1) / 12.; // isospin multiplicity contained in sumZ
      trIZ += sumZ * (ii * ii - 1) / 12.; // spin multiplicity contained in sumZ
    }
    stats.td.set("<Sz^2>", trSZ / stats.Z);
    stats.td.set("<Q^2>",  (4 * trIZ) / stats.Z);
  }
};

template<typename SC>
class SymmetryISO : public SymmetryISOcommon<SC> {
  private:
   using SymmetryISOcommon<SC>::P;
   using SymmetryISOcommon<SC>::In;
   using SymmetryISOcommon<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISO(const Params &P) : SymmetryISOcommon<SC>(P) {}

  void load() override {
    switch (P.channels) {
      case 1:
#include "iso/iso-1ch-In2.dat"
#include "iso/iso-1ch-QN.dat"
        break;
      case 2:
#include "iso/iso-2ch-In2.dat"
#include "iso/iso-2ch-QN.dat"
        break;
      case 3:
#include "iso/iso-3ch-In2.dat"
#include "iso/iso-3ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

template<typename SC>
class SymmetryISO2 : public SymmetryISOcommon<SC> {
  private:
   using SymmetryISOcommon<SC>::P;
   using SymmetryISOcommon<SC>::In;
   using SymmetryISOcommon<SC>::QN;

  public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISO2(const Params &P) : SymmetryISOcommon<SC>(P) {}

  void load() override {
    switch (P.channels) {
      case 1:
#include "iso2/iso2-1ch-In2.dat"
#include "iso2/iso2-1ch-QN.dat"
        break;
      case 2:
#include "iso2/iso2-2ch-In2.dat"
#include "iso2/iso2-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryISO<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int ss = I.get("SS");
  int ii = I.get("II");
  // nn is the index of the last site in the chain, while nn+1 is the
  // index of the site that is being added to the chain in this
  // iteration.  This is consistent with the definition in
  // isospin-1ch-automatic.nb.
  int NN = step.getnn();
  switch (P.channels) {
    case 1:
#include "iso/iso-1ch-offdiag.dat"
      break;
    case 2:
#include "iso/iso-2ch-offdiag.dat"
      break;
    case 3:
#include "iso/iso-3ch-offdiag.dat"
      break;
    default: my_assert_not_reached();
  }
}

template<typename SC>
void SymmetryISO2<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  int ss = I.get("SS");
  int ii = I.get("II");
  int NN   = step.getnn();
  switch (P.channels) {
    case 1:
#include "iso2/iso2-1ch-offdiag.dat"
      break;
    case 2:
#include "iso2/iso2-2ch-offdiag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-ISO.hpp"
#include "nrg-recalc-ISO2.hpp"
