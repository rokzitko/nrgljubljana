template<typename SC>
class SymmetryISOSZLR : public SymFieldLR<SC> {
 private:
   outfield Sz2, Sz, Q2;
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISOSZLR(const Params &P, Allfields &allfields) : SymFieldLR<SC>(P, Invar(1,0,1)),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2), Q2(P, allfields, "<Q^2>", 3) {
       initInvar({
         {"II", additive},     // isospin
         {"SSZ", additive},    // spin projection
         {"P", multiplicative} // parity
       });
     }

  // Multiplicity of the I=(II,SSZ) subspace = (2I+1) = II.
  size_t mult(const Invar &I) const override { return I.get("II"); }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    int ssz1  = I1.get("SSZ");
    int sszp  = Ip.get("SSZ");
    int sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "isoszlr/isoszlr-2ch-In2.dat"
#include "isoszlr/isoszlr-2ch-QN.dat"
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_abs_diff(Ip, I1, "SSZ", 1);
    const int iip = Ip.get("II");
    const int ii1 = I1.get("II");
    const double isofactor = (ii1 == iip + 1 ? ISO(iip) + 1.0 : ISO(iip));
    return isofactor;
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, const Stats<SC> &stats, const double factor) override {
    bucket trSZ, trSZ2, trIZ2; // Tr[S_z], Tr[S_z^2], Tr[I_z^2]
    for (const auto &[I, eig]: diag) {
      const int ii    = I.get("II");
      const int ssz  = I.get("SSZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);        // isospin multiplicity contained in sumZ
      trIZ2 += sumZ * (ii * ii - 1) / 12.; // spin multiplicity contained in sumZ
    }
    Sz  = trSZ / stats.Z;
    Sz2 = trSZ2 / stats.Z;
    Q2  = (4 * trIZ2) / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

template<typename SC>
void SymmetryISOSZLR<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) {
  int ii = I.get("II");
#include "isoszlr/isoszlr-2ch-offdiag.dat"
}

#include "nrg-recalc-ISOSZLR.h"
