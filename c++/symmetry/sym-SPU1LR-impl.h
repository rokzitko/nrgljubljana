template<typename SC>
class SymmetrySPU1LR : public SymFieldLR<SC> {
 private:
   outfield Sz2, Sz;
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPU1LR(const Params &P, Allfields &allfields) : SymFieldLR<SC>(P),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2) {
       initInvar({
         {"SSZ", additive},    // spin projection
         {"P", multiplicative} // parity
       });
       this->InvarSinglet = Invar(0, 1);
     }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1  = I1.get("SSZ");
    SZspin sszp  = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    switch (P.channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-In2.dat"
#include "spu1lr/spu1lr-1ch-QN.dat"
        break;
      case 2:
#include "spu1lr/spu1lr-2ch-In2.dat"
#include "spu1lr/spu1lr-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, const Stats<SC> &stats, const double factor) override {
    bucket trSZ, trSZ2; // Tr[S_z], Tr[S_z^2]
    for (const auto &[I, eig]: diag) {
      const SZspin ssz  = I.get("SSZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
    }
    Sz2 = trSZ2 / stats.Z;
    Sz  = trSZ / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_GLOBAL;
};

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, ch, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetrySPU1LR<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) {
  switch (P.channels) {
    case 1:
#include "spu1lr/spu1lr-1ch-offdiag.dat"
#include "spu1lr/spu1lr-1ch-anomalous.dat"
#include "spu1lr/spu1lr-1ch-diag.dat"
#include "spu1lr/spu1lr-1ch-isospinx.dat"
      break;
    case 2:
#include "spu1lr/spu1lr-2ch-diag.dat"
#include "spu1lr/spu1lr-2ch-offdiag.dat"
#include "spu1lr/spu1lr-2ch-anomalous.dat"
#include "spu1lr/spu1lr-2ch-isospinx.dat"
      break;
    default: my_assert_not_reached();
  }
}

#include "nrg-recalc-SPU1LR.h"
