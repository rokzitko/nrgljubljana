namespace NRG {

template<typename SC>
class SymmetrySPU1LR : public SymFieldLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPU1LR(const Params &P) : SymFieldLR<SC>(P, std::vector{"<Sz^2>", "<Sz>"}, Invar(0,1)) {
     initInvar({
        {"SSZ", additive},    // spin projection
        {"P", multiplicative} // parity
     });
   }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    int ssz1  = I1.get("SSZ");
    int sszp  = Ip.get("SSZ");
    int sszop = ssz1 - sszp;
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

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trSZ2; // Tr[S_z], Tr[S_z^2]
    for (const auto &[I, eig]: diag) {
      const int ssz  = I.get("SSZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
    }
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Sz>",   trSZ / stats.Z);
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
void SymmetrySPU1LR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
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

}

#include "nrg-recalc-SPU1LR.hpp"
