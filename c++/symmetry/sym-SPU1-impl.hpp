namespace NRG {

template<typename SC>
class SymmetrySPU1 : public  SymField<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetrySPU1(const Params &P) : SymField<SC>(P, std::vector{"<Sz^2>", "<Sz>"}, Invar(0)) {
     initInvar({
        {"SSZ", additive} // spin projection
     });
   }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference in Sz of both the invariant subspaces.
    int ssz1  = I1.get("SSZ");
    int sszp  = Ip.get("SSZ");
    int sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override { return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")); }

  void load() override {
    if (!P.substeps) {
      switch (P.channels) {
        case 1:
#include "spu1/spu1-1ch-In2.dat"
#include "spu1/spu1-1ch-QN.dat"
          break;
        case 2:
#include "spu1/spu1-2ch-In2.dat"
#include "spu1/spu1-2ch-QN.dat"
          break;
        default: my_assert_not_reached();
      } // switch
    } else {
#include "spu1/spu1-1ch-In2.dat"
#include "spu1/spu1-1ch-QN.dat"
    }
  }

  void make_matrix_polarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const;
  void make_matrix_nonpolarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const;

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
  HAS_SUBSTEPS;
};

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetrySPU1<SC>::make_matrix_nonpolarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  if (!P.substeps) {
    switch (P.channels) {
      case 1:
#include "spu1/spu1-1ch-offdiag.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag.dat"
#include "spu1/spu1-1ch-isospinx.dat"
        break;
      case 2:
#include "spu1/spu1-2ch-diag.dat"
#include "spu1/spu1-2ch-offdiag.dat"
#include "spu1/spu1-2ch-anomalous.dat"
#include "spu1/spu1-2ch-isospinx.dat"
        break;
      default: my_assert_not_reached();
    }
  } else {
    my_assert(P.coeffactor == 1);
    const auto [Ntrue, M] = step.NM();

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(Ntrue + 1, M), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, M, 0, t_matel(factor) * coef.kappa(Ntrue, M), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, M, 0, t_matel(factor0) * coef.xi(Ntrue, M), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(Ntrue + 1, M), h, qq)

#include "spu1/spu1-1ch-offdiag.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag.dat"
#include "spu1/spu1-1ch-isospinx.dat"
  }
}

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) this->diag_offdiag_function(step, i, j, t_matel(factor) * 2.0 * coef.delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(step, i, j, ch, 0, t_matel(factor) * coef.kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG_UP
#undef OFFDIAG_DOWN
#undef DIAG_UP
#undef DIAG_DOWN

#define OFFDIAG_UP(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiUP(step.N(), ch), h, qq, In, opch)

#define OFFDIAG_DOWN(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xiDOWN(step.N(), ch), h, qq, In, opch)

#define DIAG_UP(i, j, ch, number) this->diag_function_half(step, i, number, coef.zetaUP(step.N() + 1, ch), h, qq)

#define DIAG_DOWN(i, j, ch, number) this->diag_function_half(step, i, number, coef.zetaDOWN(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetrySPU1<SC>::make_matrix_polarized(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  my_assert(!P.substeps);
  switch (P.channels) {
  case 1:
#include "spu1/spu1-1ch-offdiag-UP.dat"
#include "spu1/spu1-1ch-offdiag-DOWN.dat"
#include "spu1/spu1-1ch-anomalous.dat"
#include "spu1/spu1-1ch-diag-UP.dat"
#include "spu1/spu1-1ch-diag-DOWN.dat"
#include "spu1/spu1-1ch-isospinx.dat"
    break;
  case 2:
#include "spu1/spu1-2ch-diag-UP.dat"
#include "spu1/spu1-2ch-diag-DOWN.dat"
#include "spu1/spu1-2ch-offdiag-UP.dat"
#include "spu1/spu1-2ch-offdiag-DOWN.dat"
#include "spu1/spu1-2ch-anomalous.dat"
#include "spu1/spu1-2ch-isospinx.dat"
    break;
  default: my_assert_not_reached();
  }
}

template<typename SC>
void SymmetrySPU1<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  if (P.polarized) {
    make_matrix_polarized(h, step, qq, I, In, opch, coef);
  } else {
    make_matrix_nonpolarized(h, step, qq, I, In, opch, coef);
  }
}

}

#include "nrg-recalc-SPU1.hpp"
