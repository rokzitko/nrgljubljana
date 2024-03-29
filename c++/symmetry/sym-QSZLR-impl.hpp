namespace NRG {

template<typename SC>
class SymmetryQSZLR : public SymFieldLR<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSZLR(const Params &P) : SymFieldLR<SC>(P, std::vector{"<Sz^2>", "<Sz>", "<Q>", "<Q^2>"}, Invar(0,0,1)) {
     initInvar({
        {"Q", additive},      // charge
        {"SSZ", additive},    // spin projection
        {"P", multiplicative} // parity
     });
   }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    const int ssz1  = I1.get("SSZ");
    const int sszp  = Ip.get("SSZ");
    const int sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "qszlr/qszlr-2ch-In2.dat"
#include "qszlr/qszlr-2ch-QN.dat"
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[(S_z)^2], etc.
    for (const auto &[I, eig]: diag) {
      const int ssz  = I.get("SSZ");
      const int q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
      trQ += sumZ * q;
      trQ2 += sumZ * pow(q,2);
    }
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Sz>",   trSZ / stats.Z);
    stats.td.set("<Q>",    trQ / stats.Z);
    stats.td.set("<Q^2>",  trQ2 / stats.Z);
  }

  DECL;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryQSZLR<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
#include "qszlr/qszlr-2ch-offdiag.dat"
#include "qszlr/qszlr-2ch-diag.dat"
}

}

#include "nrg-recalc-QSZLR.hpp"
