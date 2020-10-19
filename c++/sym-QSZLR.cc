template<typename SC>
class SymmetryQSZLR_tmpl : public SymFieldLR_tmpl<SC> {
 private:
   outfield Sz2, Sz, Q, Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   SymmetryQSZLR_tmpl(const Params &P, Allfields &allfields) : SymFieldLR_tmpl<SC>(P),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2), Q(P, allfields, "<Q>", 3), Q2(P, allfields, "<Q^2>", 4) {
       initInvar({
         {"Q", additive},      // charge
         {"SSZ", additive},    // spin projection
         {"P", multiplicative} // parity
       });
       this->InvarSinglet = Invar(0, 0, 1);
     }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    const SZspin ssz1  = I1.get("SSZ");
    const SZspin sszp  = Ip.get("SSZ");
    const SZspin sszop = ssz1 - sszp;
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

  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats &stats, const double factor) override {
    bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[(S_z)^2], etc.
    for (const auto &[I, eig]: diag) {
      const SZspin ssz  = I.get("SSZ");
      const Number q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
      trQ += sumZ * q;
      trQ2 += sumZ * pow(q,2);
    }
    Sz2 = trSZ2 / stats.Z;
    Sz  = trSZ / stats.Z;
    Q   = trQ / stats.Z;
    Q2  = trQ2 / stats.Z;
  }

  DECL;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryQSZLR_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
#include "qszlr/qszlr-2ch-offdiag.dat"
#include "qszlr/qszlr-2ch-diag.dat"
}

#include "nrg-recalc-QSZLR.cc"
