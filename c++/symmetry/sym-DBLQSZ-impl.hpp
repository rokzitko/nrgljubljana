namespace NRG {

template<typename SC>
class SymmetryDBLQSZ : public SymField<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryDBLQSZ(const Params &P) : SymField<SC>(P, std::vector{"<Sz^2>", "<Sz>", "<Q1>", "<Q1^2>", "<Q2>", "<Q2^2>"}, Invar(0, 0, 0)) {
     initInvar({
        {"Q1", additive}, // charge 1
        {"Q2", additive}, // charge 2
        {"SSZ", additive}  // spin projection
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
     return u1_equality(I1.get("Q1"),  I2.get("Q1"),  I3.get("Q1"))
         && u1_equality(I1.get("Q2"),  I2.get("Q2"),  I3.get("Q2"))
         && u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ"));
   }

  void load() override {
    switch (P.channels) {
      case 2:
#include "dblqsz/dblqsz-2ch-In2.dat"
#include "dblqsz/dblqsz-2ch-QN.dat"
        break;
      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    auto trSZ = 0.0, trSZ2 = 0.0; // Tr[S_z], Tr[S_z^2]
    auto trQ1 = 0.0, trQ12 = 0.0; // Tr[Q_1], Tr[Q_1^2]
    auto trQ2 = 0.0, trQ22 = 0.0; // Tr[Q_2], Tr[Q_2^2]
    for (const auto &[I, eig]: diag) {
      const int q1  = I.get("Q1");
      const int q2  = I.get("Q2");
      const int ssz  = I.get("SSZ");
      const auto sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
      trQ1 += sumZ * q1;
      trQ12 += sumZ * pow(q1, 2);
      trQ2 += sumZ * q2;
      trQ22 += sumZ * pow(q2, 2);
    }
    stats.td.set("<Sz>",   trSZ  / stats.Z);
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Q1>",   trQ1  / stats.Z);
    stats.td.set("<Q1^2>", trQ12 / stats.Z);
    stats.td.set("<Q2>",   trQ2  / stats.Z);
    stats.td.set("<Q2^2>", trQ22 / stats.Z);
  }

  DECL;
  HAS_DOUBLET;
  HAS_GLOBAL;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) this->diag_function(step, i, ch, number, coef.zeta(step.N() + 1, ch), h, qq)

template<typename SC>
void SymmetryDBLQSZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  switch (P.channels) {
    case 2:
#include "dblqsz/dblqsz-2ch-offdiag.dat"
#include "dblqsz/dblqsz-2ch-diag.dat"
      break;
    default: my_assert_not_reached();
  }
}

}

#include "nrg-recalc-DBLQSZ.hpp"
