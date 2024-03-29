namespace NRG {

template<typename SC>
class SymmetryQSZTZ : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSZTZ(const Params &P) : Symmetry<SC>(P, std::vector{"<Sz>", "<Sz^2>", "<Tz>", "<Tz^2>", "<Q>", "<Q^2>"}, Invar(0,0,0), Invar(1,2,1)) {
     initInvar({
        {"Q", additive},  // charge
        {"SZ", additive}, // spin
        {"TZ", additive}  // angular momentum
     });
   }
  size_t mult(const Invar &I) const override { return 1; }
  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) && u1_equality(I1.get("SZ"), I2.get("SZ"), I3.get("SZ"))
       && u1_equality(I1.get("TZ"), I2.get("TZ"), I3.get("TZ"));
  }
  bool Invar_allowed(const Invar &I) const override { return true; }
  void load() override {
    my_assert(!P.substeps);
    my_assert(P.channels == 3);
#include "qsztz/qsztz-In2.dat"
#include "qsztz/qsztz-QN.dat"
  } // load
  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trSZ2, trTZ, trTZ2, trQ, trQ2;
    for (const auto &[I, eig]: diag) {
      const int q    = I.get("Q");
      const int ssz   = I.get("SZ");
      const int tz  = I.get("TZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
      trSZ += sumZ * (ssz - 1) / 2.;
      trSZ2 += sumZ * pow((ssz - 1) / 2., 2);
      trTZ += sumZ * tz;
      trTZ2 += sumZ * tz * tz;
    }
    stats.td.set("<Sz>",   trSZ / stats.Z);
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Tz>",   trTZ / stats.Z);
    stats.td.set("<Tz^2>", trTZ2 / stats.Z);
    stats.td.set("<Q>",    trQ / stats.Z);
    stats.td.set("<Q^2>",  trQ2 / stats.Z);
  }
  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

// We take the coefficients of the first channel (indexed as 0), because all three set are exactly the same due to
// orbital symmetry. [Can be generalized for symmetry-broken cases if the need arises.]
#undef OFFDIAG
#define OFFDIAG(i, j, factor0) offdiag_function(step, i, j, 0, 0, t_matel(factor0) * coef.xi(step.N(), 0), h, qq, In, opch)

#undef DIAG
#define DIAG(i, number) this->diag_function(step, i, number, coef.zeta(step.N() + 1, 0), h, qq)

template<typename SC>
void SymmetryQSZTZ<SC>::make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, const Opch<SC> &opch, const Coef<SC> &coef) const {
  my_assert(!P.substeps);
  my_assert(P.channels == 3);
#include "qsztz/qsztz-offdiag.dat"
#include "qsztz/qsztz-diag.dat"
}

}

#include "nrg-recalc-QSZTZ.hpp"
