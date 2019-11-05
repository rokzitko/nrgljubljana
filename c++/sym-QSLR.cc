class SymmetryQSLR : public SymLR
{
private:
   outfield Sz2, Q, Q2;
   
public:
   SymmetryQSLR() : SymLR() {
      all_syms["QSLR"] = this;
   }
   
   void init() {
      Sz2.set("<Sz^2>", 1);
      Q.set("<Q>", 2);
      Q2.set("<Q^2>", 3);
      InvarStructure InvStruc[] = {
	 {"Q", additive},      // charge 
	 {"SS", additive},     // spin 
	 {"P", multiplicative} // parity 
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1, 1);
   }

  // Multiplicity of the I=(Q,SS,P) subspace = 2S+1 = SS.
  int mult(const Invar &I) {
    return I.get("SS"); // spin multiplicity
  }

  bool Invar_allowed(const Invar &I) {
    return I.get("SS") > 0;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
      z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() {
    my_assert(channels == 2);
#include "qslr/qslr-2ch-In2.dat"
#include "qslr/qslr-2ch-QN.dat"
  }

  double dynamicsusceptibility_factor(const Invar &Ip,
                                      const Invar &I1) {
    check_diff(Ip, I1, "Q", 0);
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert((abs(ss1-ssp) == 2 || ss1 == ssp));

    return switch3(ss1,
                   ssp+2, 1.+(ssp-1)/3.,
                   ssp,   ssp/3.,
                   ssp-2, (-2.+ssp)/3.);
  }

  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    check_diff(Ip, I1, "Q", 1);
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    return ( ss1 == ssp+1 ? S(ssp) + 1.0 : S(ssp) );
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ, trQ, trQ2; // Tr[S_z^2], Tr[Q], Tr[Q^2]

    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Sspin ss = I.get("SS");
      const Number q = I.get("Q");
      const double sumZ = calculate_Z(is, factor);

      trQ  += sumZ * q;
      trQ2 += sumZ * q * q;
      trSZ += sumZ * (ss*ss-1)/12.;
    }

    Sz2 = trSZ/STAT::Z;
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymQSLR = new SymmetryQSLR;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)
#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
                                 zeta(STAT::N+1, ch), h, qq)

void SymmetryQSLR::makematrix(Matrix &h, const Rmaxvals &qq,
                              const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
#include "qslr/qslr-2ch-diag.dat"
#include "qslr/qslr-2ch-offdiag.dat"
}

#include "nrg-recalc-QSLR.cc"
