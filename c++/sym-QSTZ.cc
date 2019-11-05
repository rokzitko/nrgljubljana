class SymmetryQSTZ : public Symmetry 
{
private:
   outfield Sz2, Tz2, Q, Q2;

public:
   SymmetryQSTZ() : Symmetry() {
      all_syms["QSTZ"] = this;
   }
  
   void init() {
      Sz2.set("<Sz^2>", 1);
      Tz2.set("<Tz^2>", 2);
      Q.set("<Q>", 3);
      Q2.set("<Q^2>", 4);
      InvarStructure InvStruc[] = {
	 {"Q", additive},  // charge
	 {"SS", additive}, // spin
	 {"TZ", additive}   // angular momentum
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1, 0);
      Invar_f = Invar(1, 2, 1); // (1,2,1) is correct. see triangle_inequality() below!
  }

  // Multiplicity of the (Q,SS,TZ) subspace is (2S+1 = SS).
  int mult(const Invar &I) { return I.get("SS"); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
      u1_equality(I1.get("TZ"), I2.get("TZ"), I3.get("TZ"));
  }

  bool Invar_allowed(const Invar &I) {
    return I.get("SS") > 0;
  }

  void load() {
    my_assert(!substeps);
    my_assert(channels == 3);
#include "qstz/qstz-In2.dat"
#include "qstz/qstz-QN.dat"
  } // load

  // Same as for SYMTYPE=QS, because spin operators are angular momentum singlets.
  double dynamicsusceptibility_factor(const Invar &Ip,
                                      const Invar &I1) {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "TZ", 0);

    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert((abs(ss1-ssp) == 2 || ss1 == ssp));

    return switch3(ss1,
                   ssp+2, 1.+(ssp-1)/3.,
                   ssp,   ssp/3.,
                   ssp-2, (-2.+ssp)/3.);
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    check_diff(Ip, I1, "Q", 1);
    
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert( abs(ss1-ssp) == 1 );
     
    return ( ss1 == ssp+1 ? S(ssp) + 1.0 : S(ssp) );
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ2, trTZ2, trQ, trQ2; // Tr[S_z^2], Tr[T_z^2], Tr[Q], Tr[Q^2]

    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Number q = I.get("Q");
      const Sspin ss = I.get("SS");
      const Tangmom tz = I.get("TZ");
      const double sumZ = calculate_Z(is, factor);

      trQ  += sumZ * q;
      trQ2 += sumZ * q * q;
      trSZ2 += sumZ * (ss*ss-1)/12.; // [(2S+1)(2S+1)-1]/12=S(S+1)/3
      trTZ2 += sumZ * tz * tz;
    }

    Sz2 = trSZ2/STAT::Z;
    Tz2 = trTZ2/STAT::Z;
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymQSTZ = new SymmetryQSTZ;

// We take the coefficients of the first channel (indexed as 0),
// because all three set are exactly the same due to orbital
// symmetry. [XXXX True for QST, for QSTZ we could relax this approx.]
#undef OFFDIAG
#define OFFDIAG(i, j, factor0) offdiag_function(i, j, 0, 0, \
						t_matel(factor0) * xi(STAT::N, 0), h, qq, In)

#undef DIAG
#define DIAG(i, number) diag_function(i, 0, number, \
                              zeta(STAT::N+1, 0), h, qq)

void SymmetryQSTZ::makematrix(Matrix &h, const Rmaxvals &qq,
                             const Invar &I, const InvarVec &In)
{
   Sspin ss = I.get("SS");
   
   my_assert(!substeps);
   my_assert(channels == 3);

#include "qstz/qstz-offdiag.dat"
#include "qstz/qstz-diag.dat"
}

#include "nrg-recalc-QSTZ.cc"
