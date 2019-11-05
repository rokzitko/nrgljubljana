class SymmetryQST : public Symmetry 
{
private:
   outfield Sz2, Tz2, Q, Q2;

public:
   SymmetryQST() : Symmetry() {
      all_syms["QST"] = this;
   }
  
   void init() {
      Sz2.set("<Sz^2>", 1);
      Tz2.set("<Tz^2>", 2);
      Q.set("<Q>", 3);
      Q2.set("<Q^2>", 4);
      InvarStructure InvStruc[] = {
	 {"Q", additive},  // charge
	 {"SS", additive}, // spin
	 {"T", additive}   // angular momentum
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1, 0);
      Invar_f = Invar(1, 2, 1); // (1,2,1) is correct. see triangle_inequality() below!
  }

  // Multiplicity of the (Q,SS,T) subspace is (2S+1 = SS) times (2T+1).
  int mult(const Invar &I) { return I.get("SS") * (2*I.get("T")+1); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
      su2_triangle_inequality(2*I1.get("T")+1, 2*I2.get("T")+1, 2*I3.get("T")+1);
  }

  bool Invar_allowed(const Invar &I) {
    const bool spin_ok = I.get("SS") > 0;
    const bool angmom_ok = I.get("T") >= 0;
    return spin_ok && angmom_ok;
  }

  void load() {
    my_assert(!substeps);
    my_assert(channels == 3);
#include "qst/qst-In2.dat"
#include "qst/qst-QN.dat"
  } // load

  // Same as for SYMTYPE=QS, because spin operators are angular momentum singlets.
  double dynamicsusceptibility_factor(const Invar &Ip,
                                      const Invar &I1) {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "T", 0);

    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert((abs(ss1-ssp) == 2 || ss1 == ssp));

    return switch3(ss1,
                   ssp+2, 1.+(ssp-1)/3.,
                   ssp,   ssp/3.,
                   ssp-2, (-2.+ssp)/3.);
  }

  double dynamic_orb_susceptibility_factor(const Invar &Ip,
                                           const Invar &I1) {
    check_diff(Ip, I1, "Q", 0);
    check_diff(Ip, I1, "SS", 0);

    const Sspin tp = Ip.get("T");
    const Sspin t1 = I1.get("T");
    int ttp = 2*tp+1;
    int tt1 = 2*t1+1;
    my_assert((abs(tt1-ttp) == 2 || tt1 == ttp));

    return switch3(tt1,
                   ttp+2, 1.+(ttp-1)/3.,
                   ttp,   ttp/3.,
                   ttp-2, (-2.+ttp)/3.);
  }

  // Creation operator is a spin-doublet, angular-momentum-triplet !
  // See clebsch_gordan_qst.nb
  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    check_diff(Ip, I1, "Q", 1);
    
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert( abs(ss1-ssp) == 1 );
     
    double spinfactor = ( ss1 == ssp+1 ? S(ssp) + 1.0 : S(ssp) );
     
    const Tangmom tp = Ip.get("T");
    const Tangmom t1 = I1.get("T");
    const int ttp = 2*tp+1;
    const int tt1 = 2*t1+1;
    my_assert( abs(ttp-tt1) == 2 || ttp == tt1 );

    double angmomfactor;
     
    if (tt1 == 1 && ttp == 1) {
       angmomfactor = 0; // special case, two singlets don't couple
    } else {
       angmomfactor = switch3(tt1,
			      ttp+2, 1.+(ttp-1)/3.,
			      ttp,   ttp/3.,
			      ttp-2, (-2.+ttp)/3.);
    }

    return spinfactor * angmomfactor;
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ, trTZ, trQ, trQ2; // Tr[S_z^2], Tr[T_z^2], Tr[Q], Tr[Q^2]

    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Number q = I.get("Q");
      const Sspin ss = I.get("SS");
      const Tangmom t = I.get("T");
      const double sumZ = calculate_Z(is, factor);

      trQ  += sumZ * q;
      trQ2 += sumZ * q * q;
      trSZ += sumZ * (ss*ss-1)/12.; // [(2S+1)(2S+1)-1]/12=S(S+1)/3
      trTZ += sumZ * t*(t+1)/3.;
    }

    Sz2 = trSZ/STAT::Z;
    Tz2 = trTZ/STAT::Z;
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
  }

  DECL; 
  HAS_DOUBLET; HAS_TRIPLET; HAS_ORB_TRIPLET;
};

Symmetry * SymQST = new SymmetryQST;

bool qst_exception(unsigned int i,
		   unsigned int j,
		   const Invar &I)
{
  // In these cases the subspace exists, but taking the T=2 or T=1
  // limit shows that the coefficient is actually zero, so there is
  // no contribution. (Directly computed factor is nan.)
  // This exception handling is added to avoid false positives
  // in error detection assertions.
  Tangmom T = I.get("T");
  if (i == 9 && j == 34 && T == 2) 
     return true;
  if (i == 17 && j == 23 && T == 1)
     return true;
  if (i == 52 && j == 58 && T == 1)
     return true;
   
  return false;
}

#define offdiag_qst(i, j, ch, fnr, factor0, h, qq, In, I) \
   { const bool contributes = offdiag_contributes(i, j, ch, qq); \
     if (contributes) { \
        t_matel factor; \
        if (qst_exception(i, j, I)) { factor = 0; } else { factor = factor0; } \
        offdiag_build(i, j, ch, fnr, factor, h, qq, In); \
     } \
   }; 

// We take the coefficients of the first channel (indexed as 0),
// because all three set are exactly the same due to orbital
// symmetry.
#undef OFFDIAG
#define OFFDIAG(i, j, factor0) offdiag_qst(i, j, 0, 0, \
                              t_matel(factor0) * xi(STAT::N, 0), h, qq, In, I)

#undef DIAG
#define DIAG(i, number) diag_function(i, 0, number, \
                              zeta(STAT::N+1, 0), h, qq)

void SymmetryQST::makematrix(Matrix &h, const Rmaxvals &qq,
                             const Invar &I, const InvarVec &In)
{
   Sspin ss = I.get("SS");
   Tangmom t = I.get("T");
   double T = t; // crucially important to use floating point!
   
   my_assert(!substeps);
   my_assert(channels == 3);

#include "qst/qst-offdiag.dat"
#include "qst/qst-diag.dat"
}

#include "nrg-recalc-QST.cc"
