class SymmetrySPSU2 : public Symmetry
{
private:
   outfield Sz2;

public:
   SymmetrySPSU2() : Symmetry() {
      all_syms["SPSU2"] = this;
   }
   
   void init() {
      Sz2.set("<Sz^2>", 1);
      InvarStructure InvStruc[] = {
	 {"SS", additive} // spin
      };     
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(1);
   }

  // Multiplicity of the I=(SS) subspace = 2S+1 = SS.
  int mult(const Invar &I) {
    return I.get("SS");
  }

  bool Invar_allowed(const Invar &I) {
    return I.get("SS") > 0;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"));
  }

  void load() {
   if (!substeps) {
    switch(channels) {
    case 1:
#include "spsu2/spsu2-1ch-In2.dat"
#include "spsu2/spsu2-1ch-QN.dat"
      break;

    case 2:
#include "spsu2/spsu2-2ch-In2.dat"
#include "spsu2/spsu2-2ch-QN.dat"
      break;

    case 3:
#include "spsu2/spsu2-3ch-In2.dat"
#include "spsu2/spsu2-3ch-QN.dat"
      break;

    default:
      my_assert_not_reached();
    }
   } else {
#include "spsu2/spsu2-1ch-In2.dat"
#include "spsu2/spsu2-1ch-QN.dat"
   }
  }

  double dynamicsusceptibility_factor(const Invar &ground,
                                      const Invar &I1) {
    const Sspin ssp = ground.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert( abs(ss1-ssp) == 2 || ss1 == ssp );

    return switch3(ss1,
                   ssp+2, 1.+(ssp-1)/3.,
                   ssp,   ssp/3.,
                   ssp-2, (-2.+ssp)/3.);
  }

  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert( abs(ss1-ssp) == 1 );

    return ( ss1 == ssp+1 ? S(ssp) + 1.0 : S(ssp) );
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trSZ; // Tr[S_z^2]
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Sspin ss = I.get("SS");
      const double sumZ = calculate_Z(is, factor);
      
      trSZ += sumZ * (ss*ss-1)/12.;
    }
    
    Sz2 = trSZ/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET; HAS_GLOBAL;
  HAS_SUBSTEPS;
};

Symmetry * SymSPSU2 = new SymmetrySPSU2;

// NOTE: the additional factor of 2 is due to the fact that the
// isospin is in fact defined as 1/2 of the (d^\dag d^\dag + d d)
// pairing operator.

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, \
   t_matel(factor) * 2.0 * delta(STAT::N+1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, ch, 0, \
                        t_matel(factor) * kappa(STAT::N, ch), \
                        h, qq, In)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
			  zeta(STAT::N+1, ch), h, qq)

void SymmetrySPSU2::makematrix(Matrix &h, const Rmaxvals &qq,
                               const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");

  if (!substeps) {
   switch (channels) {
    case 1:
#include "spsu2/spsu2-1ch-offdiag.dat"
#include "spsu2/spsu2-1ch-anomalous.dat"
#include "spsu2/spsu2-1ch-diag.dat"
#include "spsu2/spsu2-1ch-isospinx.dat"
     break;
    
    case 2:
#include "spsu2/spsu2-2ch-diag.dat"
#include "spsu2/spsu2-2ch-offdiag.dat"
#include "spsu2/spsu2-2ch-anomalous.dat"
#include "spsu2/spsu2-2ch-isospinx.dat"
     break;

   case 3:
#include "spsu2/spsu2-3ch-diag.dat"
#include "spsu2/spsu2-3ch-offdiag.dat"
#include "spsu2/spsu2-3ch-anomalous.dat"
#include "spsu2/spsu2-3ch-isospinx.dat"
     break;

    default:
     my_assert_not_reached();
   }
  } else {
     my_assert(P::coeffactor == 1);
     int Ntrue, M;
     tie(Ntrue, M) = get_Ntrue_M(STAT::N);
      
// Careful: for historical reasons (argh!) there is a minus sign in
// the coefficients for the second Wilson chain; see also function
// matrixopisospin[] in coefnew/spsu2/spsu2.m. That's why we need a
// (-1) factor for M==1.
// (M == 1 ? -1 : 1)
#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, M, \
  t_matel(factor) * 2.0 * (M == 1 ? -1.0 : 1.0) * \
  delta(Ntrue+1, M), h, qq)
     
#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, M, 0, \
  t_matel(factor) * kappa(Ntrue, M), h, qq, In)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, M, 0, \
  t_matel(factor0) * xi(Ntrue, M)/scale_fix(STAT::N), h, qq, In)
      
#undef DIAG
#define DIAG(i, ch, number) diag_function(i, M, number, \
  zeta(Ntrue+1, M), h, qq)
   
#include "spsu2/spsu2-1ch-offdiag.dat"
#include "spsu2/spsu2-1ch-anomalous.dat"
#include "spsu2/spsu2-1ch-diag.dat"
#include "spsu2/spsu2-1ch-isospinx.dat"
      
  }
}

#include "nrg-recalc-SPSU2.cc"
