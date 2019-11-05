class SymmetrySPSU2LR : public SymLR
{
private:
   outfield Sz2;
   
public:
   SymmetrySPSU2LR() : SymLR() {
      all_syms["SPSU2LR"] = this;
   }
   
   void init() {
      Sz2.set("<Sz^2>", 1);
      InvarStructure InvStruc[] = {
	 {"SS", additive}, // spin
	 {"P", multiplicative} // parity
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1);
   }
     
  int mult(const Invar &I) {
     return I.get("SS");
  }
   
  bool Invar_allowed(const Invar &I) {
     return I.get("SS") > 0;
  }
   
  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
         return su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
	        z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }
  
  void load() {
     my_assert(channels == 2);
#include "spsu2lr/spsu2lr-2ch-In2.dat"
#include "spsu2lr/spsu2lr-2ch-QN.dat"
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
   
  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trSZ2; // Tr[S_z^2]
  
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Sspin ss = I.get("SS");
      const double sumZ = calculate_Z(is, factor);
      
      trSZ2 += sumZ * (ss*ss-1)/12.;
    }
    
    Sz2 = trSZ2/STAT::Z;
  }  

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymSPSU2LR = new SymmetrySPSU2LR;

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, \
 t_matel(factor) * 2.0 * delta(STAT::N+1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, ch, 0, \
 t_matel(factor) * kappa(STAT::N, ch), h, qq, In)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
 zeta(STAT::N+1, ch), h, qq)

void SymmetrySPSU2LR::makematrix(Matrix &h, const Rmaxvals &qq,
                              const Invar &I, const InvarVec &In)
{
   my_assert(channels == 2);
   Sspin ss = I.get("SS");
#include "spsu2lr/spsu2lr-2ch-diag.dat"
#include "spsu2lr/spsu2lr-2ch-offdiag.dat"
#include "spsu2lr/spsu2lr-2ch-anomalous.dat"
#include "spsu2lr/spsu2lr-2ch-isospinx.dat"
}

#include "nrg-recalc-SPSU2LR.cc"
