class SymmetryQS : public Symmetry 
{
private:
   outfield Sz2, Q, Q2;

public:
   SymmetryQS() : Symmetry() {
      all_syms["QS"] = this;
   }
  
   void init() {
      Sz2.set("<Sz^2>", 1);
      Q.set("<Q>", 2);
      Q2.set("<Q^2>", 3);
      InvarStructure InvStruc[] = {
	 {"Q", additive},  // charge
	 {"SS", additive}  // spin
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1);
      Invar_f = Invar(1, 2); 
   }

  // Multiplicity of the (Q,SS) subspace is 2S+1 = SS.
  int mult(const Invar &I) { return I.get("SS"); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS"));
  }

  bool Invar_allowed(const Invar &I) {
    const bool spin_ok = I.get("SS") > 0;
    return spin_ok;
  }

  void load() {
    if (!substeps) {
       switch (channels) {
	case 1:
#include "qs/qs-1ch-In2.dat"
#include "qs/qs-1ch-QN.dat"
	  break;

	case 2:
#include "qs/qs-2ch-In2.dat"
#include "qs/qs-2ch-QN.dat"
	  break;

	case 3:
#include "qs/qs-3ch-In2.dat"
#include "qs/qs-3ch-QN.dat"
	  break;

        case 4:
#include "qs/qs-4ch-In2.dat"
#include "qs/qs-4ch-QN.dat"
	  break;

	default:
	  my_assert_not_reached();
       } // switch
    } else {
#include "qs/qs-1ch-In2.dat"
#include "qs/qs-1ch-QN.dat"
    } // if
  } // load

  double dynamicsusceptibility_factor(const Invar &Ip,
                                      const Invar &I1) {
    check_diff(Ip, I1, "Q", 0);

    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert((abs(ss1-ssp) == 2 || ss1 == ssp));

    // Ce I1 in Ip singleta, ss1=ssp=1, je rezultat 1/3. Moral bi biti 0.
    // Ce I1 in Ip dubleta, ss1=ssp=2, je resultat 2/3.
    // Ce I1 in Ip tripleta, ss1=ssp=3, je rezultat 1.

    // ss1=ssp=1 je posebni primer, vendar ni tezav, ker so itak tedaj
    // tudi ireducibilni matricni elementi enaki 0. V splosnem pa je
    // potrebna posebna obravnava tega limitnega primera.
    
    // En singlet, en triplet: ss1=1, ssp=3 => 1/3; ss1=3, ssp=1 => 1. 
    // Ni simetricno, vrstni red je pomemben!
     
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
    my_assert( abs(ss1-ssp) == 1 );
    
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

  DECL; HAS_DOUBLET; HAS_TRIPLET; HAS_GLOBAL;
  HAS_SUBSTEPS;

  void show_coefficients();
};

Symmetry * SymQS = new SymmetryQS;

// *** Helper macros for makematrix() members in matrix.cc
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
			    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

/* i - subspace index
   ch - channel (0 or 1)
   number - number of electrons added in channel 'ch' in subspace 'i' */
#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
                                 zeta(STAT::N+1, ch), h, qq)

#undef OFFDIAG_MIX
#define OFFDIAG_MIX(i, j, ch, factor) offdiag_function(i, j, ch, 0, \
						       t_matel(factor) * xiR(STAT::N, ch), h, qq, In)

#undef RUNGHOP
#define RUNGHOP(i, j, factor) diag_offdiag_function(i, j, 0, \
						    t_matel(factor) * zetaR(STAT::N+1, 0), h, qq)

void SymmetryQS::makematrix(Matrix &h, const Rmaxvals &qq,
                            const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
    
  if (!substeps) {
     switch (channels) {
      case 1:
#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"
	break;
    
      case 2:
#include "qs/qs-2ch-diag.dat"
#include "qs/qs-2ch-offdiag.dat"
	if (P::rungs) {
#include "qs/qs-2ch-offdiag-mix.dat"
#include "qs/qs-2ch-runghop.dat"
	}
	break;
    
      case 3:
#include "qs/qs-3ch-diag.dat"
#include "qs/qs-3ch-offdiag.dat"
	break;

      case 4:
#include "qs/qs-4ch-diag.dat"
#include "qs/qs-4ch-offdiag.dat"
	break;
	
      default:
	my_assert_not_reached();
     }
  } else {
     my_assert(P::coeffactor == 1);
     int Ntrue, M;
     tie(Ntrue, M) = get_Ntrue_M(STAT::N);

// Here we need scale_fix, because SCALE() function is different from
// the convention for rescaling in regular two-channel cases.
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, M, 0, \
				    t_matel(factor0) * xi(Ntrue, M)/scale_fix(STAT::N), h, qq, In)

// No scale_fix here, because SCALE() is defined as it should be.
#undef DIAG
#define DIAG(i, ch, number) diag_function(i, M, number, \
	  zeta(Ntrue+1, M), h, qq)
     
#include "qs/qs-1ch-offdiag.dat"
#include "qs/qs-1ch-diag.dat"

     if (P::rungs)
         my_error("Not implemented.");
  } 
}

void SymmetryQS::show_coefficients()
{
   if (P::rungs) {
      using namespace STAT;
      for (unsigned int i = 0; i < P::channels; i++) {
	 cout << "[" << i+1 << "]"
	    << " xi_rung(" << N << ")=" << xiR(N, i)
	       << " zeta_rung(" << N+1 << ")=" << zetaR(N+1, i) << endl;
      }
   }
}

#include "nrg-recalc-QS.cc"
