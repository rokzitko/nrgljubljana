class SymmetryANYJ : public SymField
{
private:
   outfield Sz, Sz2, Q, Q2;

public:
   SymmetryANYJ() : SymField() {
      all_syms["ANYJ"] = this;
   }
  
   void init() {
      Sz.set("<Sz>", 1);
      Sz2.set("<Sz^2>", 2);
      Q.set("<Q>", 3);
      Q2.set("<Q^2>", 4);
      InvarStructure InvStruc[] = {
	 {"Q", additive},  // charge 
	 {"SSZ", additive} // spin projection
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 0);
   }

  void load() {
    my_assert(channels == 1);

    switch(spin) {
    case 2:
#include "anyj/any-n2-In2.dat"
      break;

    case 3:
#include "anyj/any-n3-In2.dat"
      break;
      
    default:
      my_error("Spin N not implemented.");
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[S_z^2], etc.

    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const double sz = rational_cast<double>( I.get_frac("SSZ") );
      const double q = I.get("Q");
      const double sumZ = calculate_Z(is, factor);

      trSZ += sumZ * sz;
      trSZ2 += sumZ * sqr( sz );
      trQ  += sumZ * q;
      trQ2 += sumZ * sqr( q );
    }
    
    Sz2 = trSZ2/STAT::Z;
    Sz = trSZ/STAT::Z;
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
  }

  DECL; HAS_DOUBLET;
};

Symmetry * SymANYJ = new SymmetryANYJ;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
				t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
                                 zeta(STAT::N+1, ch), h, qq)

void SymmetryANYJ::makematrix(Matrix &h, const Rmaxvals &qq,
                              const Invar &I, const InvarVec &In)
{
  switch(spin) {
  case 2:
#include "anyj/any-n2-offdiag.dat"
#include "anyj/any-n2-diag.dat"
    break;
    
  case 3:
#include "anyj/any-n3-offdiag.dat"
#include "anyj/any-n3-diag.dat"
    break;
    
  default:
    my_error("Not implemented.");
  }
}

#include "nrg-recalc-ANYJ.cc"
