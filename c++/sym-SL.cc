class SymmetrySL : public Symmetry
{
private:
   outfield Q, Q2, sQ2;

public:
   SymmetrySL() : Symmetry() {
      all_syms["SL"] = this;
   }
  
   void init() {
      Q.set("<Q>", 1);
      Q2.set("<Q^2>", 2);
      sQ2.set("<sQ^2>", 3);
      InvarStructure InvStruc[] = {
	 {"Q", additive} // charge
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0);
   }
  
  void load() {
    switch (channels) {
    case 1:
#include "sl/sl-1ch-In2.dat"
#include "sl/sl-1ch-QN.dat"
      break;
      
    case 2:
#include "sl/sl-2ch-In2.dat"
#include "sl/sl-2ch-QN.dat"
      break;

    case 3:
#include "sl/sl-3ch-In2.dat"
#include "sl/sl-3ch-QN.dat"
      break;

    default:
      my_assert_not_reached();
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
    
    LOOP_const(diag, is) {
      const Number q = INVAR(is).get("Q");
      const double sumZ = calculate_Z(is, factor);

      trQ  += sumZ * q;
      trQ2 += sumZ * q * q;
    }
    
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
    // charge fluctuations -> susceptibility
    sQ2 = (trQ2/STAT::Z) - pow(trQ/STAT::Z, 2);
  }

  DECL; HAS_DOUBLET; HAS_GLOBAL;
};

Symmetry * SymSL = new SymmetrySL;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
 zeta(STAT::N+1, ch), h, qq)

void SymmetrySL::makematrix(Matrix &h, const Rmaxvals &qq,
                            const Invar &I, const InvarVec &In)
{
  switch (channels) {
  case 1:
#include "sl/sl-1ch-offdiag.dat"
#include "sl/sl-1ch-diag.dat"
    break;

  case 2:
#include "sl/sl-2ch-offdiag.dat"
#include "sl/sl-2ch-diag.dat"
    break;

  case 3:
#include "sl/sl-3ch-offdiag.dat"
#include "sl/sl-3ch-diag.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

#include "nrg-recalc-SL.cc"
