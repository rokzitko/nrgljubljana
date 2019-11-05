class SymmetrySL3 : public Symmetry
{
private:
   outfield Q1, Q12, sQ12, Q2, Q22, sQ22, Q3, Q32, sQ32;

public:
   SymmetrySL3() : Symmetry() {
      all_syms["SL3"] = this;
   }
  
   void init() {
      Q1.set("<Q1>", 1); // charge
      Q12.set("<Q1^2>", 2); // charge squared
      sQ12.set("<sQ1^2>", 3); // charge fluctuation
      Q2.set("<Q2>", 4);
      Q22.set("<Q2^2>", 5);
      sQ22.set("<sQ2^2>", 6);
      Q3.set("<Q3>", 7);
      Q32.set("<Q3^2>", 8);
      sQ32.set("<sQ3^2>", 9);
      InvarStructure InvStruc[] = {
	 {"Q1", additive}, // charge in channel 1
	 {"Q2", additive}, // charge in channel 2
	 {"Q3", additive}  // charge in channel 3
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 0, 0);
   }
  
  void load() {
    switch (channels) {
    case 3:
#include "sl3/sl3-3ch-In2.dat"
#include "sl3/sl3-3ch-QN.dat"
      break;

    default:
      my_assert_not_reached();
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trQ1, trQ12; // Tr[Q], Tr[Q^2]
    bucket trQ2, trQ22;
    bucket trQ3, trQ32;
     
    LOOP_const(diag, is) {
      const Number q1 = INVAR(is).get("Q1");
      const Number q2 = INVAR(is).get("Q2");
      const Number q3 = INVAR(is).get("Q3");
      const double sumZ = calculate_Z(is, factor);

      trQ1  += sumZ * q1;
      trQ12 += sumZ * q1 * q1;
      trQ2  += sumZ * q2;
      trQ22 += sumZ * q2 * q2;
      trQ3  += sumZ * q3;
      trQ32 += sumZ * q3 * q3;
    }
    
    Q1 = trQ1/STAT::Z;
    Q12 = trQ12/STAT::Z;
    // charge fluctuations -> susceptibility
    sQ12 = (trQ12/STAT::Z) - pow(trQ1/STAT::Z, 2);

    Q2 = trQ2/STAT::Z;
    Q22 = trQ22/STAT::Z;
    // charge fluctuations -> susceptibility
    sQ22 = (trQ22/STAT::Z) - pow(trQ2/STAT::Z, 2);

    Q3 = trQ3/STAT::Z;
    Q32 = trQ32/STAT::Z;
    // charge fluctuations -> susceptibility
    sQ32 = (trQ32/STAT::Z) - pow(trQ3/STAT::Z, 2);
}

  DECL; HAS_DOUBLET; HAS_GLOBAL;
};

Symmetry * SymSL3 = new SymmetrySL3;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
                                 zeta(STAT::N+1, ch), h, qq)

void SymmetrySL3::makematrix(Matrix &h, const Rmaxvals &qq,
                            const Invar &I, const InvarVec &In)
{
  switch (channels) {
  case 3:
#include "sl3/sl3-3ch-offdiag.dat"
#include "sl3/sl3-3ch-diag.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

#include "nrg-recalc-SL3.cc"
