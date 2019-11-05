class SymmetrySU2 : public Symmetry
{
private:
   outfield Q2;

public:
   SymmetrySU2() : Symmetry() {
      all_syms["SU2"] = this;
   }
   
   void init() {
      Q2.set("<Q^2>", 1);
      InvarStructure InvStruc[] = {
	 {"II", additive} // isospin
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(1);
   }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"));
  }

   // Multiplicity of the I=(II) subspace = (2I+1) = II.
  int mult(const Invar &I) {
    return I.get("II"); // isospin multiplicity
  }
  
  // We always must have I >= 0.
  bool Invar_allowed(const Invar &I) {
    return  I.get("II") > 0;
  }

  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");
    my_assert( abs(ii1-iip) == 1 );
    
    const double isofactor =  ( ii1 == iip+1 ? ISO(iip) + 1.0 : ISO(iip) );
    
    return isofactor;
  }

  void load() {
    switch (channels) {
    case 1:
#include "su2/su2-1ch-In2.dat"
#include "su2/su2-1ch-QN.dat"
      break;

    case 2:
#include "su2/su2-2ch-In2.dat"
#include "su2/su2-2ch-QN.dat"
      break;

    default:
      my_assert_not_reached();
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trIZ2; // Tr[I_z^2]
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Number ii = I.get("II");
      const double sumZ = calculate_Z(is, factor);

      trIZ2 += sumZ * (ii*ii-1)/12.;
    }

    Q2 = (4*trIZ2)/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_GLOBAL;
};

Symmetry * SymSU2 = new SymmetrySU2;

// For SU2, the <||f||> depend on the "type"!
// OFFDIAG_1 corresponds to the first combination, OFFDIAG_2 to the
// second combination of operators. Each of them has contributions
// in each channel.

#undef OFFDIAG_1
#define OFFDIAG_1(i, j, ch, factor) offdiag_function(i, j, ch, 0, \
						     t_matel(factor) * xi(STAT::N, ch), h, qq, In)
#undef OFFDIAG_2
#define OFFDIAG_2(i, j, ch, factor) offdiag_function(i, j, ch, 1, \
						     t_matel(factor) * xi(STAT::N, ch), h, qq, In)

void SymmetrySU2::makematrix(Matrix &h, const Rmaxvals &qq,
                            const Invar &I, const InvarVec &In)
{
  Ispin ii = I.get("II");
  int NN = getnn();

  switch (channels) {
   case 1:
#include "su2/su2-1ch-offdiag-1.dat"
#include "su2/su2-1ch-offdiag-2.dat"
     break;

  case 2:
#include "su2/su2-2ch-offdiag-1.dat"
#include "su2/su2-2ch-offdiag-2.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

#include "nrg-recalc-SU2.cc"
