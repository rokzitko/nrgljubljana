class SymmetryDBLSU2 : public Symmetry
{
private:
   outfield Q12, Q22;

public:
   SymmetryDBLSU2() : Symmetry() {
      all_syms["DBLSU2"] = this;
   }
  
  void init() {
     Q12.set("<Q1^2>", 1);
     Q22.set("<Q2^2>", 2);
     InvarStructure InvStruc[] = {
	{"II1", additive}, // isospin 1
	{"II2", additive}, // isospin 2
     };
     initInvar(InvStruc, ARRAYLENGTH(InvStruc));
     InvarSinglet = Invar(1, 1);
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return 
       su2_triangle_inequality(I1.get("II1"), I2.get("II1"), I3.get("II1")) &&
       su2_triangle_inequality(I1.get("II2"), I2.get("II2"), I3.get("II2"));
  }

   // Multiplicity of the I=(II1,II2) subspace = (2I1+1)(2I2+1) = II1 II2.
  int mult(const Invar &I) {
    return I.get("II1") * I.get("II2");
  }
  
  // We always must have I1 >= 0 and I2 >= 0.
  bool Invar_allowed(const Invar &I) {
    return (I.get("II1") > 0) && (I.get("II2") > 0);
  }

  // TO DO: support for the doublets wrt the second quantum number 
  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    const Ispin ii1p = Ip.get("II1");
    const Ispin ii11 = I1.get("II1");
    my_assert( abs(ii11-ii1p) == 1 );
    
    const double isofactor =  ( ii11 == ii1p+1 ? ISO(ii1p) + 1.0 : ISO(ii1p) );
    
    return isofactor;
  }

  void load() {
    switch (channels) {
    case 2:
#include "dblsu2/dblsu2-2ch-In2.dat"
#include "dblsu2/dblsu2-2ch-QN.dat"
      break;

    default:
      my_assert_not_reached();
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trIZ12; // Tr[I1_z^2]
    bucket trIZ22; // Tr[I2_z^2]
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Number ii1 = I.get("II1");
      const Number ii2 = I.get("II2");
      const double sumZ = calculate_Z(is, factor);

      trIZ12 += sumZ * (ii1*ii1-1)/12.;
      trIZ22 += sumZ * (ii2*ii2-1)/12.;
    }

    Q12 = (4*trIZ12)/STAT::Z;
    Q22 = (4*trIZ22)/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_GLOBAL;
};

Symmetry * SymDBLSU2 = new SymmetryDBLSU2;

// see sym-SU2.cc

#undef OFFDIAG_1
#define OFFDIAG_1(i, j, ch, factor) offdiag_function(i, j, ch, 0, \
						     t_matel(factor) * xi(STAT::N, ch), h, qq, In)
#undef OFFDIAG_2
#define OFFDIAG_2(i, j, ch, factor) offdiag_function(i, j, ch, 1, \
						     t_matel(factor) * xi(STAT::N, ch), h, qq, In)
	      
void SymmetryDBLSU2::makematrix(Matrix &h, const Rmaxvals &qq,
                                const Invar &I, const InvarVec &In)
{
  switch (channels) {
  case 2:
#include "dblsu2/dblsu2-2ch-offdiag-1.dat"
#include "dblsu2/dblsu2-2ch-offdiag-2.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

#include "nrg-recalc-DBLSU2.cc"
