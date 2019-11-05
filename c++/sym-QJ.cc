class SymmetryQJ : public Symmetry
{
private:
   outfield Jz2, Q, Q2;

public:
   SymmetryQJ() : Symmetry() {
      all_syms["QJ"] = this;
   }
  
   void init() {
      Jz2.set("<Jz^2>", 1);
      Q.set("<Q>", 2);
      Q2.set("<Q^2>", 3);
      InvarStructure InvStruc[] = {
	 {"Q", additive},  // charge
	 {"JJ", additive}  // total angular momentum
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(0, 1);
      Invar_f = Invar(1, 4);
   }

  // Multiplicity of the (Q,JJ) subspace is 2J+1 = JJ.
  int mult(const Invar &I) { return I.get("JJ"); }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      su2_triangle_inequality(I1.get("JJ"), I2.get("JJ"), I3.get("JJ"));
  }

  bool Invar_allowed(const Invar &I) {
    const bool spin_ok = I.get("JJ") > 0;
    return spin_ok;
  }

  void load() {
#include "qj/qj-In2.dat"
#include "qj/qj-QN.dat"
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trJZ2, trQ, trQ2; // Tr[J_z^2], Tr[Q], Tr[Q^2]

    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Sspin jj = I.get("JJ");
      const Number q = I.get("Q");
      const double sumZ = calculate_Z(is, factor);

      trQ  += sumZ * q;
      trQ2 += sumZ * q * q;
      trJZ2 += sumZ * (jj*jj-1)/12.;
    }

    Jz2 = trJZ2/STAT::Z;
    Q = trQ/STAT::Z;
    Q2 = trQ2/STAT::Z;
  }

  // ClebschGordan[ket (p), op, bra (1)]
  double specdens_factor(const Invar &Ip,
			 const Invar &I1) {
     check_diff(Ip, I1, "Q", 1);
     
     const Sspin jjp = Ip.get("JJ");
     const Sspin jj1 = I1.get("JJ");
     my_assert( abs(jj1-jjp) == 1 );
     
     return ( jj1 == jjp+1 ? S(jjp) + 1.0 : S(jjp) );
  }

  // See cg_factors_doublet_triplet_quadruplet.nb
  double specdensquad_factor(const Invar &Ip,
			 const Invar &I1) {
     check_diff(Ip, I1, "Q", 1);
     
     const Sspin jjp = Ip.get("JJ");
     const Sspin jj1 = I1.get("JJ");
     my_assert( abs(jj1-jjp) == 1 || abs(jj1-jjp) == 3);

     if (jj1 == jjp + 3)
	return S(jjp)/2.0 + 1.0;
     
     if (jj1 == jjp + 1) {
	my_assert(jjp >= 2); // singlet in kvadruplet se ne moreta sklopiti v dublet
	return S(jjp)/2.0 + 0.5;
     }
     
     if (jj1 == jjp - 1) {
	my_assert(jjp >= 2); // trikotniska
	return S(jjp)/2.0;
     }
     
     if (jj1 == jjp - 3) {
	my_assert(jjp >= 4); // trikotniska
	return S(jjp)/2.0 - 0.5;
     }
     
     my_assert_not_reached();
     return 0;
  }

  void offdiag_function_QJ(unsigned int i,
			    unsigned int j,
			    unsigned int ch,
			    unsigned int fnr,
			    t_matel factor,
			    Matrix &h,
			    const Rmaxvals &qq,
			    const InvarVec &In);

  HAS_DOUBLET; HAS_QUADRUPLET;
  DECL;
};

Symmetry * SymQJ = new SymmetryQJ;

#include <boost/math/special_functions/factorials.hpp>

double Factorial(double x)
{
   return boost::math::factorial<double>(round(x));
}

void SymmetryQJ::offdiag_function_QJ(unsigned int i,
                      unsigned int j,
                      unsigned int ch, // channel number
                      unsigned int fnr, // extra index for <||f||>, usually 0
                      t_matel factor, // may be complex (in principle)
                      Matrix &h,
                      const Rmaxvals &qq,
                      const InvarVec &In)
{
  const Invar Iop = (ch == 0 ? Invar(1, 2) : Invar(1, 4));
  const Invar I1 = In[i];
  const Invar I2 = In[j];
  const bool triangle = triangle_inequality(I1, I2, Iop); // I1 = I2+Iop

  const bool contributes = offdiag_contributes(i, j, ch, qq);

  if (triangle && contributes) {
    offdiag_build(i, j, ch, fnr, factor, h, qq, In);
  }
}

// *** Helper macros for makematrix() members in matrix.cc
// Jndx = 0 for doublet, Jndx = 1 for quadruplet
#undef OFFDIAG
#define OFFDIAG(i, j, Jndx, factor0) offdiag_function_QJ(i, j, Jndx, 0, \
							 t_matel(factor0) * xi(STAT::N, 0), h, qq, In)

#undef DIAG
#define DIAG(i, number) diag_function(i, 0, number, \
                                 zeta(STAT::N+1, 0), h, qq)

inline double J(int JJ)
{
   return (JJ-1.0)/2.0; // JJ=2J+1
}

void SymmetryQJ::makematrix(Matrix &h, const Rmaxvals &qq,
                            const Invar &I, const InvarVec &In)
{
  Sspin jj = I.get("JJ");
    
#include "qj/qj-offdiag.dat"
#include "qj/qj-diag.dat"
}

#include "nrg-recalc-QJ.cc"
