class SymmetryISOSZLR : public SymFieldLR
{
private:
  outfield Sz2, Sz, Q2;

public:
   SymmetryISOSZLR() : SymFieldLR() {
       all_syms["ISOSZLR"] = this;
  }

  void init() {
     Sz2.set("<Sz^2>", 1);
     Sz.set("<Sz>", 2);
     Q2.set("<Q^2>", 3);
     InvarStructure InvStruc[] = {
       {"II", additive},     // isospin 
       {"SSZ", additive},    // spin projection 
       {"P", multiplicative} // parity 
     };
     initInvar(InvStruc, ARRAYLENGTH(InvStruc));
     InvarSinglet = Invar(1, 0, 1);
  }

  // Multiplicity of the I=(II,SSZ) subspace = (2I+1) = II.
  int mult(const Invar &I) {
    return I.get("II");
  }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1 = I1.get("SSZ");
    SZspin sszp = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) &&
      su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II")) &&
      z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() {
     my_assert(channels == 2);
	
#include "isoszlr/isoszlr-2ch-In2.dat"
#include "isoszlr/isoszlr-2ch-QN.dat"
  }


  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    check_abs_diff(Ip, I1, "SSZ", 1);

    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");

    const double isofactor =  ( ii1 == iip+1 ? ISO(iip) + 1.0 : ISO(iip) );
    return isofactor;
  }

  void calculate_TD(const DiagInfo &diag, double factor)
  {
    bucket trSZ, trSZ2, trIZ2; // Tr[S_z], Tr[S_z^2], Tr[I_z^2]
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Ispin ii = I.get("II");
      const SZspin ssz = I.get("SSZ");
      const double sumZ = calculate_Z(is, factor);

      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * sqr( SZ(ssz) ); // isospin multiplicity contained in sumZ
      trIZ2 += sumZ * (ii*ii-1)/12.; // spin multiplicity contained in sumZ
    }
    
    Sz = trSZ/STAT::Z;
    Sz2 = trSZ2/STAT::Z;
    Q2 = (4*trIZ2)/STAT::Z;
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymISOSZLR = new SymmetryISOSZLR;
     
#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

void SymmetryISOSZLR::makematrix(Matrix &h, const Rmaxvals &qq,
                                 const Invar &I, const InvarVec &In)
{
  Ispin ii = I.get("II");

#include "isoszlr/isoszlr-2ch-offdiag.dat"
}

#include "nrg-recalc-ISOSZLR.cc"
