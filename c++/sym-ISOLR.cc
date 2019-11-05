class SymmetryISOLRcommon : public SymLR
{
private:
   outfield Sz2, Q2;
  
public:
   SymmetryISOLRcommon() : SymLR() {} // see below: ISOLR and ISO2LR

   void init() {
      Sz2.set("<Sz^2>", 1);
      Q2.set("<Q^2>", 2);
      InvarStructure InvStruc[] = {
	 {"II", additive},     // isospin
	 {"SS", additive},     // spin
	 {"P", multiplicative} // parity
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));    
      InvarSinglet = Invar(1, 1, 1);
      //     Invar_f = Invar(?, ?, ?)
   }
   
  // Multiplicity of the I=(II,SS,P) subspace = (2I+1)(2S+1) = II SS.
  int mult(const Invar &I) {
    int mi = I.get("II"); // isospin multiplicity
    int ms = I.get("SS"); // spin multiplicity
    return mi * ms;
  }

  // We always must have S >= 0 and I >= 0.
  bool Invar_allowed(const Invar &I) {
    const bool isospin_ok = I.get("II") > 0;
    const bool spin_ok = I.get("SS") > 0;
    return isospin_ok && spin_ok;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II")) &&
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
      z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  double specdens_factor(const Invar &Ip,
                         const Invar &I1) {
    const Sspin ssp = Ip.get("SS");
    const Sspin ss1 = I1.get("SS");
    my_assert( abs(ss1-ssp) == 1 );
    
    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");
    my_assert( abs(ii1-iip) == 1 );
    
    const double spinfactor = ( ss1 == ssp+1 ? S(ssp) + 1.0   : S(ssp) );
    const double isofactor =  ( ii1 == iip+1 ? ISO(iip) + 1.0 : ISO(iip) );
    
    return spinfactor * isofactor;
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ, trIZ; // Tr[S_z^2], Tr[I_z^2]
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const Ispin ii = I.get("II");
      const Sspin ss = I.get("SS");
      const double sumZ = calculate_Z(is, factor);

      trSZ += sumZ * (ss*ss-1)/12.; // isospin multiplicity contained in sumZ
      trIZ += sumZ * (ii*ii-1)/12.; // spin multiplicity contained in sumZ
    }

    Sz2 = trSZ/STAT::Z;
    Q2 = (4*trIZ)/STAT::Z;
  }
};

class SymmetryISOLR : public SymmetryISOLRcommon
{
public:
  SymmetryISOLR() : SymmetryISOLRcommon() {
       all_syms["ISOLR"] = this;
  };

  void load()
  {
    my_assert(channels == 2);
#include "isolr/isolr-2ch-In2.dat"
#include "isolr/isolr-2ch-QN.dat"
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

class SymmetryISO2LR : public SymmetryISOLRcommon
{
public:
  SymmetryISO2LR() : SymmetryISOLRcommon() {
      all_syms["ISO2LR"] = this;
  };

  void load()
  {
    my_assert(channels == 2);
#include "iso2lr/iso2lr-2ch-In2.dat"
#include "iso2lr/iso2lr-2ch-QN.dat"
  }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymISOLR = new SymmetryISOLR;
Symmetry * SymISO2LR = new SymmetryISO2LR;

// *** Helper macros for makematrix() members in matrix.cc
#undef OFFIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)
                              
void SymmetryISOLR::makematrix(Matrix &h, const Rmaxvals &qq,
                               const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");

#include "isolr/isolr-2ch-offdiag.dat"
}

void SymmetryISO2LR::makematrix(Matrix &h, const Rmaxvals &qq,
                                const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");
  
#include "iso2lr/iso2lr-2ch-offdiag.dat"
}

#include "nrg-recalc-ISOLR.cc"
#include "nrg-recalc-ISO2LR.cc"
