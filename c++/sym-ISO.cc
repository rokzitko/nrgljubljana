class SymmetryISOcommon : public Symmetry  
{
private:
   outfield Sz2, Q2;

public:
   SymmetryISOcommon() : Symmetry() {} // see below: ISO, ISO2
   
   void init() {
      Sz2.set("<Sz^2>", 1);
      Q2.set("<Q^2>", 2);
      InvarStructure InvStruc[] = {
	 {"II", additive}, // isospin
	 {"SS", additive}  // spin
      };
      initInvar(InvStruc, ARRAYLENGTH(InvStruc));
      InvarSinglet = Invar(1, 1);
      Invar_f = Invar(2, 2);
   }

  // Multiplicity of the I=(II,SS) subspace = (2I+1)(2S+1) = II SS.
  int mult(const Invar &I) {
    int mi = I.get("II"); // isospin multiplicity
    int ms = I.get("SS"); // spin multiplicity
    return mi * ms;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      su2_triangle_inequality(I1.get("SS"), I2.get("SS"), I3.get("SS")) &&
      su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"));
  }

  // We always must have S >= 0 and I >= 0.
  bool Invar_allowed(const Invar &I) {
    const bool isospin_ok = I.get("II") > 0;
    const bool spin_ok = I.get("SS") > 0;
    return isospin_ok && spin_ok;
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

class SymmetryISO : public SymmetryISOcommon {
public:
 SymmetryISO() : SymmetryISOcommon() {
     all_syms["ISO"] = this;
 };

 void load()
 {
   switch (channels) {
    case 1:
#include "iso/iso-1ch-In2.dat"
#include "iso/iso-1ch-QN.dat"
    break;

  case 2:
#include "iso/iso-2ch-In2.dat"
#include "iso/iso-2ch-QN.dat"
    break;

  case 3:
#include "iso/iso-3ch-In2.dat"
#include "iso/iso-3ch-QN.dat"
    break;

  default:
    my_assert_not_reached();
  }
 }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

class SymmetryISO2 : public SymmetryISOcommon {
public:
 SymmetryISO2() : SymmetryISOcommon() {
        all_syms["ISO2"] = this;
  };

 void load()
 {
   switch (channels) {
    case 1:
#include "iso2/iso2-1ch-In2.dat"
#include "iso2/iso2-1ch-QN.dat"
    break;

  case 2:
#include "iso2/iso2-2ch-In2.dat"
#include "iso2/iso2-2ch-QN.dat"
    break;

  default:
    my_assert_not_reached();
  }
 }

  DECL; HAS_DOUBLET; HAS_TRIPLET;
};

Symmetry * SymISO = new SymmetryISO;
Symmetry * SymISO2 = new SymmetryISO2;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
					  zeta(STAT::N+1, ch), h, qq)

void SymmetryISO::makematrix(Matrix &h, const Rmaxvals &qq,
                             const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");
  // nn is the index of the last site in the chain, while nn+1 is the
  // index of the site that is being added to the chain in this
  // iteration.  This is consistent with the definition in
  // isospin-1ch-automatic.nb.
  int NN = getnn();

  switch (channels) {
  case 1:
#include "iso/iso-1ch-offdiag.dat"
    break;

  case 2:
#include "iso/iso-2ch-offdiag.dat"
    break;

  case 3:
#include "iso/iso-3ch-offdiag.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

void SymmetryISO2::makematrix(Matrix &h, const Rmaxvals &qq,
                              const Invar &I, const InvarVec &In)
{
  Sspin ss = I.get("SS");
  Ispin ii = I.get("II");
  int NN = getnn();

  switch (channels) {
  case 1:
#include "iso2/iso2-1ch-offdiag.dat"
    break;

  case 2:
#include "iso2/iso2-2ch-offdiag.dat"
    break;

  default:
    my_assert_not_reached();
  }
}

#include "nrg-recalc-ISO.cc"
#include "nrg-recalc-ISO2.cc"
