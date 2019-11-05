class SymmetryQSZLR : public SymFieldLR
{
private:
   outfield Sz2, Sz, Q, Q2;
   
public:
   SymmetryQSZLR() : SymFieldLR() {
      all_syms["QSZLR"] = this;
   }
   
  void init() { 
     Sz2.set("<Sz^2>", 1);
     Sz.set("<Sz>", 2);
     Q.set("<Q>", 3);
     Q2.set("<Q^2>", 4);
     InvarStructure InvStruc[] = {
	{"Q", additive},      // charge 
	{"SSZ", additive},    // spin projection 
	{"P", multiplicative} // parity 
     };
     initInvar(InvStruc, ARRAYLENGTH(InvStruc));
     InvarSinglet = Invar(0, 0, 1);
  }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    const SZspin ssz1 = I1.get("SSZ");
    const SZspin sszp = Ip.get("SSZ");
    const SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }
  
  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) {
    return
      u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")) &&
      u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) &&
      z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() {
    my_assert(channels == 2);
#include "qszlr/qszlr-2ch-In2.dat"
#include "qszlr/qszlr-2ch-QN.dat"
  }

  void calculate_TD(const DiagInfo &diag, double factor) {
    bucket trSZ, trSZ2, trQ, trQ2; // Tr[S_z], Tr[(S_z)^2], etc.
    
    LOOP_const(diag, is) {
      const Invar I = INVAR(is);
      const SZspin ssz = I.get("SSZ");
      const Number q = I.get("Q");
      const double sumZ = calculate_Z(is, factor);

      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * sqr( SZ(ssz) );
      trQ  += sumZ * q;
      trQ2 += sumZ * sqr( q );
    }

    Sz2= trSZ2/STAT::Z;
    Sz = trSZ/STAT::Z;
    Q = trQ/STAT::Z;
    Q2= trQ2/STAT::Z;
  }

  DECL;
};

Symmetry * SymQSZLR = new SymmetryQSZLR;

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, \
						    t_matel(factor0) * xi(STAT::N, ch), h, qq, In)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, \
                                 zeta(STAT::N+1, ch), h, qq)

void SymmetryQSZLR::makematrix(Matrix &h, const Rmaxvals &qq,
                               const Invar &I, const InvarVec &In)
{
#include "qszlr/qszlr-2ch-offdiag.dat"
#include "qszlr/qszlr-2ch-diag.dat"
}

#include "nrg-recalc-QSZLR.cc"
