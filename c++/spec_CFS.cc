// Choose one of the following two! OLD is the non-optimized code by
// Rok Zitko, OPTIMIZED is the hand-tuned code contributed by Markus
// Greger. The optimized code is faster by an order of magnitude!

class SPEC_CFSls : virtual public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "CFSls"; }
   string merge() { return "CFS"; }
};

class SPEC_CFSgt : virtual public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "CFSgt"; }
   string merge() { return "CFS"; }
};

class SPEC_CFS : public SPEC_CFSls, public SPEC_CFSgt
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen & a1, const Eigen & a2, const Matrix & a3, const Matrix & a4,
	     const BaseSpectrum & a5, t_factor a6, ChainSpectrum *a7, const Invar & a8, const Invar & a9)
   {
      SPEC_CFSgt::calc(a1,a2,a3,a4,a5,a6,a7,a8,a9);
      SPEC_CFSls::calc(a1,a2,a3,a4,a5,a6,a7,a8,a9);
   }
   string name() { return "CFS"; }
   string merge() { return "CFS"; }
};

//#define SPEC_CFS_OLD
#define SPEC_CFS_OPTIMIZED

// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).

#if defined(NRG_COMPLEX) || defined(SPEC_CFS_OLD)
void SPEC_CFSls::calc(const Eigen & diagIp, const Eigen & diagI1, 
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs,
		      t_factor spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ip, const Invar &I1)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   using namespace STAT;
   const Matrix & rhoNIp = rho[Ip]; // hand optimised out of the loops
   const Matrix & rhoNI1 = rho[I1];
   auto dimp = rhoNIp.size1();
   auto dim1 = rhoNI1.size1();
   // Convention: k-loops over retained states, l-loop over discarded
   // states.
   // i-term, Eq. (11). This part is analogous to that for SPEC_FT,
   // i.e., it has the form of the usual Lehmann representation. In a
   // typical calculation only a small contribution to the total
   // weight comes from this term. (This is the case both for T=0 and
   // T!=0 calculation.)
   if (LAST_ITERATION()) {
      dim1 = diagI1.getnr(); // override
      dimp = diagIp.getnr();
      for (size_t r1 = 0; r1 < dim1; r1++) {
	 const t_eigen E1 = diagI1.value(r1);
	 for (size_t rp = 0; rp < dimp; rp++) {
	    const t_eigen Ep = diagIp.value(rp);
	    DELTA d;
	    d.energy = E1 - Ep;
	    // Zft is Z_N, i.e. sum_i exp(-beta E_i) at the last NRG
	    // iteration.
	    d.weight = (spinfactor/Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp)
	      * exp(-E1 * scT) * (-sign); // (***)
	    cs->add(scale * d.energy, d.weight);
	 }
      }
   } else {
      // iii-term, Eq. (16), positive frequency excitations
      const size_t dimA = diagI1.getnr();
      for (size_t rl = dim1; rl < dimA; rl++) {
	 const t_eigen El = diagI1.value(rl);
	 for (size_t rk = 0; rk < dimp; rk++) {
	    const t_eigen Ek = diagIp.value(rk);
	    DELTA d;
	    d.energy = El - Ek;
	    my_assert(d.energy >= 0.0); // always positive!
	    weight_bucket sum;
	    for (size_t rkp = 0; rkp < dimp; rkp++) 
	      sum += op2II(rl, rkp) * rhoNIp(rkp, rk); // no sign here!
	    d.weight = spinfactor * CONJ_ME(op1II(rl, rk)) * t_weight(sum) * (-sign); // (***)
	    cs->add(scale * d.energy, d.weight);
	 }
      }
   } // if (last)
}

void SPEC_CFSgt::calc(const Eigen & diagIp, const Eigen & diagI1, 
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs,
		      t_factor spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ip, const Invar &I1)
{
   using namespace STAT;
   const Matrix & rhoNIp = rho[Ip]; // hand optimised out of the loops
   const Matrix & rhoNI1 = rho[I1];
   auto dimp = rhoNIp.size1();
   auto dim1 = rhoNI1.size1();
   // Convention: k-loops over retained states, l-loop over discarded
   // states.
   // i-term, Eq. (11). This part is analogous to that for SPEC_FT,
   // i.e., it has the form of the usual Lehmann representation. In a
   // typical calculation only a small contribution to the total
   // weight comes from this term. (This is the case both for T=0 and
   // T!=0 calculation.)
   if (LAST_ITERATION()) {
      dim1 = diagI1.getnr();
      dimp = diagIp.getnr();
      for (size_t r1 = 0; r1 < dim1; r1++) {
	 const t_eigen E1 = diagI1.value(r1);
	 for (size_t rp = 0; rp < dimp; rp++) {
	    const t_eigen Ep = diagIp.value(rp);
	    DELTA d;
	    d.energy = E1 - Ep;
	    // Zft is Z_N, i.e. sum_i exp(-beta E_i) at the last NRG
	    // iteration.
	    d.weight = (spinfactor/Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp)
	      * exp(-Ep * scT); // (***) removed (-sign)
	    cs->add(scale * d.energy, d.weight);
	 }
      }
   } else {
      // ii-term, Eq. (15), negative frequency excitations
      for (size_t rk = 0; rk < dim1; rk++) {
	 const t_eigen Ek = diagI1.value(rk);
	 const size_t dimB = diagIp.getnr();
	 for (size_t rl = dimp; rl < dimB; rl++) {
	    const t_eigen El = diagIp.value(rl);
	    DELTA d;
	    d.energy = Ek - El;
	    my_assert(d.energy <= 0.0); // always negative!
	    weight_bucket sum;
	    for (size_t rkp = 0; rkp < dim1; rkp++)
	      sum += CONJ_ME( op1II(rkp, rl) ) * rhoNI1(rkp, rk);
	    d.weight = spinfactor * t_weight(sum) * op2II(rk, rl); // (***) removed (-sign)
	    cs->add(scale * d.energy, d.weight);
	 }
      }
   } // if (last)
}
#endif

// Based on the implementation by Markus Greger.
#if defined(NRG_REAL) && defined(SPEC_CFS_OPTIMIZED)
void SPEC_CFSls::calc(const Eigen & diagIp, const Eigen & diagI1, 
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs, double spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ip, const Invar &I1)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   using namespace STAT;
   const Matrix & rhoNIp = rho[Ip]; // hand optimised out of the loops
   const Matrix & rhoNI1 = rho[I1];
   auto dimp = rhoNIp.size1();
   auto dim1 = rhoNI1.size1();
   // Convention: k-loops over retained states, l-loop over discarded
   // states.
   // i-term, Eq. (11). This part is analogous to that for SPEC_FT,
   // i.e. it has the form of the usual Lehmann representation. In a
   // typical calculation only a small contribution to the total
   // weight comes from this term.
   if (LAST_ITERATION()) {
      dim1 = diagI1.getnr();
      dimp = diagIp.getnr();
      for (size_t r1 = 0; r1 < dim1; r1++) {
	 const double E1 = diagI1.value(r1);
	 for (size_t rp = 0; rp < dimp; rp++) {
	    const double Ep = diagIp.value(rp);	    
	    double d_energy = E1 - Ep;
	    // Zft is Z_N, i.e. sum_i exp(-beta E_i) at the last NRG
	    // iteration.
	    double d_weight = (spinfactor/Zft) * op1II(r1, rp) * op2II(r1, rp)
	      * exp(-E1 * scT) * (-sign); // (***)
	    cs->add(scale * d_energy, d_weight);
	 }
      }
   } else { 
      // iii-term, Eq. (16), positive frequency excitations
      const size_t dimA = diagI1.getnr();
      DVEC::const_iterator energies_beginIp = begin(diagIp.value);
      DVEC::const_iterator energies_beginI1 = begin(diagI1.value);
      if (dimA && dimp) {
	 Matrix op2II_m_rho;
	 const matrix_range<const Matrix> op2II_TK(op2II, range(0, op2II.size1()), range(0, rhoNIp.size2()));
	 op2II_m_rho=Matrix(op2II_TK.size1(),rhoNIp.size2());
	 atlas::detail::symm(CblasRight, CblasUpper, 1.0, rhoNIp, op2II_TK,0.0, op2II_m_rho); // rhoNEW <- rhoNEW + factor T U
	 for (size_t rl = dim1; rl < dimA; rl++) {
	    const double El = *(energies_beginI1 + rl);
	    for (size_t rk = 0; rk < dimp; rk++) {
	       const double Ek = *(energies_beginIp + rk);
	       const double d_energy = El - Ek;
	       const double sum = op2II_m_rho(rl,rk);
	       const double d_weight = spinfactor * op1II(rl, rk) * sum * (-sign); // (***)
	       cs->add(scale * d_energy, d_weight);
	    }
	 }
      }
   } // if (last)
}

void SPEC_CFSgt::calc(const Eigen & diagIp, const Eigen & diagI1, 
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs, double spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ip, const Invar &I1)
{
   using namespace STAT;
   const Matrix & rhoNIp = rho[Ip]; // hand optimised out of the loops
   const Matrix & rhoNI1 = rho[I1];
   auto dimp = rhoNIp.size1();
   auto dim1 = rhoNI1.size1();
   // Convention: k-loops over retained states, l-loop over
   // discarded states.
   // i-term, Eq. (11). This part is analogous to that for SPEC_FT,
   // i.e. it has the form of the usual Lehmann representation. In a
   // typical calculation only a small contribution to the total
   // weight comes from this term.
   if (LAST_ITERATION()) {
      dim1 = diagI1.getnr();
      dimp = diagIp.getnr();
      for (size_t r1 = 0; r1 < dim1; r1++) {
	 const double E1 = diagI1.value(r1);
	 for (size_t rp = 0; rp < dimp; rp++) {
	    const double Ep = diagIp.value(rp);	    
	    double d_energy = E1 - Ep;
	    // Zft is Z_N, i.e. sum_i exp(-beta E_i) at the last NRG
	    // iteration.
	    double d_weight = (spinfactor/Zft) * op1II(r1, rp) * op2II(r1, rp)
	      * exp(-Ep * scT); // (***) removed (-sign)
	    cs->add(scale * d_energy, d_weight);
	 }
      }
   } else { 
      const size_t dimB = diagIp.getnr();
      DVEC::const_iterator energies_beginIp = begin(diagIp.value);
      DVEC::const_iterator energies_beginI1 = begin(diagI1.value);
      if (dim1 && dimB) {
	 const matrix_range<const Matrix> op1II_KT(op1II, range(0, rhoNI1.size1()), range(0, op1II.size2()));
	 Matrix op1II_m_rho(rhoNI1.size2(),op1II_KT.size2());
	 atlas::detail::symm(CblasLeft, CblasUpper, 1.0, rhoNI1, op1II_KT,0.0, op1II_m_rho); // rhoNEW <- rhoNEW + factor T U
	 for (size_t rk = 0; rk < dim1; rk++) { // ii-term, Eq. (15), negative frequency excitations
	    const double Ek = *(energies_beginI1+rk);	 
	    for (size_t rl = dimp; rl < dimB; rl++) {
	       const double El = *(energies_beginIp + rl);	    
	       const double d_energy = Ek - El;
	       const double sum = op1II_m_rho(rk,rl);
	       double d_weight = spinfactor * sum; // (***) removed (-sign)
	       d_weight *= op2II(rk, rl);
	       cs->add(scale * d_energy, d_weight);
	    }
	 }
      }
   } // if (last)
}
#endif
