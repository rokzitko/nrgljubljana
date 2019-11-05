class SPEC_FT : public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "FT"; }
   string merge() { return "NN2"; }
};

const double WEIGHT_TOL = 1e-8; // where to switch to l'Hospital rule form

// The first matrix element is conjugated!
// This is <rp|OP1^dag|r1> <r1|OP2|rp> (wp - s*w1)/(z+Ep-E1)
// s=1 for bosons, s=-1 for fermions
// See Eq.(9) in Peters, Pruschke, Anders, PRB (74) 245114 (2006)
void SPEC_FT::calc(const Eigen & diagIp, const Eigen & diagI1, 
		   const Matrix & op1II, const Matrix & op2II,
		   const BaseSpectrum &bs, t_factor spinfactor,
		   ChainSpectrum *cs, 	
		   const Invar &Ip, const Invar &I1)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   using namespace STAT;
   for (size_t r1 = 0; r1 < diagI1.getnr(); r1++) {
      const t_eigen E1 = diagI1.value(r1);
      for (size_t rp = 0; rp < diagIp.getnr(); rp++) {
	 const t_eigen Ep = diagIp.value(rp);
	 DELTA d;
	 d.weight = (spinfactor/Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp)
	    * ( (-sign) * exp(-E1 * scT) + exp(-Ep * scT) );
	 d.energy = E1 - Ep;
	 cs->add(scale * d.energy, d.weight);
      }
   }
}

class SPEC_FTmats : public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &bs) { return new ChainSpectrumMatsubara(bs.mt); }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "FTmats"; }
};

void SPEC_FTmats::calc(const Eigen & diagIp, const Eigen & diagI1, 
		       const Matrix & op1II, const Matrix & op2II,
		       const BaseSpectrum &bs, t_factor spinfactor,
		       ChainSpectrum *cs,		
		       const Invar &Ip, const Invar &I1)
{
   const size_t cutoff = P::mats;
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   auto *csm = dynamic_cast<ChainSpectrumMatsubara*>(cs);
   using namespace STAT;
   for (size_t r1 = 0; r1 < diagI1.getnr(); r1++) {
      const t_eigen E1 = diagI1.value(r1);
      for (size_t rp = 0; rp < diagIp.getnr(); rp++) {
	 const t_eigen Ep = diagIp.value(rp);
	 DELTA d;
	 d.weight = (spinfactor/Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp)
	   * ( (-sign) * exp(-E1 * scT) + exp(-Ep * scT) ); // sign!
	 d.energy = E1 - Ep;
#pragma omp parallel for schedule(static)
	 for (size_t n = 1; n < cutoff; n++)
	    csm->add(n, d.weight/(cmpl(0,w(n,bs.mt))-scale*d.energy));
	 if (abs(d.energy) > WEIGHT_TOL || bs.mt == matstype::fermionic)
	    csm->add(size_t(0), d.weight/(cmpl(0,w(0,bs.mt))-scale*d.energy));
	 else // bosonic w=0 && E1=Ep case
	    csm->add(size_t(0), (spinfactor/Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp) * (-exp(-E1*scT)/P::T));
      }
   }
}

class SPEC_GT_generic : public SPEC
{
protected:
   int power;
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumTemp; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "ERROR"; }
};

class SPEC_GT : public SPEC_GT_generic
{
public:
   SPEC_GT() { power = 0; }
   string name() { return "GT"; }
};

class SPEC_I1T : public SPEC_GT_generic
{
public:
   SPEC_I1T() { power = 1; }
   string name() { return "I1T"; }
};

class SPEC_I2T : public SPEC_GT_generic
{
public:
   SPEC_I2T() { power = 2; }
   string name() { return "I2T"; }
};

// Calculation of the temperature-dependent linear conductrance G(T) using
// the linear response theory & impurity-level spectral density. Note that
// we are not calculating a spectral function (i.e. a collection of delta
// peaks), but rather a tabulated G(T), so binning needs to be turned off.
// See Yoshida, Seridonio, Oliveira, arxiv:0906.4289, Eq. (8).
void SPEC_GT_generic::calc(const Eigen & diagIp, const Eigen & diagI1, 
			   const Matrix & op1II, const Matrix & op2II,
			   const BaseSpectrum &bs, t_factor spinfactor,
			   ChainSpectrum * cs,
			   const Invar &Ip, const Invar &I1)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   my_assert(sign == S_FERMIONIC); // restricted implementation
   const double temperature = P::gtp * STAT::scale; // in absolute units!
   const double beta = 1.0/temperature; 
   weight_bucket value;
   for (size_t r1 = 0; r1 < diagI1.getnr(); r1++) {
      const t_eigen E1 = diagI1.value(r1);
      for (size_t rp = 0; rp < diagIp.getnr(); rp++) {
	 const t_eigen Ep = diagIp.value(rp);
 	 // Note that Zgt needs to be calculated with the same
 	 // 'temperature' parameter that we use for the exponential
 	 // functions in the following equation.
	 value += beta * (spinfactor/STAT::Zgt) * CONJ_ME( op1II(r1, rp) ) * op2II(r1, rp) 
	   / (exp(+E1 * STAT::scale * beta) + exp(+Ep * STAT::scale * beta)) 
	     * pow((E1-Ep)*STAT::scale, power);
      } // loop over r1
   } // loop over rp
   cs->add(temperature, value);
}

// weight=(exp(-beta Em)-exp(-beta En))/(beta En-beta Em).
// NOTE: arguments En, Em are order omega_N, while beta is order
// 1/omega_N, thus the combinations betaEn and betaEm are order 1.
// Also En>0, Em>0, since these are excitation energies !
inline t_weight chit_weight(double En, double Em, double beta)
{
   const double betaEn = beta * En;
   const double betaEm = beta * Em;
   const double x = betaEn-betaEm;
   if (abs(x) > WEIGHT_TOL) {
      // If one of {betaEm,betaEn} is small, one of exp() will have a
      // value around 1, the other around 0, thus the overall result
      // will be approximately +-1/x.
      return (exp(-betaEm)-exp(-betaEn))/x;
   } else {
      // Special case for Em~En. In this case, we are integrating a
      // constant over tau\in{0,\beta}, and dividing this by beta we
      // get 1. What remains is the Boltzmann weight exp(-betaEm).
      return exp(-betaEm);
   }
}
 
class SPEC_CHIT : public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumTemp; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "CHIT"; }
};

// Calculation of the temperature-dependent susceptibility chi_AB(T)
// using the linear response theory and the matrix elements of global
// operators. Binning needs to be turned off.
// Note that Zchit needs to be calculated with the same 'temperature'
// parameter that we use for the exponential functions in the
// following equation. 
// NOTE: The output is chi/beta = k_B T chi, as we prefer.
void SPEC_CHIT::calc(const Eigen & diagIp, const Eigen & diagI1, 
		     const Matrix & op1II, const Matrix & op2II,
		     const BaseSpectrum &bs, t_factor spinfactor,
		     ChainSpectrum * cs,
		     const Invar &Ip, const Invar &I1)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   my_assert(sign == S_BOSONIC); // restricted implementation
   const size_t dimp = diagIp.getnr();
   const size_t dim1 = diagI1.getnr();
   my_assert(dim1 == dimp);
   const double temperature = P::chitp * STAT::scale; // in absolute units!
   const double beta = 1.0/temperature; 
   weight_bucket w;
   for (size_t r1 = 0; r1 < dim1; r1++) {
      const t_eigen E1 = diagI1.value(r1);
      for (size_t rp = 0; rp < dimp; rp++) {
	 const t_eigen Ep = diagIp.value(rp);
	 const t_weight weight = chit_weight(STAT::scale * E1, STAT::scale * Ep, beta);
	 w += op1II(r1, rp) * op2II(rp, r1) * weight;
      }
   }
   cs->add(temperature, spinfactor/STAT::Zchit * t_weight(w));
}
