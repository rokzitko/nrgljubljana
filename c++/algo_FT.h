#ifndef _algo_FT_h_
#define _algo_FT_h_

template<typename S>
class Algo_FT : public Algo<S> {
 private:
   inline static const std::string algoname = "FT";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo<S>::P;
   Algo_FT(const std::string &name, const std::string &prefix, const gf_type &gt, const Params &P) :
     Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname), P), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   // The first matrix element is conjugated! This is <rp|OP1^dag|r1> <r1|OP2|rp> (wp - s*w1)/(z+Ep-E1)
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             const t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     for (const auto r1: diagI1.kept()) {
       const auto E1 = diagI1.value_zero(r1);
       for (const auto rp: diagIp.kept()) {
         const auto Ep = diagIp.value_zero(rp);
         const auto weight = (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * ((-sign) * exp(-E1 * step.scT()) + exp(-Ep * step.scT()));
         const auto energy = E1 - Ep;
         cb->add(step.scale() * energy, weight);
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeNN2(*cb.get(), step);
     cb.reset();
   }
   ~Algo_FT() { spec.save(); }
};

template<typename S>
class Algo_FTmats : public Algo<S> {
 private:
   inline static const std::string algoname = "FTmats";
   GFMatsubara<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara<S>;
   std::unique_ptr<CM> cm;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo<S>::P;
   Algo_FTmats(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     const size_t cutoff = P.mats;
     for (const auto r1: diagI1.kept()) {
       const auto E1 = diagI1.value_zero(r1);
       for (const auto rp: diagIp.kept()) {
         const auto Ep = diagIp.value_zero(rp);
         const auto weight = (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * ((-sign) * exp(-E1 * step.scT()) + exp(-Ep * step.scT()));
         const auto energy = E1 - Ep;
#pragma omp parallel for schedule(static)
         for (size_t n = 1; n < cutoff; n++) cm->add(n, weight / (cmpl(0, ww(n, gt, P.T)) - step.scale() * energy));
         if (abs(energy) > WEIGHT_TOL || gt == gf_type::fermionic)
           cm->add(size_t(0), weight / (cmpl(0, ww(0, gt, P.T)) - step.scale() * energy));
         else // bosonic w=0 && E1=Ep case
           cm->add(size_t(0), (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * (-exp(-E1 * step.scT()) / P.T));
       }
     }
   }
   void end(const Step &step) override {
     gf.merge(*cm.get());
     cm.reset();
   }
   ~Algo_FTmats() { gf.save(); }
};

// Calculation of the temperature-dependent linear conductrance G(T) using the linear response theory &
// impurity-level spectral density.  See Yoshida, Seridonio, Oliveira, arxiv:0906.4289, Eq. (8).
template<typename S, int n>
class Algo_GT : public Algo<S> {
 private:
   inline static const std::string algoname = n == 0 ? "GT" : (n == 1 ? "I1T" : "I2T");
   TempDependence<S> td;
   using CT = ChainTempDependence<S>;
   std::unique_ptr<CT> ct;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo<S>::P;
   Algo_GT(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) : 
     Algo<S>(P), td(name, algoname, spec_fn(name, prefix, algoname), P) {
     my_assert(gt == gf_type::fermionic);
     static_assert(n ==0 || n == 1 || n == 2);
   }
   void begin(const Step &) override { ct = std::make_unique<CT>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override 
   {
     const double temperature = P.gtp * step.scale(); // in absolute units!
     const auto beta          = 1.0 / temperature;
     typename traits<S>::t_weight value{};
     for (const auto r1: diagI1.kept()) {
       const auto E1 = diagI1.value_zero(r1);
       for (const auto rp: diagIp.kept()) {
         const auto Ep = diagIp.value_zero(rp);
         // Note that Zgt needs to be calculated with the same 'temperature' parameter that we use for the
         // exponential functions in the following equation.
         value += beta * (factor / stats.Zgt) * conj_me(op1(r1, rp)) * op2(r1, rp)
           / (exp(+E1 * step.scale() * beta) + exp(+Ep * step.scale() * beta)) * pow((E1 - Ep) * step.scale(), n);
       } // loop over r1
     }   // loop over rp
     ct->add(temperature, value);
   }
   void end(const Step &) override {
     td.merge(*ct.get());
   }
   ~Algo_GT() { td.save(); }
};

// weight=(exp(-beta Em)-exp(-beta En))/(beta En-beta Em). NOTE: arguments En, Em are order omega_N, while beta is
// order 1/omega_N, thus the combinations betaEn and betaEm are order 1. Also En>0, Em>0, since these are excitation
// energies !
inline auto chit_weight(const double En, const double Em, const double beta) {
  const auto betaEn = beta * En;
  const auto betaEm = beta * Em;
  const auto x      = betaEn - betaEm;
  if (abs(x) > WEIGHT_TOL) {
    // If one of {betaEm,betaEn} is small, one of exp() will have a value around 1, the other around 0, thus the
    // overall result will be approximately +-1/x.
    return (exp(-betaEm) - exp(-betaEn)) / x;
  } else {
    // Special case for Em~En. In this case, we are integrating a constant over tau\in{0,\beta}, and dividing this by
    // beta we get 1. What remains is the Boltzmann weight exp(-betaEm).
    return exp(-betaEm);
  }
}

// Calculation of the temperature-dependent susceptibility chi_AB(T) using the linear response theory and the matrix
// elements of global operators. Binning needs to be turned off. Note that Zchit needs to be calculated with the same
// 'temperature' parameter that we use for the exponential functions in the following equation. The output is
// chi/beta = k_B T chi, as we prefer.
template<typename S>
class Algo_CHIT : public Algo<S> {
 private:
   inline static const std::string algoname = "CHIT";
   TempDependence<S> td;
   using CT = ChainTempDependence<S>;
   std::unique_ptr<CT> ct;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo<S>::P;
   Algo_CHIT(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) : 
     Algo<S>(P), td(name, algoname, spec_fn(name, prefix, algoname), P) {
     my_assert(gt == gf_type::bosonic);
   }
   void begin(const Step &) override { ct = std::make_unique<CT>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     const double temperature = P.chitp * step.scale(); // in absolute units!
     const auto beta          = 1.0 / temperature;
     typename traits<S>::t_weight value{};
     for (const auto r1: diagI1.kept()) {
       for (const auto rp: diagIp.kept()) {
         const auto E1 = diagI1.value_zero(r1);
         const auto Ep = diagIp.value_zero(rp);
         value += (factor/stats.Zchit) * chit_weight(step.scale()*E1, step.scale()*Ep, beta) * conj_me(op1(r1,rp)) * op2(r1,rp);
       }
     }
     ct->add(temperature, value);
   }
   void end(const Step &) override {
     td.merge(*ct.get());
   }
   ~Algo_CHIT() { td.save(); }
};

#endif
