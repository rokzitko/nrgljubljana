#ifndef _algo_FT_hpp_
#define _algo_FT_hpp_

#include <complex>
#include "traits.hpp"
#include "algo.hpp"
#include "spectrum.hpp"

// The first matrix element is conjugated! This is <rp|OP1^dag|r1> <r1|OP2|rp> (wp - s*w1)/(z+Ep-E1)

namespace NRG {

using namespace std::complex_literals;

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FT : public Algo<S> {
 private:
   inline static const std::string algoname = "FT";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Algo<S>::P;
   Algo_FT(const std::string &name, const std::string &prefix, const gf_type &gt, const Params &P) :
     Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname), P), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             const t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     auto stat_factor = [beta = step.scT(), Z = stats.Zft, this](const auto E1, const auto Ep) {
       return ((-sign) * exp(-beta*E1) + exp(-beta*Ep))/Z;
     };
     auto term = [&diagI1, &diagIp, &op1, &op2, &stat_factor](const auto r1, const auto rp) {
       const auto E1 = diagI1.values.rel_zero(r1);
       const auto Ep = diagIp.values.rel_zero(rp);
       return std::make_pair(E1 - Ep, conj_me(op1(r1, rp)) * op2(r1, rp) * stat_factor(E1,Ep));
     };
     for (const auto r1: diagI1.kept()) {
       for (const auto rp: diagIp.kept()) {
         const auto [energy, weight] = term(r1, rp);
         cb->add(step.scale() * energy, factor * weight);
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeNN2(*cb.get(), step);
     cb.reset();
   }
   ~Algo_FT() { spec.save(); }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>, typename t_weight = weight_traits<S>>
class Algo_FTmats : public Algo<S> {
 private:
   inline static const std::string algoname = "FTmats";
   GFMatsubara<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara<S>;
   std::unique_ptr<CM> cm;
 public:
   using Algo<S>::P;
   Algo_FTmats(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     auto stat_factor = [beta = step.scT(), scale = step.scale(), Z = stats.Zft, T = P.T, this](const auto E1, const auto Ep, const auto n) -> t_weight {
       const auto energy = E1-Ep;
       if (gt == gf_type::fermionic || n>0 || abs(energy) > WEIGHT_TOL) // [[likely]]
         return ((-sign) * exp(-beta*E1) + exp(-beta*Ep)) / (Z * (ww(n, gt, T)*1i - scale*energy));
       else // bosonic w=0 && E1=Ep case
         return -exp(-beta*E1) / (Z * T);
     };
     auto term = [&diagI1, &diagIp, &op1, &op2, &stat_factor](const auto r1, const auto rp, const auto n) {
       const auto E1 = diagI1.values.rel_zero(r1);
       const auto Ep = diagIp.values.rel_zero(rp);
       return conj_me(op1(r1, rp)) * op2(r1, rp) * stat_factor(E1,Ep,n);
     };
     const size_t cutoff = P.mats;
     for (const auto r1: diagI1.kept()) {
       for (const auto rp: diagIp.kept()) {
#pragma omp parallel for schedule(static)
         for (size_t n = 0; n < cutoff; n++) {
           const auto weight = term(r1, rp, n);
           cm->add(n, factor * weight);
         }
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
template<scalar S, int n, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_GT : public Algo<S> {
 private:
   inline static const std::string algoname = n == 0 ? "GT" : (n == 1 ? "I1T" : "I2T");
   TempDependence<S> td;
   using CT = ChainTempDependence<S>;
   std::unique_ptr<CT> ct;
 public:
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
     const double temperature = P.gtp * step.scale(); // in absolute units! stats.Zgt is evaluated for this temperature.
     auto stat_factor = [beta = 1.0/P.gtp, scale = step.scale(), Z = stats.Zgt](const auto E1, const auto Ep) {
       return (beta/scale) / (exp(+beta*E1) + exp(+beta*Ep)) * pow((E1 - Ep) * scale, n)/Z; // n is template parameter
     };
     auto term = [&diagI1, &diagIp, &op1, &op2, &stat_factor](const auto r1, const auto rp) {
       const auto E1 = diagI1.values.rel_zero(r1);
       const auto Ep = diagIp.values.rel_zero(rp);
       return conj_me(op1(r1, rp)) * op2(r1, rp) * stat_factor(E1,Ep);
     };
     weight_traits<S> value{};
     for (const auto r1: diagI1.kept())
       for (const auto rp: diagIp.kept())
         value += term(r1, rp); 
     ct->add(temperature, factor * value);
   }
   void end(const Step &) override {
     td.merge(*ct.get());
   }
   ~Algo_GT() { td.save(); }
};

// Calculation of the temperature-dependent susceptibility chi_AB(T) using the linear response theory and the matrix
// elements of global operators. Binning needs to be turned off. Note that Zchit needs to be calculated with the same
// 'temperature' parameter that we use for the exponential functions in the following equation. The output is
// chi/beta = k_B T chi, as we prefer.
template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_CHIT : public Algo<S> {
 private:
   inline static const std::string algoname = "CHIT";
   TempDependence<S> td;
   using CT = ChainTempDependence<S>;
   std::unique_ptr<CT> ct;
 public:
   using Algo<S>::P;
   Algo_CHIT(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) : 
     Algo<S>(P), td(name, algoname, spec_fn(name, prefix, algoname), P) {
     my_assert(gt == gf_type::bosonic);
   }
   void begin(const Step &) override { ct = std::make_unique<CT>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &, const Invar &, const DensMatElements<S> &, const Stats<S> &stats) override
   {
     const double temperature = P.chitp * step.scale(); // in absolute units! stats.Zchit is evaluated for this temperature.
     auto stat_factor = [temperature, beta = 1.0/P.chitp, scale = step.scale(), Z = stats.Zchit](const auto E1, const auto Ep) {
       return chit_weight(scale*E1, scale*Ep, 1.0/temperature)/Z;
     };
     auto term = [&diagI1, &diagIp, &op1, &op2, &stat_factor](const auto r1, const auto rp) {
       const auto E1 = diagI1.values.rel_zero(r1);
       const auto Ep = diagIp.values.rel_zero(rp);
       return conj_me(op1(r1, rp)) * op2(r1, rp) * stat_factor(E1,Ep);
     };
     weight_traits<S> value{};
     for (const auto r1: diagI1.kept())
       for (const auto rp: diagIp.kept())
         value += term(r1, rp);
     ct->add(temperature, factor * value);
   }
   void end(const Step &) override {
     td.merge(*ct.get());
   }
   ~Algo_CHIT() { td.save(); }
};

} // namespace

#endif
