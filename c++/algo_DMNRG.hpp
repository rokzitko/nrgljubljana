#ifndef _algo_DMNRG_hpp_
#define _algo_DMNRG_hpp_

#include <complex>
#include <tuple>
#include "traits.hpp"
#include "algo.hpp"
#include "spectrum.hpp"

namespace NRG {

using namespace std::complex_literals;

// OPTIMIZATION NOTE: the inner loop should involve the last index.

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>, typename t_weight = weight_traits<S>>
class Algo_DMNRG : public Algo<S> {
 private:
   inline static const std::string algoname = "DMNRG";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Algo<S>::P;
   Algo_DMNRG(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname), P), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     auto weights = [Emin = step.scale() * P.getEmin(), Emax = step.scale() * P.getEmax(), &rhoNIp = rho.at(Ip), &rhoNI1 = rho.at(I1), &diagIp, &diagI1, &op1, &op2](const auto rm, const auto rj) { 
       const auto Em = diagIp.values.abs_zero(rm);
       const auto Ej = diagI1.values.abs_zero(rj);
       const auto energy = Ej-Em;
       if (abs(energy) < Emin || abs(energy) > Emax) return std::make_tuple(energy, t_weight{}, t_weight{}); // does not contribute
       t_weight sumA{};
       for (const auto ri: diagIp.kept()) sumA += op2(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
       const auto weightA = sumA * conj_me(op1(rj, rm));
       t_weight sumB{};
       for (const auto ri: diagI1.kept()) sumB += conj_me(op1(ri, rm)) * rhoNI1(rj, ri); // non-optimal
       const auto weightB = sumB * op2(rj, rm);
       return std::make_tuple(energy, weightA, weightB);
     };
     auto term = [&weights, this](const auto rm, const auto rj) {
       const auto [energy, weightA, weightB] = weights(rm, rj);
       return std::make_pair(energy, weightA + (-sign) * weightB);
     };
     for (const auto rm: diagIp.kept()) {
       for (const auto rj: diagI1.kept()) {
         const auto [energy, weight] = term(rm, rj);
         cb->add(energy, factor * weight);
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeNN2(*cb.get(), step);
     cb.reset();
   }
   ~Algo_DMNRG() { spec.save(); }
   std::string rho_type() override { return "rho"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>, typename t_weight = weight_traits<S>>
class Algo_DMNRGmats : public Algo<S> {
 private:
   inline static const std::string algoname = "DMNRGmats";
   GFMatsubara<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara<S>;
   std::unique_ptr<CM> cm;
 public:
   using Algo<S>::P;
   Algo_DMNRGmats(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
      auto weights = [&rhoNIp = rho.at(Ip), &rhoNI1 = rho.at(I1), &diagIp, &diagI1, &op1, &op2, this](const auto rm, const auto rj) { 
         const auto Em = diagIp.values.abs_zero(rm);
         const auto Ej = diagI1.values.abs_zero(rj);
         t_weight sumA{};
         for (const auto ri: diagIp.kept()) sumA += op2(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
         const auto weightA = sumA * conj_me(op1(rj, rm));
         t_weight sumB{};
         for (const auto ri: diagI1.kept()) sumB += conj_me(op1(ri, rm)) * rhoNI1(rj, ri); // non-optimal
         const auto weightB = sumB * op2(rj, rm);
         return std::make_tuple(Ej-Em, weightA, weightB);
     };
     auto term = [&weights, this](const auto rm, const auto rj, const auto n) {
       const auto [energy, weightA, weightB] = weights(rm, rj);
       if (gt == gf_type::fermionic || n>0 || abs(energy) > WEIGHT_TOL) // [[likely]]
         return (weightA + (-sign) * weightB) / (ww(n, gt, P.T)*1i - energy);
       else // bosonic w=0 && Em=Ej case
         return -weightA / t_weight(P.T);
     };
     for (const auto rm: diagIp.kept())
       for (const auto rj: diagI1.kept())
         for (size_t n = 0; n < P.mats; n++) {
           const auto weight = term(rm, rj, n);
           cm->add(n, factor * weight); 
         }
   }
   void end(const Step &step) override {
          gf.merge(*cm.get());
          cm.reset();
   }
   ~Algo_DMNRGmats() { gf.save(); }
   std::string rho_type() override { return "rho"; }
};

} // namespace

#endif
