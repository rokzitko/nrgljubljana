#ifndef _algo_DMNRG_hpp_
#define _algo_DMNRG_hpp_

#include <complex>
using namespace std::complex_literals;
#include "algo.hpp"
#include "spectrum.hpp"

namespace NRG {

// OPTIMIZATION NOTE: the inner loop should involve the last index.

template<typename S>
class Algo_DMNRG : public Algo<S> {
 private:
   inline static const std::string algoname = "DMNRG";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using t_weight = typename traits<S>::t_weight;
   using Algo<S>::P;
   Algo_DMNRG(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname), P), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     const double Emin = P.ZBW ? 0 : P.getEmin();
     const double Emax = P.ZBW ? std::numeric_limits<double>::max() : P.getEmax();
     const Matrix &rhoNIp = rho.at(Ip);
     const Matrix &rhoNI1 = rho.at(I1);
     for (const auto rm: diagIp.kept()) {
       for (const auto rj: diagI1.kept()) {
         const auto Em = diagIp.value_zero(rm);
         const auto Ej = diagI1.value_zero(rj);
         const auto energy = Ej - Em;
         const auto absE = abs(energy);
         if (absE < Emin || absE > Emax) // does not contribute
           continue;
         t_weight sumA{};
         for (const auto ri: diagIp.kept()) sumA += op2(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
         const auto weightA = sumA * conj_me(op1(rj, rm));
         t_weight sumB{};
         for (const auto ri: diagI1.kept()) sumB += conj_me(op1(ri, rm)) * rhoNI1(rj, ri); // non-optimal
         const auto weightB = sumB * op2(rj, rm);
         const auto weight  = factor * (weightA + (-sign) * weightB);
         cb->add(step.scale() * energy, weight);
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

template<typename S>
class Algo_DMNRGmats : public Algo<S> {
 private:
   inline static const std::string algoname = "DMNRGmats";
   GFMatsubara<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara<S>;
   std::unique_ptr<CM> cm;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using t_weight = typename traits<S>::t_weight;
   using Algo<S>::P;
   Algo_DMNRGmats(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     const Matrix &rhoNIp = rho.at(Ip);
     const Matrix &rhoNI1 = rho.at(I1);
     for (const auto rm: diagIp.kept()) {
       for (const auto rj: diagI1.kept()) {
         const auto Em = diagIp.value_zero(rm);
         const auto Ej = diagI1.value_zero(rj);
         const auto energy = Ej - Em;
         t_weight sumA{};
         for (const auto ri: diagIp.kept()) sumA += op2(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
         const auto weightA = sumA * conj_me(op1(rj, rm));
         t_weight sumB{};
         for (const auto ri: diagI1.kept()) sumB += conj_me(op1(ri, rm)) * rhoNI1(rj, ri); // non-optimal
         const auto weightB = sumB * op2(rj, rm);
         const auto weight  = factor * (weightA + (-sign) * weightB);
         for (size_t n = 1; n < P.mats; n++) cm->add(n, weight / (ww(n, gt, P.T)*1i - step.scale() * energy));
         if (abs(energy) > WEIGHT_TOL || gt == gf_type::fermionic)
           cm->add(size_t(0), weight / (ww(0, gt, P.T)*1i - step.scale() * energy));
         else // bosonic w=0 && E1=Ep case
           cm->add(size_t(0), factor * (-weightA / t_weight(P.T)));
       }
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
