#ifndef _algo_CFS_hpp_
#define _algo_CFS_hpp_

#include "traits.hpp"
#include "algo.hpp"
#include "spectrum.hpp"

namespace NRG {

// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).
// Based on the implementation by Markus Greger.

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_CFSls : virtual public Algo<S> {
 private:
   inline static const std::string algoname = "CFSls";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   const bool save; // if true, save spectral function to a file in destructor
 protected: 
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Algo<S>::P;
   Algo_CFSls(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
     : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, [[maybe_unused]] const Invar &Ip, [[maybe_unused]] const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     // Convention: k-loops over retained states, l-loop over discarded states.
     if (step.last()) {
       //  i-term, Eq. (11).
       const auto term1 = [&diagI1, &diagIp, Z=stats.Zft, &op1, &op2, T = P.T.value(), this](const auto r1, const auto rp) {
         const auto E1     = diagI1.values.abs_zero(r1);
         const auto Ep     = diagIp.values.abs_zero(rp);
         const auto weight = conj_me(op1(r1, rp)) * op2(r1, rp) * exp(-E1/T) * (-sign)/Z;
         return std::make_pair(E1-Ep, weight);
       };
       for (const auto r1: diagI1.kept())
         for (const auto rp: diagIp.kept())
           cb->add(term1(r1, rp), factor);
     } else {
       // iii-term, Eq. (16), positive frequency excitations
       const auto op2_rho = prod_fit_left(op2, rho.at(Ip));
       const auto term3 = [&diagI1, &diagIp, &op1, &op2_rho, this](const auto rl, const auto rk) {
         const auto El     = diagI1.values.abs_zero(rl);
         const auto Ek     = diagIp.values.abs_zero(rk);
         const auto weight = conj_me(op1(rl, rk)) * op2_rho(rl, rk) * (-sign);
         return std::make_pair(El-Ek, weight);
       };
       for (const auto rl: diagI1.discarded())
         for (const auto rk: diagIp.kept())
           cb->add(term3(rl, rk), factor);
     }
   }
   void end([[maybe_unused]] const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_CFSls() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_CFSgt : virtual public Algo<S> {
 private:
   inline static const std::string algoname = "CFSgt";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   const bool save;
 protected:
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Algo<S>::P;
   Algo_CFSgt(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
     : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     // Convention: k-loops over retained states, l-loop over discarded states.
     if (step.last()) {
        // i-term, Eq. (11).
        const auto term1 = [&diagI1, &diagIp, Z=stats.Zft, &op1, &op2, T = P.T.value()](const auto r1, const auto rp) {
         const auto E1     = diagI1.values.abs_zero(r1);
         const auto Ep     = diagIp.values.abs_zero(rp);
         const auto weight = conj_me(op1(r1, rp)) * op2(r1, rp) * exp(-Ep/T)/Z;
         return std::make_pair(E1-Ep, weight);
       };
       for (const auto r1: diagI1.kept())
         for (const auto rp: diagIp.kept())
           cb->add(term1(r1, rp), factor);
     } else {
       // ii-term, Eq. (15), negative frequency excitations
       const auto op1_rho = prod_adj_fit_left(op1, rho.at(I1));
       const auto term2 = [&diagI1, &diagIp, &op1_rho, &op2](const auto rk, const auto rl) {
         const auto Ek     = diagI1.values.abs_zero(rk);
         const auto El     = diagIp.values.abs_zero(rl);
         const auto weight = op1_rho(rl, rk) * op2(rk, rl);
         return std::make_pair(Ek-El, weight);
       };
       for (const auto rk: diagI1.kept())
         for (const auto rl: diagIp.discarded())
           cb->add(term2(rk, rl), factor);
     }
   }
   void end([[maybe_unused]] const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_CFSgt() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_CFS : public Algo_CFSls<S>, public Algo_CFSgt<S> {
 private:
   inline static const std::string algoname2 = "CFS";
   SpectrumRealFreq<S> spec_tot;
 public:
   using Algo<S>::P;
   Algo_CFS(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), Algo_CFSls<S>(name, prefix, gt, P, false), Algo_CFSgt<S>(name, prefix, gt, P, false), spec_tot(name, algoname2, spec_fn(name, prefix, algoname2), P) {}
   void begin(const Step &step) override {
     Algo_CFSgt<S>::begin(step);
     Algo_CFSls<S>::begin(step);
   }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho, const Stats<S> &stats) override
   {
     Algo_CFSgt<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_CFSls<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end([[maybe_unused]] const Step &step) override {
     spec_tot.mergeCFS(*Algo_CFSgt<S>::cb.get());
     spec_tot.mergeCFS(*Algo_CFSls<S>::cb.get());
     Algo_CFSgt<S>::cb.reset();
     Algo_CFSls<S>::cb.reset();
   }
   ~Algo_CFS() { spec_tot.save(); }
   std::string rho_type() override { return "rho"; }
};

} // namespace

#endif
