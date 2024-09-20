#ifndef _algo_FDM_hpp_
#define _algo_FDM_hpp_

#include <complex>
#include "traits.hpp"
#include "algo.hpp"
#include "spectrum.hpp"

namespace NRG {

using namespace std::complex_literals;

// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FDMls : virtual public Algo<S> {
 private:
   inline static const std::string algoname = "FDMls";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   const bool save;
 protected:
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
    using Algo<S>::P;
   Algo_FDMls(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
     : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
   {
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto rho_op2 = prod_fit(rhoFDM.at(Ij), op2);
     const auto energies = [&diagIi, &diagIj](const auto i, const auto j) {
       return std::make_pair(diagIi.values.abs_G(i), diagIj.values.abs_G(j));
     };
     const auto term1 = [&energies, &op1, &op2, T = P.T.value(), wnf, this](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * op2(j, i) * (-sign) * exp(-Ej/T) * wnf);
     };
     const auto term2 = [&energies, &op1, &rho_op2, this](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * rho_op2(j, i) * (-sign));
     };
     const auto term3 = [&energies, &op1, &op2, T = P.T.value(), wnf, this](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * op2(j, i) * (-sign) * exp(-Ej/T) * wnf);
     };
     for (const auto i : diagIi.Drange()) // XXX: change order for sequential memory access!
       for (const auto j : diagIj.Drange())
         cb->add(term1(i,j), factor);
     for (const auto i : diagIi.Drange())
       for (const auto j : diagIj.Krange())
         cb->add(term2(i,j), factor);
     for (const auto i : diagIi.Krange())
       for (const auto j : diagIj.Drange())
         cb->add(term3(i,j), factor);
   }
   void end([[maybe_unused]] const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_FDMls() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FDMgt : virtual public Algo<S> {
 private:
   inline static const std::string algoname = "FDMgt";
   SpectrumRealFreq<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   const bool save;
 protected:
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
 public:
   using Algo<S>::P;
   Algo_FDMgt(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
     : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
   {
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto op2_rho = prod_fit(op2, rhoFDM.at(Ii));
     const auto energies = [&diagIi, &diagIj](const auto i, const auto j) {
       return std::make_pair(diagIi.values.abs_G(i), diagIj.values.abs_G(j));
     };
     const auto term1 = [&energies, &op1, &op2, T = P.T.value(), wnf](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * op2(j, i) * exp(-Ei/T) * wnf);
     };
     const auto term2 = [&energies, &op1, &op2, T = P.T.value(), wnf](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * op2(j, i) * exp(-Ei/T) * wnf);
     };
     const auto term3 = [&energies, &op1, &op2_rho](const auto i, const auto j) {
       const auto [Ei, Ej] = energies(i, j);
       return std::make_pair(Ej-Ei, conj_me(op1(j, i)) * op2_rho(j, i));
     };
     for (const auto i : diagIi.Drange())
       for (const auto j : diagIj.Drange())
         cb->add(term1(i,j), factor);
     for (const auto i : diagIi.Drange())
       for (const auto j : diagIj.Krange())
         cb->add(term2(i,j), factor);
     for (const auto i : diagIi.Krange())
       for (const auto j : diagIj.Drange())
         cb->add(term3(i,j), factor);
   }
   void end([[maybe_unused]] const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_FDMgt() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FDM : public Algo_FDMls<S>, public Algo_FDMgt<S> {
 private:
   inline static const std::string algoname2 = "FDM";
   SpectrumRealFreq<S> spec_tot;
 public:
   using Algo<S>::P;
   Algo_FDM(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), Algo_FDMls<S>(name, prefix, gt, P, false), Algo_FDMgt<S>(name, prefix, gt, P, false), spec_tot(name, algoname2, spec_fn(name, prefix, algoname2), P) {}
   void begin(const Step &step) override {
     Algo_FDMgt<S>::begin(step);
     Algo_FDMls<S>::begin(step);
   }
   void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho,
             const Stats<S> &stats) override
   {
     Algo_FDMgt<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_FDMls<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end([[maybe_unused]] const Step &step) override {
     spec_tot.mergeCFS(*Algo_FDMgt<S>::cb.get());
     spec_tot.mergeCFS(*Algo_FDMls<S>::cb.get());
     Algo_FDMgt<S>::cb.reset();
     Algo_FDMls<S>::cb.reset();
   }
   ~Algo_FDM() { spec_tot.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>, typename t_weight = weight_traits<S>>
class Algo_FDMmats : public Algo<S> {
 private:
   inline static const std::string algoname = "FDMmats";
   GFMatsubara<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara<S>;
   std::unique_ptr<CM> cm;
 public:
   using Algo<S>::P;
   Algo_FDMmats(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
   {
     const size_t cutoff = P.mats;
     const auto wnf      = stats.wnfactor[step.ndx()];
     const auto rho_op2  = prod_fit(rhoFDM.at(Ij), op2);
     const auto op2_rho  = prod_fit(op2, rhoFDM.at(Ii));
     const auto energies = [&diagIi, &diagIj](const auto i, const auto j) { 
       return std::make_pair(diagIi.values.abs_G(i), diagIj.values.abs_G(j));
     };
     const auto term1 = [&energies, &op1, &op2, T = P.T.value(), wnf, this](const auto i, const auto j, const auto n) -> t_weight {
       const auto [Ei, Ej] = energies(i, j);
       const auto energy = Ej - Ei;
       const auto weightA = conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei/T); // a[ij] b[ji] exp(-beta e[i])
       const auto weightB = conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej/T); // a[ij] b[ji] sign exp(-beta e[j])
       if (gt == gf_type::fermionic || n>0 || abs(energy) > WEIGHT_TOL)
          return (weightA + weightB) / (ww(n, gt, T)*1i - energy);
       else // bosonic w=0 && Ei=Ej case
          return -weightA/T;
     };
     const auto term2 = [&energies, &op1, &op2, &rho_op2, T = P.T.value(), wnf, this](const auto i, const auto j, const auto n) {
       const auto [Ei, Ej] = energies(i, j);
         const auto energy = Ej - Ei;
         const auto weightA = conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei/T);
         const auto weightB = conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
         return (weightA + weightB) / (ww(n, gt, T)*1i - energy);
     };
     const auto term3 = [&energies, &op1, &op2, &op2_rho, T = P.T.value(), wnf, this](const auto i, const auto j, const auto n) {
       const auto [Ei, Ej] = energies(i, j);
       const auto energy = Ej - Ei;
       const auto weightA = conj_me(op1(j, i)) * op2_rho(j, i);
       const auto weightB = (-sign) * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ej/T);
       return (weightA + weightB) / (ww(n, gt, T)*1i - energy);
     };
     for (const auto i : diagIi.Drange())
       for (const auto j : diagIj.Drange())
         #pragma omp parallel for schedule(static)
         for (size_t n = 0; n < cutoff; n++)
           cm->add(n, term1(i,j,n) * factor);
     for (const auto i : diagIi.Drange())
       for (const auto j : diagIj.Krange())
         #pragma omp parallel for schedule(static)
         for (size_t n = 0; n < cutoff; n++)
           cm->add(n, term2(i,j,n) * factor);
     for (const auto i : diagIi.Krange())
       for (const auto j : diagIj.Drange())
         #pragma omp parallel for schedule(static)
         for (size_t n = 0; n < cutoff; n++)
           cm->add(n, term3(i,j,n) * factor);
   }
   void end([[maybe_unused]] const Step &step) override {
     gf.merge(*cm.get());
     cm.reset();
   }
   ~Algo_FDMmats() { gf.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

} // namespace

#endif
