#ifndef _algo_FDM_hpp_
#define _algo_FDM_hpp_

#include <complex>
#include <vector>
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
    const bool should_save;
 protected:
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
  public:
     using Algo<S>::P;
    Algo_FDMls(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
      : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), should_save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, [[maybe_unused]] const Invar &Ii, [[maybe_unused]] const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
    {
      const auto wnf     = stats.wnfactor[step.ndx()];
      const auto T       = P.T.value();
      const auto rho_op2 = prod_fit(rhoFDM.at(Ij), op2);
      const auto absGi   = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGj = absGj | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;
      const auto add_exp_weight = [this, &absGi, &absGj, &boltzGj, &op1, &op2, factor, weight_scale = (-sign) * wnf](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2(j, i) * weight_scale * boltzGj[j]), factor);
      };
      const auto add_rho_weight = [this, &absGi, &absGj, &op1, &rho_op2, factor, sign = sign](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * rho_op2(j, i) * (-sign)), factor);
      };
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Drange())
          add_exp_weight(i, j);
      for (const auto j : diagIj.Krange())
        for (const auto i : diagIi.Drange())
          add_rho_weight(i, j);
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Krange())
          add_exp_weight(i, j);
    }
    void end([[maybe_unused]] const Step &step) override {
      spec.mergeCFS(*cb.get());
      cb.reset();
    }
    void save() override { if (should_save) spec.save(); }
    std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FDMgt : virtual public Algo<S> {
  private:
    inline static const std::string algoname = "FDMgt";
    SpectrumRealFreq<S> spec;
    const int sign; // 1 for bosons, -1 for fermions
    const bool should_save;
 protected:
   using CB = ChainBinning<S>;
   std::unique_ptr<CB> cb;
  public:
    using Algo<S>::P;
    Algo_FDMgt(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P, const bool save = true)
      : Algo<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname, save), P), sign(gf_sign(gt)), should_save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, [[maybe_unused]] const Invar &Ii, [[maybe_unused]] const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
    {
      const auto wnf     = stats.wnfactor[step.ndx()];
      const auto T       = P.T.value();
      const auto op2_rho = prod_fit(op2, rhoFDM.at(Ii));
      const auto absGi   = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi = absGi | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;
      const auto add_exp_weight = [this, &absGi, &absGj, &boltzGi, &op1, &op2, factor, wnf](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2(j, i) * wnf * boltzGi[i]), factor);
      };
      const auto add_rho_weight = [this, &absGi, &absGj, &op1, &op2_rho, factor](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2_rho(j, i)), factor);
      };
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Drange())
          add_exp_weight(i, j);
      for (const auto j : diagIj.Krange())
        for (const auto i : diagIi.Drange())
          add_exp_weight(i, j);
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Krange())
          add_rho_weight(i, j);
    }
    void end([[maybe_unused]] const Step &step) override {
      spec.mergeCFS(*cb.get());
      cb.reset();
    }
    void save() override { if (should_save) spec.save(); }
    std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
class Algo_FDM : public Algo_FDMls<S>, public Algo_FDMgt<S> {
 private:
    inline static const std::string algoname2 = "FDM";
    SpectrumRealFreq<S> spec_tot;
    const int sign;
  public:
    using Algo<S>::P;
    Algo_FDM(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
      Algo<S>(P), Algo_FDMls<S>(name, prefix, gt, P, false), Algo_FDMgt<S>(name, prefix, gt, P, false), spec_tot(name, algoname2, spec_fn(name, prefix, algoname2), P), sign(gf_sign(gt)) {}
    void begin(const Step &step) override {
      Algo_FDMgt<S>::begin(step);
      Algo_FDMls<S>::begin(step);
   }
    void calc(const Step &step, const Eigen<S> &diagIp, const Eigen<S> &diagI1, const Matrix &op1, const Matrix &op2,
              t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements<S> &rho,
              const Stats<S> &stats) override
    {
      auto &gt_cb = *Algo_FDMgt<S>::cb;
      auto &ls_cb = *Algo_FDMls<S>::cb;
      const auto wnf     = stats.wnfactor[step.ndx()];
      const auto T       = P.T.value();
      const auto rho_op2 = prod_fit(rho.at(I1), op2);
      const auto op2_rho = prod_fit(op2, rho.at(Ip));
      const auto absGi   = diagIp.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagI1.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi = absGi | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;
      const auto boltzGj = absGj | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;

      const auto add_exp_weight = [&gt_cb, &ls_cb, &absGi, &absGj, &boltzGi, &boltzGj, &op1, &op2, factor, wnf, sign = sign](const auto i, const auto j) {
        const auto energy = absGj[j] - absGi[i];
        const auto op1ji = conj_me(op1(j, i));
        const auto op2ji = op2(j, i);
        gt_cb.add(std::make_pair(energy, op1ji * op2ji * wnf * boltzGi[i]), factor);
        ls_cb.add(std::make_pair(energy, op1ji * op2ji * (-sign) * wnf * boltzGj[j]), factor);
      };
      const auto add_ls_rho_weight = [&ls_cb, &absGi, &absGj, &op1, &rho_op2, factor, sign = sign](const auto i, const auto j) {
        ls_cb.add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * rho_op2(j, i) * (-sign)), factor);
      };
      const auto add_gt_rho_weight = [&gt_cb, &absGi, &absGj, &op1, &op2_rho, factor](const auto i, const auto j) {
        gt_cb.add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2_rho(j, i)), factor);
      };

      for (const auto j : diagI1.Drange())
        for (const auto i : diagIp.Drange())
          add_exp_weight(i, j);
      for (const auto j : diagI1.Krange())
        for (const auto i : diagIp.Drange()) {
          const auto energy = absGj[j] - absGi[i];
          const auto op1ji = conj_me(op1(j, i));
          const auto op2ji = op2(j, i);
          gt_cb.add(std::make_pair(energy, op1ji * op2ji * wnf * boltzGi[i]), factor);
          add_ls_rho_weight(i, j);
        }
      for (const auto j : diagI1.Drange())
        for (const auto i : diagIp.Krange()) {
          add_gt_rho_weight(i, j);
          const auto energy = absGj[j] - absGi[i];
          const auto op1ji = conj_me(op1(j, i));
          const auto op2ji = op2(j, i);
          ls_cb.add(std::make_pair(energy, op1ji * op2ji * (-sign) * wnf * boltzGj[j]), factor);
        }
    }
    void end([[maybe_unused]] const Step &step) override {
      spec_tot.mergeCFS(*Algo_FDMgt<S>::cb.get());
      spec_tot.mergeCFS(*Algo_FDMls<S>::cb.get());
      Algo_FDMgt<S>::cb.reset();
      Algo_FDMls<S>::cb.reset();
    }
    void save() override { spec_tot.save(); }
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
      const auto T        = P.T.value();
      const auto rho_op2  = prod_fit(rhoFDM.at(Ij), op2);
      const auto op2_rho  = prod_fit(op2, rhoFDM.at(Ii));
      const auto absGi    = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj    = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi  = absGi | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;
      const auto boltzGj  = absGj | ranges::views::transform([T](const auto E) { return exp(-E / T); }) | ranges::to_vector;
      std::vector<t_weight> mats_freq(cutoff);
      for (size_t n = 0; n < cutoff; n++) mats_freq[n] = ww(static_cast<short>(n), gt, T) * 1i;
      const auto term1_factors = [&absGi, &absGj, &boltzGi, &boltzGj, &op1, &op2, wnf, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto op2ji = op2(j, i);
        const auto weightA = op1ji * op2ji * wnf * boltzGi[i];
        const auto weightB = op1ji * op2ji * (-sign) * wnf * boltzGj[j];
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto term1 = [T, this, &mats_freq](const auto energy, const auto weightA, const auto weightB, const auto n) -> t_weight {
         if (gt == gf_type::fermionic || n>0 || abs(energy) > WEIGHT_TOL)
           return (weightA + weightB) / (mats_freq[n] - energy);
         else // bosonic w=0 && Ei=Ej case
           return -weightA/T;
       };
      const auto term2_factors = [&absGi, &absGj, &boltzGi, &op1, &op2, &rho_op2, wnf, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto weightA = op1ji * op2(j, i) * wnf * boltzGi[i];
        const auto weightB = op1ji * rho_op2(j, i) * (-sign);
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto term2 = [&mats_freq](const auto energy, const auto weightA, const auto weightB, const auto n) {
        return (weightA + weightB) / (mats_freq[n] - energy);
      };
      const auto term3_factors = [&absGi, &absGj, &boltzGj, &op1, &op2, &op2_rho, wnf, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto weightA = op1ji * op2_rho(j, i);
        const auto weightB = (-sign) * op1ji * op2(j, i) * wnf * boltzGj[j];
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto term3 = [&mats_freq](const auto energy, const auto weightA, const auto weightB, const auto n) {
        return (weightA + weightB) / (mats_freq[n] - energy);
      };
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Drange()) {
          const auto factors = term1_factors(i, j);
          const auto energy = std::get<0>(factors);
          const auto weightA = std::get<1>(factors);
          const auto weightB = std::get<2>(factors);
#if NRG_ENABLE_APP_OPENMP
# pragma omp parallel for schedule(static)
#endif
          for (size_t n = 0; n < cutoff; n++)
            cm->add(n, term1(energy, weightA, weightB, n) * factor);
        }
      for (const auto j : diagIj.Krange())
        for (const auto i : diagIi.Drange()) {
          const auto factors = term2_factors(i, j);
          const auto energy = std::get<0>(factors);
          const auto weightA = std::get<1>(factors);
          const auto weightB = std::get<2>(factors);
#if NRG_ENABLE_APP_OPENMP
# pragma omp parallel for schedule(static)
#endif
          for (size_t n = 0; n < cutoff; n++)
            cm->add(n, term2(energy, weightA, weightB, n) * factor);
        }
      for (const auto j : diagIj.Drange())
        for (const auto i : diagIi.Krange()) {
          const auto factors = term3_factors(i, j);
          const auto energy = std::get<0>(factors);
          const auto weightA = std::get<1>(factors);
          const auto weightB = std::get<2>(factors);
#if NRG_ENABLE_APP_OPENMP
# pragma omp parallel for schedule(static)
#endif
          for (size_t n = 0; n < cutoff; n++)
            cm->add(n, term3(energy, weightA, weightB, n) * factor);
        }
    }
    void end([[maybe_unused]] const Step &step) override {
      gf.merge(*cm.get());
      cm.reset();
    }
    void save() override { gf.save(); }
    std::string rho_type() override { return "rhoFDM"; }
};

} // namespace

#endif
