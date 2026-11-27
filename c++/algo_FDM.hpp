#ifndef _algo_FDM_hpp_
#define _algo_FDM_hpp_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include "traits.hpp"
#include "algo.hpp"
#include "fdm_zero_dd.hpp"
#include "spectrum.hpp"

namespace NRG {

using namespace std::complex_literals;

// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

namespace detail {

// Only the kept-row, discarded-column block of rho * op enters the spectrum.
template<scalar S>
auto fdm_rho_times_op(const DensMatElements<S> &rho, const Invar &row, const Matrix_traits<S> &op,
                      const Eigen<S> &row_diag, const Eigen<S> &col_diag) {
  const auto kept_rows = row_diag.Krange().size();
  const auto kept_cols = col_diag.Krange().size();
  const auto cols = col_diag.getnrall();
  if (!kept_rows || kept_cols == cols) return Matrix_traits<S>{};
  const auto &density = rho.at(row);
  const auto inner = std::min(size2(density), size1(op));
  return matrix_prod<S>(submatrix_const(density, {0, kept_rows}, {0, inner}),
                        submatrix_const(op, {0, inner}, {kept_cols, cols}));
}

// Only the discarded-row, kept-column block of op * rho enters the spectrum.
template<scalar S>
auto fdm_op_times_rho(const Matrix_traits<S> &op, const DensMatElements<S> &rho, const Invar &col,
                      const Eigen<S> &row_diag, const Eigen<S> &col_diag) {
  const auto kept_rows = row_diag.Krange().size();
  const auto kept_cols = col_diag.Krange().size();
  const auto rows = row_diag.getnrall();
  if (kept_rows == rows || !kept_cols) return Matrix_traits<S>{};
  const auto &density = rho.at(col);
  const auto inner = std::min(size2(op), size1(density));
  return matrix_prod<S>(submatrix_const(op, {kept_rows, rows}, {0, inner}),
                        submatrix_const(density, {0, inner}, {0, kept_cols}));
}

} // namespace detail

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
             t_coef factor, const Invar &, const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
    {
      const auto wnf     = stats.wnfactor[step.ndx()];
      const auto T       = P.T.value();
      const auto rho_op2 = detail::fdm_rho_times_op(rhoFDM, Ij, op2, diagIj, diagIi);
      const auto ki = diagIi.Krange().size();
      const auto absGi   = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGj = absGj | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;
      const auto add_exp_weight = [this, &absGi, &absGj, &boltzGj, &op1, &op2, factor, weight_scale = (-sign) * wnf](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2(j, i) * weight_scale * boltzGj[j]), factor);
      };
      const auto add_rho_weight = [this, &absGi, &absGj, &op1, &rho_op2, factor, ki, sign_ = sign](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * rho_op2(j, i - ki) * (-sign_)), factor);
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
             t_coef factor, const Invar &Ii, const Invar &, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
    {
      const auto wnf     = stats.wnfactor[step.ndx()];
      const auto T       = P.T.value();
      const auto op2_rho = detail::fdm_op_times_rho(op2, rhoFDM, Ii, diagIj, diagIi);
      const auto kj = diagIj.Krange().size();
      const auto absGi   = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi = absGi | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;
      const auto add_exp_weight = [this, &absGi, &absGj, &boltzGi, &op1, &op2, factor, wnf](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2(j, i) * wnf * boltzGi[i]), factor);
      };
      const auto add_rho_weight = [this, &absGi, &absGj, &op1, &op2_rho, factor, kj](const auto i, const auto j) {
        cb->add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2_rho(j - kj, i)), factor);
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
      const auto rho_op2 = detail::fdm_rho_times_op(rho, I1, op2, diagI1, diagIp);
      const auto op2_rho = detail::fdm_op_times_rho(op2, rho, Ip, diagI1, diagIp);
      const auto ki = diagIp.Krange().size(), kj = diagI1.Krange().size();
      const auto absGi   = diagIp.values.all_abs_G() | ranges::to_vector;
      const auto absGj   = diagI1.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi = absGi | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;
      const auto boltzGj = absGj | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;

      const auto add_exp_weight = [&gt_cb, &ls_cb, &absGi, &absGj, &boltzGi, &boltzGj, &op1, &op2, factor, wnf, sign_ = sign](const auto i, const auto j) {
        const auto energy = absGj[j] - absGi[i];
        const auto op1ji = conj_me(op1(j, i));
        const auto op2ji = op2(j, i);
        gt_cb.add(std::make_pair(energy, op1ji * op2ji * wnf * boltzGi[i]), factor);
        ls_cb.add(std::make_pair(energy, op1ji * op2ji * (-sign_) * wnf * boltzGj[j]), factor);
      };
      const auto add_ls_rho_weight = [&ls_cb, &absGi, &absGj, &op1, &rho_op2, factor, ki, sign_ = sign](const auto i, const auto j) {
        ls_cb.add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * rho_op2(j, i - ki) * (-sign_)), factor);
      };
      const auto add_gt_rho_weight = [&gt_cb, &absGi, &absGj, &op1, &op2_rho, factor, kj](const auto i, const auto j) {
        gt_cb.add(std::make_pair(absGj[j] - absGi[i], conj_me(op1(j, i)) * op2_rho(j - kj, i)), factor);
      };

      const auto di = diagIp.Drange(), dj = diagI1.Drange();
      if (!fdm_zero_dd_block(di, dj, boltzGi, boltzGj, absGi, absGj, op1, op2, T, wnf, factor)) {
        for (const auto j : dj)
          for (const auto i : di)
            add_exp_weight(i, j);
      }
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
   Algo_FDMmats(const std::string &name, const std::string &prefix, const gf_type gt_, const Params &P_) :
     Algo<S>(P_), gf(name, algoname, spec_fn(name, prefix, algoname), gt_, P_), sign(gf_sign(gt_)), gt(gt_) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen<S> &diagIi, const Eigen<S> &diagIj, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements<S> &rhoFDM,
             const Stats<S> &stats) override
    {
      const size_t cutoff = P.mats;
      if (!cutoff) return;
      const auto wnf      = stats.wnfactor[step.ndx()];
      const auto T        = P.T.value();
      const auto rho_op2  = detail::fdm_rho_times_op(rhoFDM, Ij, op2, diagIj, diagIi);
      const auto op2_rho  = detail::fdm_op_times_rho(op2, rhoFDM, Ii, diagIj, diagIi);
      const auto ki = diagIi.Krange().size(), kj = diagIj.Krange().size();
      const auto di = diagIi.Drange().size(), dj = diagIj.Drange().size();
      const auto absGi    = diagIi.values.all_abs_G() | ranges::to_vector;
      const auto absGj    = diagIj.values.all_abs_G() | ranges::to_vector;
      const auto boltzGi  = absGi | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;
      const auto boltzGj  = absGj | ranges::views::transform([T](const auto E) { return std::exp(-E / T); }) | ranges::to_vector;
      const auto dd_factors = [&absGi, &absGj, &boltzGi, &boltzGj, &op1, &op2, wnf, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto op2ji = op2(j, i);
        const auto weightA = op1ji * op2ji * wnf * boltzGi[i];
        const auto weightB = op1ji * op2ji * (-sign) * wnf * boltzGj[j];
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto kd_factors = [&absGi, &absGj, &boltzGi, &op1, &op2, &rho_op2, wnf, ki, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto weightA = op1ji * op2(j, i) * wnf * boltzGi[i];
        const auto weightB = op1ji * rho_op2(j, i - ki) * (-sign);
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto dk_factors = [&absGi, &absGj, &boltzGj, &op1, &op2, &op2_rho, wnf, kj, this](const auto i, const auto j) {
        const auto Ei = absGi[i];
        const auto Ej = absGj[j];
        const auto energy = Ej - Ei;
        const auto op1ji = conj_me(op1(j, i));
        const auto weightA = op1ji * op2_rho(j - kj, i);
        const auto weightB = (-sign) * op1ji * op2(j, i) * wnf * boltzGj[j];
        return std::make_tuple(energy, weightA, weightB);
      };
      const auto dd_count = dj * di, kd_count = kj * di, dk_count = dj * ki;
      // Traverse DD, KD, DK in order, with i varying fastest in each sector.
      detail::accumulate_matsubara_transitions(cutoff, dd_count + kd_count + dk_count,
        [&](size_t pair) {
          if (pair < dd_count) return dd_factors(ki + pair % di, kj + pair / di);
          pair -= dd_count;
          if (pair < kd_count) return kd_factors(ki + pair % di, pair / di);
          pair -= kd_count;
          return dk_factors(pair % ki, kj + pair / ki);
        },
        [&](const size_t n, const auto &transition) {
          const auto &[energy, weightA, weightB] = transition;
          const t_weight term = gt == gf_type::fermionic || n > 0 || std::abs(energy) > WEIGHT_TOL
                              ? (weightA + weightB) / (cm->frequency(n) * 1i - energy) : t_weight(-weightA / T);
          cm->add(n, term * factor);
        });
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
