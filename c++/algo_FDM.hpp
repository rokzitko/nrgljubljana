#ifndef _algo_FDM_hpp_
#define _algo_FDM_hpp_

#include <complex>
#include "traits.hpp"
#include "algo.hpp"
#include "spectrum.hpp"

namespace NRG {

using namespace std::complex_literals;

// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

#define LOOP_D(n)                                                                                                                                    \
  for (size_t n = ret##n; n < all##n; n++) {                                                                                                         \
    const t_eigen E##n = diagI##n.values.abs_G(n);
#define LOOP_K(n)                                                                                                                                    \
  for (size_t n = 0; n < ret##n; n++) {                                                                                                              \
    const t_eigen E##n = diagI##n.values.abs_G(n);

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
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = step.last() ? 0 : diagIi.getnrkept();
     const auto retj  = step.last() ? 0 : diagIj.getnrkept();
     const auto alli  = diagIi.getnrall();
     const auto allj  = diagIj.getnrall();
     LOOP_D(i)
       LOOP_D(j) // B3
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T);
         cb->add(energy, weight);
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [K=retj, D=alli]
       const auto op2cut = submatrix(op2, {0, retj}, {0, alli});
       Matrix rho_op2(retj, alli);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
       LOOP_D(i)
         LOOP_K(j)
           const auto energy = Ej - Ei;
           const auto weight = factor * conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
           cb->add(energy, weight);
         }
       }
     }
     LOOP_K(i)
       LOOP_D(j)
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T);
         cb->add(energy, weight);
       }
     }
   }
   void end(const Step &step) override {
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
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = step.last() ? 0 : diagIi.getnrkept();
     const auto retj  = step.last() ? 0 : diagIj.getnrkept();
     const auto alli  = diagIi.getnrall();
     const auto allj  = diagIj.getnrall();
     LOOP_D(i)
       LOOP_D(j)
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
         cb->add(energy, weight);
       }
     }
     LOOP_D(i)
       LOOP_K(j)
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
         cb->add(energy, weight);
       }
     }
     if (allj > 0 && reti > 0) {
       // B [D=allj,K=reti] x rho [reti,reti]
       const auto op2cut = submatrix(op2, {0, allj}, {0, reti});
       Matrix op2_rho(allj, reti);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
       LOOP_K(i)
         LOOP_D(j)
           const auto energy = Ej - Ei;
           const auto weight = factor * conj_me(op1(j, i)) * op2_rho(j, i);
           cb->add(energy, weight);
         }
       }
     }
   }
   void end(const Step &step) override {
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
   void end(const Step &step) override {
     spec_tot.mergeCFS(*Algo_FDMgt<S>::cb.get());
     spec_tot.mergeCFS(*Algo_FDMls<S>::cb.get());
     Algo_FDMgt<S>::cb.reset();
     Algo_FDMls<S>::cb.reset();
   }
   ~Algo_FDM() { spec_tot.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>, typename t_eigen = eigen_traits<S>>
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
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = (step.last() ? 0 : rhoi.size1());
     const auto retj  = (step.last() ? 0 : rhoj.size1());
     const auto alli  = diagIi.getnrall();
     const auto allj  = diagIj.getnrall();
     LOOP_D(i)
       LOOP_D(j)
         const auto energy = Ej - Ei;
         const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T); // a[ij] b[ji] exp(-beta e[i])
         const auto weightB = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T); // a[ij] b[ji] sign exp(-beta e[j])
         #pragma omp parallel for schedule(static)
         for (size_t n = 1; n < cutoff; n++) cm->add(n, (weightA + weightB) / (ww(n, gt, P.T)*1i - energy));
         if (gt == gf_type::fermionic || abs(energy) > WEIGHT_TOL)
           cm->add(size_t(0), (weightA + weightB) / (ww(0, gt, P.T)*1i - energy));
         else // bosonic w=0 && Ei=Ej case
           cm->add(size_t(0), (-weightA / (double)P.T));
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [retj, alli]
       const auto op2cut = submatrix(op2, {0, retj}, {0, alli});
       Matrix rho_op2(retj, alli);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
       LOOP_D(i)
         LOOP_K(j)
           const auto energy = Ej - Ei;
           const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
           const auto weightB = factor * conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
           #pragma omp parallel for schedule(static)
           for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (ww(n, gt, P.T)*1i - energy));
         }
       }
     }
     if (allj > 0 && reti > 0) {
        // B [allj,reti] x rho [reti,reti]
        const auto op2cut = submatrix(op2, {0, allj}, {0, reti});
        Matrix op2_rho(allj, reti);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
        LOOP_K(i)
          LOOP_D(j)
            const auto energy = Ej - Ei;
            const auto weightA = factor * conj_me(op1(j, i)) * op2_rho(j, i);
            const auto weightB = factor * (-sign) * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ej / P.T);
            #pragma omp parallel for schedule(static)
            for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (ww(n, gt, P.T)*1i - energy));
          }
        }
     }
   }
   void end(const Step &step) override {
     gf.merge(*cm.get());
     cm.reset();
   }
   ~Algo_FDMmats() { gf.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

} // namespace

#endif
