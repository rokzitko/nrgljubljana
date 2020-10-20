// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

#define LOOP_D(n)                                                                                                                                    \
  for (size_t n = ret##n; n < all##n; n++) {                                                                                                         \
    const t_eigen E##n = diagI##n.absenergyG(n);
#define LOOP_K(n)                                                                                                                                    \
  for (size_t n = 0; n < ret##n; n++) {                                                                                                              \
    const t_eigen E##n = diagI##n.absenergyG(n);

template<typename S>
class Algo_FDMls_tmpl : virtual public Algo_tmpl<S> {
 protected:
   SpectrumRealFreq_tmpl<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning_tmpl<S>;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo_tmpl<S>::P;
   explicit Algo_FDMls_tmpl(SpectrumRealFreq_tmpl<S> spec, gf_type gt, const Params &P, const bool save = true) 
     : Algo_tmpl<S>(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIi, const Eigen_tmpl<S> &diagIj, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements_tmpl<S> &rhoFDM, 
             const Stats_tmpl<S> &stats) override
   {
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = step.last() ? 0 : diagIi.getnrkept();
     const auto retj  = step.last() ? 0 : diagIj.getnrkept();
     const auto alli  = diagIi.getnrstored();
     const auto allj  = diagIj.getnrstored();
     LOOP_D(i)
       LOOP_D(j) // B3
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T);
         cb->add(energy, weight);
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [K=retj, D=alli]
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, retj), ublas::range(0, alli));
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
   ~Algo_FDMls_tmpl() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};
// AAA using Algo_FDMls = Algo_FDMls_tmpl<scalar>;

template<typename S>
class Algo_FDMgt_tmpl : virtual public Algo_tmpl<S> {
 protected:
   SpectrumRealFreq_tmpl<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning_tmpl<S>;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo_tmpl<S>::P;
   explicit Algo_FDMgt_tmpl(SpectrumRealFreq_tmpl<S> spec, const gf_type gt, const Params &P, const bool save = true) 
     : Algo_tmpl<S>(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIi, const Eigen_tmpl<S> &diagIj, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements_tmpl<S> &rhoFDM, 
             const Stats_tmpl<S> &stats) override
   {
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = step.last() ? 0 : diagIi.getnrkept();
     const auto retj  = step.last() ? 0 : diagIj.getnrkept();
     const auto alli  = diagIi.getnrstored();
     const auto allj  = diagIj.getnrstored();
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
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, allj), ublas::range(0, reti));
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
   ~Algo_FDMgt_tmpl() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};
// AAA using Algo_FDMgt = Algo_FDMgt_tmpl<scalar>;

template<typename S>
class Algo_FDM_tmpl : public Algo_FDMls_tmpl<S>, public Algo_FDMgt_tmpl<S> {
 private:
   SpectrumRealFreq_tmpl<S> spec_tot;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo_tmpl<S>::P;
   explicit Algo_FDM_tmpl(SpectrumRealFreq_tmpl<S> spec, const gf_type gt, const Params &P) :
     Algo_tmpl<S>(P), Algo_FDMls_tmpl<S>(spec, gt, P, false), Algo_FDMgt_tmpl<S>(spec, gt, P, false), spec_tot(spec) {}
   void begin(const Step &step) override {
     Algo_FDMgt_tmpl<S>::begin(step);
     Algo_FDMls_tmpl<S>::begin(step);
   }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, 
             const Stats_tmpl<S> &stats) override
   {
     Algo_FDMgt_tmpl<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_FDMls_tmpl<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end(const Step &step) override {
     spec_tot.mergeCFS(*Algo_FDMgt_tmpl<S>::cb.get());
     spec_tot.mergeCFS(*Algo_FDMls_tmpl<S>::cb.get());
     Algo_FDMgt_tmpl<S>::cb.reset();
     Algo_FDMls_tmpl<S>::cb.reset();
   }
   ~Algo_FDM_tmpl() { spec_tot.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};
// AAA using Algo_FDM = Algo_FDM_tmpl<scalar>;

template<typename S>
class Algo_FDMmats_tmpl : public Algo_tmpl<S> {
 private:
   GFMatsubara_tmpl<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara_tmpl<S>;
   std::unique_ptr<CM> cm;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo_tmpl<S>::P;
   explicit Algo_FDMmats_tmpl(GFMatsubara_tmpl<S> gf, const gf_type gt, const Params &P) : 
     Algo_tmpl<S>(P), gf(gf), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIi, const Eigen_tmpl<S> &diagIj, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ii, const Invar &Ij, const DensMatElements_tmpl<S> &rhoFDM, 
             const Stats_tmpl<S> &stats) override
   {
     const size_t cutoff = P.mats;
     const auto wnf   = stats.wnfactor[step.ndx()];
     const auto &rhoi = rhoFDM.at(Ii);
     const auto &rhoj = rhoFDM.at(Ij);
     const auto reti  = (step.last() ? 0 : rhoi.size1());
     const auto retj  = (step.last() ? 0 : rhoj.size1());
     const auto alli  = diagIi.getnrstored();
     const auto allj  = diagIj.getnrstored();
     LOOP_D(i)
       LOOP_D(j)
         const auto energy = Ej - Ei;
         const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T); // a[ij] b[ji] exp(-beta e[i])
         const auto weightB = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T); // a[ij] b[ji] sign exp(-beta e[j])
         #pragma omp parallel for schedule(static)
         for (size_t n = 1; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
         if (gt == gf_type::fermionic || abs(energy) > WEIGHT_TOL)
           cm->add(size_t(0), (weightA + weightB) / (cmpl(0, ww(0, gt, P.T)) - energy));
         else // bosonic w=0 && Ei=Ej case
           cm->add(size_t(0), (-weightA / P.T));
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [retj, alli]
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, retj), ublas::range(0, alli));
       Matrix rho_op2(retj, alli);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
       LOOP_D(i)
         LOOP_K(j)
           const auto energy = Ej - Ei;
           const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
           const auto weightB = factor * conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
           #pragma omp parallel for schedule(static)
           for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
         }
       }
     }
     if (allj > 0 && reti > 0) {
        // B [allj,reti] x rho [reti,reti]
        const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, allj), ublas::range(0, reti));
        Matrix op2_rho(allj, reti);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
        LOOP_K(i)
          LOOP_D(j)
            const auto energy = Ej - Ei;
            const auto weightA = factor * conj_me(op1(j, i)) * op2_rho(j, i);
            const auto weightB = factor * (-sign) * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ej / P.T);
            #pragma omp parallel for schedule(static)
            for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
          }
        }
     }
   }
   void end(const Step &step) override {
     gf.merge(*cm.get());
     cm.reset();
   }
   ~Algo_FDMmats_tmpl() { gf.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};
// AAA using Algo_FDMmats = Algo_FDMmats_tmpl<scalar>;
