// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).
// Based on the implementation by Markus Greger.

template<typename S>
class Algo_CFSls_tmpl : virtual public Algo_tmpl<S> {
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
   explicit Algo_CFSls_tmpl(SpectrumRealFreq_tmpl<S> spec, const gf_type gt, const Params &P, const bool save = true)
     : Algo_tmpl<S>(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, const Stats_tmpl<S> &stats) override
   {
     const auto &rhoNIp = rho.at(Ip);
     const auto &rhoNI1 = rho.at(I1);
     // Convention: k-loops over retained states, l-loop over discarded states. i-term, Eq. (11).
     if (step.last()) {
       for (const auto r1: diagI1.kept()) {
         for (const auto rp: diagIp.kept()) {
           const auto E1 = diagI1.value_zero(r1);
           const auto Ep = diagIp.value_zero(rp);
           const auto weight = (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * exp(-E1 * step.scT()) * (-sign);
           cb->add(step.scale() * (E1-Ep), weight);
         }
       }
     } else {
       // iii-term, Eq. (16), positive frequency excitations
       if (op2.size1() && rhoNIp.size1()) {
         Matrix op2_m_rho;
         const ublas::matrix_range<const Matrix> op2_TK(op2, ublas::range(0, op2.size1()), ublas::range(0, rhoNIp.size1()));
         op2_m_rho = Matrix(op2_TK.size1(), rhoNIp.size2());
         atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2_TK, rhoNIp, 0.0, op2_m_rho); // rhoNEW <- rhoNEW + factor T U
         for (const auto rl: diagI1.discarded()) {
           for (const auto rk: diagIp.kept()) {
             const auto El       = diagI1.value_zero(rl);
             const auto Ek       = diagIp.value_zero(rk);
             const auto weight = factor * conj_me(op1(rl, rk)) * op2_m_rho(rl, rk) * (-sign);
             cb->add(step.scale() * (El-Ek), weight);
           }
         }
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_CFSls_tmpl() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};
// AAA using Algo_CFSls = Algo_CFSls_tmpl<scalar>;

template<typename S>
class Algo_CFSgt_tmpl : virtual public Algo_tmpl<S> {
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
   explicit Algo_CFSgt_tmpl(SpectrumRealFreq_tmpl<S> spec, const gf_type gt, const Params &P, const bool save = true)
     : Algo_tmpl<S>(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, const Stats_tmpl<S> &stats) override
   {
     const auto &rhoNIp = rho.at(Ip);
     const auto &rhoNI1 = rho.at(I1);
     // Convention: k-loops over retained states, l-loop over discarded states. i-term, Eq. (11).
     if (step.last()) {
       for (const auto r1: diagI1.kept()) {
         for (const auto rp: diagIp.kept()) {
           const auto E1 = diagI1.value_zero(r1);
           const auto Ep = diagIp.value_zero(rp);
           const auto weight = (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * exp(-Ep * step.scT());
           cb->add(step.scale() * (E1-Ep), weight);
         }
       }
     } else {
       if (rhoNI1.size1() && op1.size2()) {
         const ublas::matrix_range<const Matrix> op1_KT(op1, ublas::range(0, rhoNI1.size1()), ublas::range(0, op1.size2()));
         Matrix op1_m_rho(rhoNI1.size2(), op1_KT.size2());
         if constexpr (std::is_same_v<S, double>) {
           atlas::gemm(CblasTrans, CblasNoTrans, 1.0, rhoNI1, op1_KT, 0.0, op1_m_rho); // rhoNEW <- rhoNEW + factor T U
         } else {
           const ublas::matrix<std::complex<double>> conj_op1_KT = conj(op1_KT);
           atlas::gemm(CblasTrans, CblasNoTrans, 1.0, rhoNI1, conj_op1_KT, 0.0, op1_m_rho); // rhoNEW <- rhoNEW + factor T U
         }
         for (const auto rk: diagI1.kept()) {                                          // ii-term, Eq. (15), negative frequency excitations
           for (const auto rl: diagIp.discarded()) {
             const auto Ek       = diagI1.value_zero(rk);
             const auto El       = diagIp.value_zero(rl);
             const auto weight = factor * op1_m_rho(rk, rl) * op2(rk, rl);
             cb->add(step.scale() * (Ek-El), weight);
           }
         }
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_CFSgt_tmpl() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};
// AAA using Algo_CFSgt = Algo_CFSgt_tmpl<scalar>;

template<typename S>
class Algo_CFS_tmpl : public Algo_CFSls_tmpl<S>, public Algo_CFSgt_tmpl<S> {
 private:
    SpectrumRealFreq_tmpl<S> spec_tot;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using Algo_tmpl<S>::P;
   explicit Algo_CFS_tmpl(SpectrumRealFreq_tmpl<S> spec, gf_type gt, const Params &P) :
     Algo_tmpl<S>(P), Algo_CFSls_tmpl<S>(spec, gt, P, false), Algo_CFSgt_tmpl<S>(spec, gt, P, false), spec_tot(spec) {}
   void begin(const Step &step) override {
     Algo_CFSgt_tmpl<S>::begin(step);
     Algo_CFSls_tmpl<S>::begin(step);
   }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2, 
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, const Stats_tmpl<S> &stats) override
   {
     Algo_CFSgt_tmpl<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_CFSls_tmpl<S>::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end(const Step &step) override {
     spec_tot.mergeCFS(*Algo_CFSgt_tmpl<S>::cb.get());
     spec_tot.mergeCFS(*Algo_CFSls_tmpl<S>::cb.get());
     Algo_CFSgt_tmpl<S>::cb.reset();
     Algo_CFSls_tmpl<S>::cb.reset();
   }
   ~Algo_CFS_tmpl() { spec_tot.save(); }
   std::string rho_type() override { return "rho"; }
};
// AAA using Algo_CFS = Algo_CFS_tmpl<scalar>;
