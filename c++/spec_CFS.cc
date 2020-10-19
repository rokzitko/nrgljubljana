// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).
// Based on the implementation by Markus Greger.

class Algo_CFSls : virtual public Algo {
 protected:
   SpectrumRealFreq spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   explicit Algo_CFSls(SpectrumRealFreq spec, gf_type gt, const Params &P, const bool save = true)
     : Algo(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
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
   ~Algo_CFSls() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};

class Algo_CFSgt : virtual public Algo {
 protected:
   SpectrumRealFreq spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   explicit Algo_CFSgt(SpectrumRealFreq spec, gf_type gt, const Params &P, const bool save = true)
     : Algo(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
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
         if constexpr (std::is_same_v<t_matel, double>) {
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
   ~Algo_CFSgt() { if (save) spec.save(); }
   std::string rho_type() override { return "rho"; }
};

class Algo_CFS : public Algo_CFSls, public Algo_CFSgt {
 private:
    SpectrumRealFreq spec_tot;
 public:
   explicit Algo_CFS(SpectrumRealFreq spec, gf_type gt, const Params &P) :
     Algo(P), Algo_CFSls(spec, gt, P, false), Algo_CFSgt(spec, gt, P, false), spec_tot(spec) {}
   void begin(const Step &step) override {
     Algo_CFSgt::begin(step);
     Algo_CFSls::begin(step);
   }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
   {
     Algo_CFSgt::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_CFSls::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end(const Step &step) override {
     spec_tot.mergeCFS(*Algo_CFSgt::cb.get());
     spec_tot.mergeCFS(*Algo_CFSls::cb.get());
     Algo_CFSgt::cb.reset();
     Algo_CFSls::cb.reset();
   }
   ~Algo_CFS() { spec_tot.save(); }
   std::string rho_type() override { return "rho"; }
};

