// OPTIMIZATION NOTE: the inner loop should involve the last index.

template<typename S>
class Algo_DMNRG_tmpl : public Algo_tmpl<S> {
 private:
   inline static const std::string algoname = "DMNRG";
   SpectrumRealFreq_tmpl<S> spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning_tmpl<S>;
   std::unique_ptr<CB> cb;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using t_weight = typename traits<S>::t_weight;
   using Algo_tmpl<S>::P;
   Algo_DMNRG_tmpl(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo_tmpl<S>(P), spec(name, algoname, spec_fn(name, prefix, algoname), P), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, const Stats_tmpl<S> &stats) override
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
   ~Algo_DMNRG_tmpl() { spec.save(); }
   std::string rho_type() override { return "rho"; }
};

template<typename S>
class Algo_DMNRGmats_tmpl : public Algo_tmpl<S> {
 private:
   inline static const std::string algoname = "DMNRGmats";
   GFMatsubara_tmpl<S> gf;
   const int sign;
   const gf_type gt;
   using CM = ChainMatsubara_tmpl<S>;
   std::unique_ptr<CM> cm;
 public:
   using Matrix = typename traits<S>::Matrix;
   using t_coef = typename traits<S>::t_coef;
   using t_eigen = typename traits<S>::t_eigen;
   using t_weight = typename traits<S>::t_weight;
   using Algo_tmpl<S>::P;
   Algo_DMNRGmats_tmpl(const std::string &name, const std::string &prefix, const gf_type gt, const Params &P) :
     Algo_tmpl<S>(P), gf(name, algoname, spec_fn(name, prefix, algoname), gt, P), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen_tmpl<S> &diagIp, const Eigen_tmpl<S> &diagI1, const Matrix &op1, const Matrix &op2,
             t_coef factor, const Invar &Ip, const Invar &I1, const DensMatElements_tmpl<S> &rho, const Stats_tmpl<S> &stats) override
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
         for (size_t n = 1; n < P.mats; n++) cm->add(n, weight / (cmpl(0, ww(n, gt, P.T)) - step.scale() * energy));
         if (abs(energy) > WEIGHT_TOL || gt == gf_type::fermionic)
           cm->add(size_t(0), weight / (cmpl(0, ww(0, gt, P.T)) - step.scale() * energy));
         else // bosonic w=0 && E1=Ep case
           cm->add(size_t(0), factor * (-weightA / t_weight(P.T)));
       }
     }
   }
   void end(const Step &step) override {
          gf.merge(*cm.get());
          cm.reset();
   }
   ~Algo_DMNRGmats_tmpl() { gf.save(); }
   std::string rho_type() override { return "rho"; }
};
