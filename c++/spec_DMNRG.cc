// OPTIMIZATION NOTE: the inner loop should involve the last index.

class Algo_DMNRG : public Algo {
 private:
   SpectrumRealFreq spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning;
   std::unique_ptr<CB> cb;
 public:
   explicit Algo_DMNRG(SpectrumRealFreq spec, gf_type gt, const Params &P) : Algo(P), spec(spec), sign(gf_sign(gt)) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
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
         const auto weightA = t_weight(sumA) * conj_me(op1(rj, rm));
         t_weight sumB{};
         for (const auto ri: diagI1.kept()) sumB += conj_me(op1(ri, rm)) * rhoNI1(rj, ri); // non-optimal
         const auto weightB = t_weight(sumB) * op2(rj, rm);
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

class Algo_DMNRGmats : public Algo {
 private:
   GFMatsubara gf;
   const int sign;
   gf_type gt;
   using CM = ChainMatsubara;
   std::unique_ptr<CM> cm;
 public:
   explicit Algo_DMNRGmats(GFMatsubara gf, gf_type gt, const Params &P) : Algo(P), gf(gf), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
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
   ~Algo_DMNRGmats() { gf.save(); }
   std::string rho_type() override { return "rho"; }
};
