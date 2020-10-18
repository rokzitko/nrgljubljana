// OPTIMIZATION NOTE: the inner loop should involve the last index.

class Algo_DMNRG : public Algo {
 public:
   explicit Algo_DMNRG(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_coef, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "DMNRG"; }
   string merge() override { return "NN2"; }
   string rho_type() override { return "rho"; }
};

void Algo_DMNRG::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_coef spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const auto sign = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
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
      for (const auto ri: diagIp.kept()) sumA += op2II(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
      const auto weightA = t_weight(sumA) * conj_me(op1II(rj, rm));
      t_weight sumB{};
      for (const auto ri: diagI1.kept()) sumB += conj_me(op1II(ri, rm)) * rhoNI1(rj, ri); // non-optimal
      const auto weightB = t_weight(sumB) * op2II(rj, rm);
      const auto weight  = spinfactor * (weightA + (-sign) * weightB);
      cs->add(step.scale() * energy, weight);
    }
  }
}

class Algo_DMNRGmats : public Algo {
 public:
   explicit Algo_DMNRGmats(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &bs) override { return make_shared<ChainSpectrumMatsubara>(P, bs.mt); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_coef, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "DMNRGmats"; }
   string rho_type() override { return "rho"; }
};

void Algo_DMNRGmats::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs,
                          t_coef spinfactor, spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const auto sign = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
  auto csm   = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  for (const auto rm: diagIp.kept()) {
    for (const auto rj: diagI1.kept()) {
      const auto Em = diagIp.value_zero(rm);
      const auto Ej = diagI1.value_zero(rj);
      const auto energy = Ej - Em;
      t_weight sumA{};
      for (const auto ri: diagIp.kept()) sumA += op2II(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
      const auto weightA = sumA * conj_me(op1II(rj, rm));
      t_weight sumB{};
      for (const auto ri: diagI1.kept()) sumB += conj_me(op1II(ri, rm)) * rhoNI1(rj, ri); // non-optimal
      const auto weightB = sumB * op2II(rj, rm);
      const auto weight  = spinfactor * (weightA + (-sign) * weightB);
      for (size_t n = 1; n < P.mats; n++) csm->add(n, weight / (cmpl(0, ww(n, bs.mt, P.T)) - step.scale() * energy));
      if (abs(energy) > WEIGHT_TOL || bs.mt == matstype::fermionic)
        csm->add(size_t(0), weight / (cmpl(0, ww(0, bs.mt, P.T)) - step.scale() * energy));
      else // bosonic w=0 && E1=Ep case
        csm->add(size_t(0), spinfactor * (-weightA / t_weight(P.T)));
    }
  }
}
