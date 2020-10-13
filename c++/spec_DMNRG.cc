// OPTIMIZATION NOTE: the inner loop should involve the last index.

class Algo_DMNRG : public Algo {
 public:
   Algo_DMNRG(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "DMNRG"; }
   string merge() override { return "NN2"; }
   string rho_type() override { return "rho"; }
};

void Algo_DMNRG::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
  double Emin = P.getEmin(); // used in optimization
  double Emax = P.getEmax();
  if (P.ZBW) {
    Emin = 0;
    Emax = std::numeric_limits<double>::max(); // infinity
  }
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  auto dimp            = min(rhoNIp.size1(), diagIp.getnrstored());
  auto dim1            = min(rhoNI1.size1(), diagI1.getnrstored());
  for (const auto rm: range0(dimp)) {
    const t_eigen Em = diagIp.value_zero(rm);
    for (const auto rj: range0(dim1)) {
      const t_eigen Ej = diagI1.value_zero(rj);
      DELTA d;
      d.energy          = Ej - Em;
      const double absE = abs(d.energy);
      if (absE < Emin || absE > Emax) // does not contribute
        continue;
      weight_bucket sumA;
      for (const auto ri: range0(dimp)) sumA += op2II(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
      t_weight weightA = t_weight(sumA) * CONJ_ME(op1II(rj, rm));
      weight_bucket sumB;
      for (const auto ri: range0(dim1)) sumB += CONJ_ME(op1II(ri, rm)) * rhoNI1(rj, ri); // non-optimal
      t_weight weightB = t_weight(sumB) * op2II(rj, rm);
      d.weight         = spinfactor * (weightA + (-sign) * weightB);
      cs->add(step.scale() * d.energy, d.weight);
    }
  }
}

class Algo_DMNRGmats : public Algo {
 public:
   Algo_DMNRGmats(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &bs) override { return make_shared<ChainSpectrumMatsubara>(P, bs.mt); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "DMNRGmats"; }
   string rho_type() override { return "rho"; }
};

void Algo_DMNRGmats::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs,
                          t_factor spinfactor, spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
  auto csm   = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  auto dimp            = min(rhoNIp.size1(), diagIp.getnrstored());
  auto dim1            = min(rhoNI1.size1(), diagI1.getnrstored());
  for (const auto rm: range0(dimp)) {
    const t_eigen Em = diagIp.value_zero(rm);
    for (const auto rj: range0(dim1)) {
      const t_eigen Ej = diagI1.value_zero(rj);
      DELTA d;
      d.energy = Ej - Em;
      weight_bucket sumA;
      for (const auto ri: range0(dimp)) sumA += op2II(rj, ri) * rhoNIp(rm, ri); // rm <-> ri, rho symmetric
      t_weight weightA = t_weight(sumA) * CONJ_ME(op1II(rj, rm));
      weight_bucket sumB;
      for (const auto ri: range0(dim1)) sumB += CONJ_ME(op1II(ri, rm)) * rhoNI1(rj, ri); // non-optimal
      t_weight weightB = t_weight(sumB) * op2II(rj, rm);
      d.weight         = spinfactor * (weightA + (-sign) * weightB);
      for (size_t n = 1; n < P.mats; n++) csm->add(n, d.weight / (cmpl(0, ww(n, bs.mt, P.T)) - step.scale() * d.energy));
      if (abs(d.energy) > WEIGHT_TOL || bs.mt == matstype::fermionic)
        csm->add(size_t(0), d.weight / (cmpl(0, ww(0, bs.mt, P.T)) - step.scale() * d.energy));
      else // bosonic w=0 && E1=Ep case
        csm->add(size_t(0), spinfactor * (-weightA / t_weight(P.T)));
    }
  }
}
