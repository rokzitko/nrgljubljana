class Algo_CFSls : virtual public Algo {
 public:
   explicit Algo_CFSls(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "CFSls"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rho"; }
};

class Algo_CFSgt : virtual public Algo {
 public:
   explicit Algo_CFSgt(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "CFSgt"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rho"; }
};

class Algo_CFS : public Algo_CFSls, public Algo_CFSgt {
 public:
   explicit Algo_CFS(const Params &P) : Algo(P), Algo_CFSls(P), Algo_CFSgt(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
             spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const override {
               Algo_CFSgt::calc(step, diagIp, diagI1, op1II, op2II, bs, spinfactor, cs, Ip, I1, rho, stats);
               Algo_CFSls::calc(step, diagIp, diagI1, op1II, op2II, bs, spinfactor, cs, Ip, I1, rho, stats);
             }
   string name() override { return "CFS"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rho"; }
};

// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).

// Based on the implementation by Markus Greger.
void Algo_CFSls::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const auto sign = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
  const auto &rhoNIp = rho.at(Ip);
  const auto &rhoNI1 = rho.at(I1);
  // Convention: k-loops over retained states, l-loop over discarded states.
  // i-term, Eq. (11).
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      for (const auto rp: diagIp.kept()) {
        const auto E1 = diagI1.value_zero(r1);
        const auto Ep = diagIp.value_zero(rp);
        const auto weight = (spinfactor / stats.Zft) * conj_me(op1II(r1, rp)) * op2II(r1, rp) * exp(-E1 * step.scT()) * (-sign);
        cs->add(step.scale() * (E1-Ep), weight);
      }
    }
  } else {
    // iii-term, Eq. (16), positive frequency excitations
    if (op2II.size1() && rhoNIp.size1()) {
      Matrix op2II_m_rho;
      const ublas::matrix_range<const Matrix> op2II_TK(op2II, ublas::range(0, op2II.size1()), ublas::range(0, rhoNIp.size1()));
      op2II_m_rho = Matrix(op2II_TK.size1(), rhoNIp.size2());
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2II_TK, rhoNIp, 0.0, op2II_m_rho); // rhoNEW <- rhoNEW + factor T U
      for (const auto rl: diagI1.discarded()) {
        for (const auto rk: diagIp.kept()) {
          const auto El       = diagI1.value_zero(rl);
          const auto Ek       = diagIp.value_zero(rk);
          const auto weight = spinfactor * conj_me(op1II(rl, rk)) * op2II_m_rho(rl, rk) * (-sign);
          cs->add(step.scale() * (El-Ek), weight);
        }
      }
    }
  }
}

void Algo_CFSgt::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const auto &rhoNIp = rho.at(Ip);
  const auto &rhoNI1 = rho.at(I1);
  // Convention: k-loops over retained states, l-loop over discarded states.
  // i-term, Eq. (11).
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      for (const auto rp: diagIp.kept()) {
        const auto E1 = diagI1.value_zero(r1);
        const auto Ep = diagIp.value_zero(rp);
        const auto weight = (spinfactor / stats.Zft) * conj_me(op1II(r1, rp)) * op2II(r1, rp) * exp(-Ep * step.scT());
        cs->add(step.scale() * (E1-Ep), weight);
      }
    }
  } else {
    if (rhoNI1.size1() && op1II.size2()) {
      const ublas::matrix_range<const Matrix> op1II_KT(op1II, ublas::range(0, rhoNI1.size1()), ublas::range(0, op1II.size2()));
      Matrix op1II_m_rho(rhoNI1.size2(), op1II_KT.size2());
      if constexpr (std::is_same_v<t_matel, double>) {
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, rhoNI1, op1II_KT, 0.0, op1II_m_rho); // rhoNEW <- rhoNEW + factor T U
      } else {
        const ublas::matrix<std::complex<double>> conj_op1II_KT = conj(op1II_KT);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, rhoNI1, conj_op1II_KT, 0.0, op1II_m_rho); // rhoNEW <- rhoNEW + factor T U
      }
      for (const auto rk: diagI1.kept()) {                                          // ii-term, Eq. (15), negative frequency excitations
        for (const auto rl: diagIp.discarded()) {
          const auto Ek       = diagI1.value_zero(rk);
          const auto El       = diagIp.value_zero(rl);
          const auto weight = spinfactor * op1II_m_rho(rk, rl) * op2II(rk, rl);
          cs->add(step.scale() * (Ek-El), weight);
        }
      }
    }
  }
}

