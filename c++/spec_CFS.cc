// Choose one of the following two! OLD is the non-optimized code by Rok Zitko, OPTIMIZED is the hand-tuned code
// contributed by Markus Greger. The optimized code is faster by an order of magnitude!

class Algo_CFSls : virtual public Algo {
 public:
   Algo_CFSls(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "CFSls"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rho"; }
};

class Algo_CFSgt : virtual public Algo {
 public:
   Algo_CFSgt(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "CFSgt"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rho"; }
};

class Algo_CFS : public Algo_CFSls, public Algo_CFSgt {
 public:
   Algo_CFS(const Params &P) : Algo(P), Algo_CFSls(P), Algo_CFSgt(P) {}
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

//#define Algo_CFS_OLD
#define Algo_CFS_OPTIMIZED

// Cf. Peters, Pruschke, Anders, Phys. Rev. B 74, 245113 (2006).

#if defined(NRG_COMPLEX) || defined(Algo_CFS_OLD)
void Algo_CFSls::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  // Convention: k-loops over retained states, l-loop over discarded states.
  // i-term, Eq. (11). This part is analogous to that for Algo_FT, i.e., it has the form of the usual Lehmann
  // representation.
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      const t_eigen E1 = diagI1.value_zero(r1);
      for (const auto rp: diagIp.kept()) {
        const t_eigen Ep = diagIp.value_zero(rp);
        DELTA d;
        d.energy = E1 - Ep;
        d.weight = (spinfactor / stats.Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp) * exp(-E1 * step.scT()) * (-sign);
        cs->add(step.scale() * d.energy, d.weight);
      }
    }
  } else {
    // iii-term, Eq. (16), positive frequency excitations
    const auto dimA = diagI1.getnrstored();
    for (const auto rl: diagI1.discarded()) {
      const t_eigen El = diagI1.value_zero(rl);
      for (const auto rk: diagIp.kept()) {
        const t_eigen Ek = diagIp.value_zero(rk);
        DELTA d;
        d.energy = El - Ek;
        my_assert(d.energy >= 0.0); // always positive!
        weight_bucket sum;
        for (const auto rkp: diagIp.kept()) sum += op2II(rl, rkp) * rhoNIp(rkp, rk);
        d.weight = spinfactor * CONJ_ME(op1II(rl, rk)) * t_weight(sum) * (-sign);
        cs->add(step.scale() * d.energy, d.weight);
      }
    }
  } // if (last)
}

void Algo_CFSgt::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  // Convention: k-loops over retained states, l-loop over discarded states.
  // i-term, Eq. (11).
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      const t_eigen E1 = diagI1.value_zero(r1);
      for (const auto rp: diagIp.kept()) {
        const t_eigen Ep = diagIp.value_zero(rp);
        DELTA d;
        d.energy = E1 - Ep;
        d.weight = (spinfactor / stats.Zft) * CONJ_ME(op1II(r1, rp)) * op2II(r1, rp) * exp(-Ep * step.scT());
        cs->add(step.scale() * d.energy, d.weight);
      }
    }
  } else {
    // ii-term, Eq. (15), negative frequency excitations
    for (const auto rk: diagI1.kept()) {
      const t_eigen Ek  = diagI1.value_zero(rk);
      const auto dimB = diagIp.getnrstored();
      for (const auto rl: diagIp.discarded()) {
        const t_eigen El = diagIp.value_zero(rl);
        DELTA d;
        d.energy = Ek - El;
        my_assert(d.energy <= 0.0); // always negative!
        weight_bucket sum;
        for (const auto rkp: diagI1.kept()) sum += CONJ_ME(op1II(rkp, rl)) * rhoNI1(rkp, rk);
        d.weight = spinfactor * t_weight(sum) * op2II(rk, rl);
        cs->add(step.scale() * d.energy, d.weight);
      }
    }
  } // if (last)
}
#endif

// Based on the implementation by Markus Greger.
#if defined(NRG_REAL) && defined(Algo_CFS_OPTIMIZED)
void Algo_CFSls::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, double spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  auto dimp            = rhoNIp.size1();
  // Convention: k-loops over retained states, l-loop over discarded
  // states.
  // i-term, Eq. (11).
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      const double E1 = diagI1.value_zero(r1);
      for (const auto rp: diagIp.kept()) {
        const double Ep = diagIp.value_zero(rp);
        double d_energy = E1 - Ep;
        double d_weight = (spinfactor / stats.Zft) * op1II(r1, rp) * op2II(r1, rp) * exp(-E1 * step.scT()) * (-sign);
        cs->add(step.scale() * d_energy, d_weight);
      }
    }
  } else {
    // iii-term, Eq. (16), positive frequency excitations
    const auto dimA     = diagI1.getnrstored();
    if (dimA && dimp) {
      Matrix op2II_m_rho;
      const ublas::matrix_range<const Matrix> op2II_TK(op2II, ublas::range(0, op2II.size1()), ublas::range(0, rhoNIp.size1()));
      op2II_m_rho = Matrix(op2II_TK.size1(), rhoNIp.size2());
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2II_TK, rhoNIp, 0.0, op2II_m_rho); // rhoNEW <- rhoNEW + factor T U
      for (const auto rl: diagI1.discarded()) {
        const double El = diagI1.value_zero(rl);
        for (const auto rk: diagIp.kept()) {
          const double Ek       = diagIp.value_zero(rk);
          const double d_energy = El - Ek;
          const double sum      = op2II_m_rho(rl, rk);
          const double d_weight = spinfactor * op1II(rl, rk) * sum * (-sign);
          cs->add(step.scale() * d_energy, d_weight);
        }
      }
    }
  } // if (last)
}

void Algo_CFSgt::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, double spinfactor,
                      spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const {
  const Matrix &rhoNIp = rho.at(Ip);
  const Matrix &rhoNI1 = rho.at(I1);
  auto dim1            = rhoNI1.size1();
  // Convention: k-loops over retained states, l-loop over discarded states.
  // i-term, Eq. (11).
  if (step.last()) {
    for (const auto r1: diagI1.kept()) {
      const double E1 = diagI1.value_zero(r1);
      for (const auto rp: diagIp.kept()) {
        const double Ep = diagIp.value_zero(rp);
        double d_energy = E1 - Ep;
        double d_weight = (spinfactor / stats.Zft) * op1II(r1, rp) * op2II(r1, rp) * exp(-Ep * step.scT());
        cs->add(step.scale() * d_energy, d_weight);
      }
    }
  } else {
    const auto dimB     = diagIp.getnrstored();
    if (dim1 && dimB) {
      const ublas::matrix_range<const Matrix> op1II_KT(op1II, ublas::range(0, rhoNI1.size1()), ublas::range(0, op1II.size2()));
      Matrix op1II_m_rho(rhoNI1.size2(), op1II_KT.size2());
      atlas::gemm(CblasTrans, CblasNoTrans, 1.0, rhoNI1, op1II_KT, 0.0, op1II_m_rho); // rhoNEW <- rhoNEW + factor T U
      for (const auto rk: diagI1.kept()) {                                          // ii-term, Eq. (15), negative frequency excitations
        const double Ek = diagI1.value_zero(rk);
        for (const auto rl: diagIp.discarded()) {
          const double El       = diagIp.value_zero(rl);
          const double d_energy = Ek - El;
          const double sum      = op1II_m_rho(rk, rl);
          double d_weight       = spinfactor * sum;
          d_weight *= op2II(rk, rl);
          cs->add(step.scale() * d_energy, d_weight);
        }
      }
    }
  } // if (last)
}
#endif
