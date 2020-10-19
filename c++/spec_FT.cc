class Algo_FT : public Algo {
 private:
   SpectrumRealFreq spec;
   const int sign;
   using CS = ChainBinning;
   std::unique_ptr<CS> cs;
 public:
   explicit Algo_FT(SpectrumRealFreq spec, const int sign, const Params &P) : Algo(P), spec(std::move(spec)), sign(sign) {}
   void begin(const Step &) override {
     cs = std::make_unique<CS>(P);
   }
   // The first matrix element is conjugated!
   // This is <rp|OP1^dag|r1> <r1|OP2|rp> (wp - s*w1)/(z+Ep-E1)
   // s=1 for bosons, s=-1 for fermions
   // See Eq.(9) in Peters, Pruschke, Anders, PRB (74) 245114 (2006)
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, const t_coef factor,
             const Invar &, const Invar &, const DensMatElements &, const Stats &stats) override
   {
     for (const auto r1: diagI1.kept()) {
       const auto E1 = diagI1.value_zero(r1);
       for (const auto rp: diagIp.kept()) {
         const auto Ep = diagIp.value_zero(rp);
         const auto weight = (factor / stats.Zft) * conj_me(op1(r1, rp)) * op2(r1, rp) * ((-sign) * exp(-E1 * step.scT()) + exp(-Ep * step.scT()));
         const auto energy = E1 - Ep;
         cs->add(step.scale() * energy, weight);
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeNN2(*cs.get(), step);
     cs.reset();
   }
};

/*
 class Algo_FTmats : public Algo {
 private:
   SpectrumMatsubara spec;
 public:
   explicit Algo_FTmats(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &bs) override { return make_shared<ChainSpectrumMatsubara>(P, bs.mt); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_coef, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "FTmats"; }
};

void Algo_FTmats::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs,
                       t_coef spinfactor, spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &, const Stats &stats) const {
  const size_t cutoff = P.mats;
  const auto sign     = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
  auto csm            = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
  for (const auto r1: diagI1.kept()) {
    const auto E1 = diagI1.value_zero(r1);
    for (const auto rp: diagIp.kept()) {
      const auto Ep = diagIp.value_zero(rp);
      const auto weight = (spinfactor / stats.Zft) * conj_me(op1II(r1, rp)) * op2II(r1, rp) * ((-sign) * exp(-E1 * step.scT()) + exp(-Ep * step.scT())); // sign!
      const auto energy = E1 - Ep;
#pragma omp parallel for schedule(static)
      for (size_t n = 1; n < cutoff; n++) csm->add(n, weight / (cmpl(0, ww(n, bs.mt, P.T)) - step.scale() * energy));
      if (abs(energy) > WEIGHT_TOL || bs.mt == matstype::fermionic)
        csm->add(size_t(0), weight / (cmpl(0, ww(0, bs.mt, P.T)) - step.scale() * energy));
      else // bosonic w=0 && E1=Ep case
        csm->add(size_t(0), (spinfactor / stats.Zft) * conj_me(op1II(r1, rp)) * op2II(r1, rp) * (-exp(-E1 * step.scT()) / P.T));
    }
  }
}

class Algo_GT_generic : public Algo {
 protected:
   int power{};
 public:
   explicit Algo_GT_generic(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumTemp>(P); }
   void calc(const Step &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_coef, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "ERROR"; }
};

class Algo_GT : public Algo_GT_generic {
 public:
   explicit Algo_GT(const Params &P) : Algo_GT_generic(P) { power = 0; }
   string name() override { return "GT"; }
};

class Algo_I1T : public Algo_GT_generic {
 public:
   explicit Algo_I1T(const Params &P) : Algo_GT_generic(P) { power = 1; }
  string name() override { return "I1T"; }
};

class Algo_I2T : public Algo_GT_generic {
 public:
   explicit Algo_I2T(const Params &P) : Algo_GT_generic(P) { power = 2; }
   string name() override { return "I2T"; }
};

// Calculation of the temperature-dependent linear conductrance G(T) using
// the linear response theory & impurity-level spectral density. Note that
// we are not calculating a spectral function (i.e. a collection of delta
// peaks), but rather a tabulated G(T), so binning needs to be turned off.
// See Yoshida, Seridonio, Oliveira, arxiv:0906.4289, Eq. (8).
void Algo_GT_generic::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs,
                           t_coef spinfactor, spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &, const Stats &stats) const {
  my_assert(bs.mt == matstype::fermionic);                  // restricted implementation
  const double temperature = P.gtp * step.scale(); // in absolute units!
  const double beta        = 1.0 / temperature;
  weight_bucket value;
  for (const auto r1: diagI1.kept()) {
    const t_eigen E1 = diagI1.value_zero(r1);
    for (const auto rp: diagIp.kept()) {
      const t_eigen Ep = diagIp.value_zero(rp);
      // Note that Zgt needs to be calculated with the same
      // 'temperature' parameter that we use for the exponential
      // functions in the following equation.
      value += beta * (spinfactor / stats.Zgt) * conj_me(op1II(r1, rp)) * op2II(r1, rp)
         / (exp(+E1 * step.scale() * beta) + exp(+Ep * step.scale() * beta)) * pow((E1 - Ep) * step.scale(), power);
    } // loop over r1
  }   // loop over rp
  cs->add(temperature, value);
}

// weight=(exp(-beta Em)-exp(-beta En))/(beta En-beta Em).
// NOTE: arguments En, Em are order omega_N, while beta is order
// 1/omega_N, thus the combinations betaEn and betaEm are order 1.
// Also En>0, Em>0, since these are excitation energies !
inline t_weight chit_weight(double En, double Em, double beta) {
  const double betaEn = beta * En;
  const double betaEm = beta * Em;
  const double x      = betaEn - betaEm;
  if (abs(x) > WEIGHT_TOL) {
    // If one of {betaEm,betaEn} is small, one of exp() will have a
    // value around 1, the other around 0, thus the overall result
    // will be approximately +-1/x.
    return (exp(-betaEm) - exp(-betaEn)) / x;
  } else {
    // Special case for Em~En. In this case, we are integrating a
    // constant over tau\in{0,\beta}, and dividing this by beta we
    // get 1. What remains is the Boltzmann weight exp(-betaEm).
    return exp(-betaEm);
  }
}

class Algo_CHIT : public Algo {
 public:
   explicit Algo_CHIT(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumTemp>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_coef, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "CHIT"; }
};

// Calculation of the temperature-dependent susceptibility chi_AB(T) using the linear response theory and the matrix
// elements of global operators. Binning needs to be turned off. Note that Zchit needs to be calculated with the same
// 'temperature' parameter that we use for the exponential functions in the following equation. The output is
// chi/beta = k_B T chi, as we prefer.
void Algo_CHIT::calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_coef spinfactor,
                     spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &, const Stats &stats) const {
  my_assert(bs.mt == matstype::bosonic); // restricted implementation
  const double temperature = P.chitp * step.scale(); // in absolute units!
  const double beta        = 1.0 / temperature;
  weight_bucket w;
  for (const auto r1: diagI1.kept()) {
    const t_eigen E1 = diagI1.value_zero(r1);
    for (const auto rp: diagIp.kept()) {
      const t_eigen Ep      = diagIp.value_zero(rp);
      const t_weight weight = chit_weight(step.scale() * E1, step.scale() * Ep, beta);
      w += op1II(r1, rp) * op2II(rp, r1) * weight;
    }
  }
  cs->add(temperature, spinfactor / stats.Zchit * t_weight(w));
}
*/
