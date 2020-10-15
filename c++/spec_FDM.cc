// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

class Algo_FDMls : virtual public Algo {
 public:
   explicit Algo_FDMls(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "FDMls"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rhoFDM"; }
};

class Algo_FDMgt : virtual public Algo {
 public:
   explicit Algo_FDMgt(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "FDMgt"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rhoFDM"; }
};

class Algo_FDM : public Algo_FDMls, public Algo_FDMgt {
 public:
   explicit Algo_FDM(const Params &P) : Algo(P), Algo_FDMls(P), Algo_FDMgt(P) {}
   spCS_t make_cs(const BaseSpectrum &) override { return make_shared<ChainSpectrumBinning>(P); }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
             spCS_t cs, const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) const override {
               Algo_FDMgt::calc(step, diagIp, diagI1, op1II, op2II, bs, spinfactor, cs, Ip, I1, rho, stats);
               Algo_FDMls::calc(step, diagIp, diagI1, op1II, op2II, bs, spinfactor, cs, Ip, I1, rho, stats);
             }
   string name() override { return "FDM"; }
   string merge() override { return "CFS"; }
   string rho_type() override { return "rhoFDM"; }
};

#define LOOP_D(n)                                                                                                                                    \
  for (size_t n = ret##n; n < all##n; n++) {                                                                                                         \
    const t_eigen E##n = diagI##n.absenergyG(n);
#define LOOP_K(n)                                                                                                                                    \
  for (size_t n = 0; n < ret##n; n++) {                                                                                                              \
    const t_eigen E##n = diagI##n.absenergyG(n);

// *********** Greater correlation function ***********
void Algo_FDMgt::calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) const 
{
  const double wnf   = stats.wnfactor[step.ndx()];
  const Matrix &rhoi = rhoFDM.at(Ii);
  const Matrix &rhoj = rhoFDM.at(Ij);
  const auto reti  = step.last() ? 0 : diagIi.getnrkept();
  const auto retj  = step.last() ? 0 : diagIj.getnrkept();
  const auto alli  = diagIi.getnrstored();
  const auto allj  = diagIj.getnrstored();
  LOOP_D(i)
   LOOP_D(j)
    const auto energy = Ej - Ei;
    const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei / P.T);
    cs->add(energy, weight);
   }
  }
  LOOP_D(i)
   LOOP_K(j)
    const auto energy = Ej - Ei;
    const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei / P.T);
    cs->add(energy, weight);
   }
  }
  if (allj > 0 && reti > 0) {
    // B [D=allj,K=reti] x rho [reti,reti]
    const ublas::matrix_range<const Matrix> op2cut(op2II, ublas::range(0, allj), ublas::range(0, reti));
    Matrix op2_rho(allj, reti);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
    LOOP_K(i)
     LOOP_D(j)
      const auto energy = Ej - Ei;
      const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * op2_rho(j, i);
      cs->add(energy, weight);
     }
    }
  }
}

// ************ Lesser correlation functions ***************
void Algo_FDMls::calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs, t_factor spinfactor,
                      spCS_t cs, const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) const {
  const auto sign    = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
  const double wnf   = stats.wnfactor[step.ndx()];
  const Matrix &rhoi = rhoFDM.at(Ii);
  const Matrix &rhoj = rhoFDM.at(Ij);
  const auto reti  = step.last() ? 0 : diagIi.getnrkept();
  const auto retj  = step.last() ? 0 : diagIj.getnrkept();
  const auto alli  = diagIi.getnrstored();
  const auto allj  = diagIj.getnrstored();
  LOOP_D(i)
   LOOP_D(j) // B3
    const auto energy = Ej - Ei;
    const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej / P.T);
    cs->add(energy, weight);
   }
  }
  if (retj > 0 && alli > 0) {
    // rho [retj, retj] x B [K=retj, D=alli]
    const ublas::matrix_range<const Matrix> op2cut(op2II, ublas::range(0, retj), ublas::range(0, alli));
    Matrix rho_op2(retj, alli);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
    LOOP_D(i)
     LOOP_K(j)
      const auto energy = Ej - Ei;
      const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * rho_op2(j, i) * (-sign);
      cs->add(energy, weight);
     }
    }
  }
  LOOP_K(i)
   LOOP_D(j)
    const auto energy = Ej - Ei;
    const auto weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej / P.T);
    cs->add(energy, weight);
   }
  }
}

class Algo_FDMmats : public Algo {
 public:
   Algo_FDMmats(const Params &P) : Algo(P) {}
   spCS_t make_cs(const BaseSpectrum &bs) override { return make_shared<ChainSpectrumMatsubara>(P, bs.mt); }
   void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
             const Invar &, const DensMatElements &, const Stats &stats) const override;
   string name() override { return "FDMmats"; }
   string rho_type() override { return "rhoFDM"; }
};

// *********** Matsubara axis version  ***********

void Algo_FDMmats::calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1II, const Matrix &op2II, const BaseSpectrum &bs,
                        t_factor spinfactor, spCS_t cs, const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) const {
  const size_t cutoff = P.mats;
  const auto sign    = bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC;
  auto csm           = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
  const double wnf   = stats.wnfactor[step.ndx()];
  const Matrix &rhoi = rhoFDM.at(Ii);
  const Matrix &rhoj = rhoFDM.at(Ij);
  const auto reti  = (step.last() ? 0 : rhoi.size1());
  const auto retj  = (step.last() ? 0 : rhoj.size1());
  const auto alli  = diagIi.getnrstored();
  const auto allj  = diagIj.getnrstored();
  LOOP_D(i)
   LOOP_D(j) 
    const auto energy = Ej - Ei;
    const auto weightA = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei / P.T); // a[ij] b[ji] exp(-beta e[i])
    const auto weightB = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej / P.T); // a[ij] b[ji] sign exp(-beta e[j])
    #pragma omp parallel for schedule(static)
    for (size_t n = 1; n < cutoff; n++) csm->add(n, (weightA + weightB) / (cmpl(0, ww(n, bs.mt, P.T)) - energy));
    if (bs.mt == matstype::fermionic || abs(energy) > WEIGHT_TOL)
      csm->add(size_t(0), (weightA + weightB) / (cmpl(0, ww(0, bs.mt, P.T)) - energy));
    else // bosonic w=0 && Ei=Ej case
      csm->add(size_t(0), (-weightA / t_weight(P.T)));
    }
  }
  if (retj > 0 && alli > 0) {
     // rho [retj, retj] x B [retj, alli]
    const ublas::matrix_range<const Matrix> op2cut(op2II, ublas::range(0, retj), ublas::range(0, alli));
    Matrix rho_op2(retj, alli);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
    LOOP_D(i) 
     LOOP_K(j) 
      const auto energy = Ej - Ei;
      const auto weightA = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei / P.T);
      const auto weightB = spinfactor * CONJ_ME(op1II(j, i)) * rho_op2(j, i) * (-sign);
      #pragma omp parallel for schedule(static)
      for (size_t n = 0; n < cutoff; n++) csm->add(n, (weightA + weightB) / (cmpl(0, ww(n, bs.mt, P.T)) - energy));
     }
    }
  }
  if (allj > 0 && reti > 0) {
    // B [allj,reti] x rho [reti,reti]
    const ublas::matrix_range<const Matrix> op2cut(op2II, ublas::range(0, allj), ublas::range(0, reti));
    Matrix op2_rho(allj, reti);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
    LOOP_K(i) 
     LOOP_D(j)
      const auto energy = Ej - Ei;
      const auto weightA = spinfactor * CONJ_ME(op1II(j, i)) * op2_rho(j, i);
      const auto weightB = spinfactor * (-sign) * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ej / P.T);
      #pragma omp parallel for schedule(static)
      for (size_t n = 0; n < cutoff; n++) csm->add(n, (weightA + weightB) / (cmpl(0, ww(n, bs.mt, P.T)) - energy));
     }
    }
  }
}
