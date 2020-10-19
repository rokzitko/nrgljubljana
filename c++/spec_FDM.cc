// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

#define LOOP_D(n)                                                                                                                                    \
  for (size_t n = ret##n; n < all##n; n++) {                                                                                                         \
    const t_eigen E##n = diagI##n.absenergyG(n);
#define LOOP_K(n)                                                                                                                                    \
  for (size_t n = 0; n < ret##n; n++) {                                                                                                              \
    const t_eigen E##n = diagI##n.absenergyG(n);

class Algo_FDMls : virtual public Algo {
 protected:
   SpectrumRealFreq spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   explicit Algo_FDMls(SpectrumRealFreq spec, gf_type gt, const Params &P, const bool save = true) 
     : Algo(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) override
   {
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
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T);
         cb->add(energy, weight);
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [K=retj, D=alli]
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, retj), ublas::range(0, alli));
       Matrix rho_op2(retj, alli);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
       LOOP_D(i)
         LOOP_K(j)
           const auto energy = Ej - Ei;
           const auto weight = factor * conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
           cb->add(energy, weight);
         }
       }
     }
     LOOP_K(i)
       LOOP_D(j)
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T);
         cb->add(energy, weight);
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_FDMls() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

class Algo_FDMgt : virtual public Algo {
 protected:
   SpectrumRealFreq spec;
   const int sign; // 1 for bosons, -1 for fermions
   using CB = ChainBinning;
   std::unique_ptr<CB> cb;
   const bool save;
 public:
   explicit Algo_FDMgt(SpectrumRealFreq spec, gf_type gt, const Params &P, const bool save = true) 
     : Algo(P), spec(spec), sign(gf_sign(gt)), save(save) {}
   void begin(const Step &) override { cb = std::make_unique<CB>(P); }
   void calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) override
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
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
         cb->add(energy, weight);
       }
     }
     LOOP_D(i)
       LOOP_K(j)
         const auto energy = Ej - Ei;
         const auto weight = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
         cb->add(energy, weight);
       }
     }
     if (allj > 0 && reti > 0) {
       // B [D=allj,K=reti] x rho [reti,reti]
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, allj), ublas::range(0, reti));
       Matrix op2_rho(allj, reti);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
       LOOP_K(i)
         LOOP_D(j)
           const auto energy = Ej - Ei;
           const auto weight = factor * conj_me(op1(j, i)) * op2_rho(j, i);
           cb->add(energy, weight);
         }
       }
     }
   }
   void end(const Step &step) override {
     spec.mergeCFS(*cb.get());
     cb.reset();
   }
   ~Algo_FDMgt() { if (save) spec.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

class Algo_FDM : public Algo_FDMls, public Algo_FDMgt {
 private:
   SpectrumRealFreq spec_tot;
 public:
   explicit Algo_FDM(SpectrumRealFreq spec, gf_type gt, const Params &P) :
     Algo(P), Algo_FDMls(spec, gt, P, false), Algo_FDMgt(spec, gt, P, false), spec_tot(spec) {}
   void begin(const Step &step) override {
     Algo_FDMgt::begin(step);
     Algo_FDMls::begin(step);
   }
   void calc(const Step &step, const Eigen &diagIp, const Eigen &diagI1, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ip, const Invar &I1, const DensMatElements &rho, const Stats &stats) override
   {
     Algo_FDMgt::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
     Algo_FDMls::calc(step, diagIp, diagI1, op1, op2, factor, Ip, I1, rho, stats);
   }
   void end(const Step &step) override {
     spec_tot.mergeCFS(*Algo_FDMgt::cb.get());
     spec_tot.mergeCFS(*Algo_FDMls::cb.get());
     Algo_FDMgt::cb.reset();
     Algo_FDMls::cb.reset();
   }
   ~Algo_FDM() { spec_tot.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};

class Algo_FDMmats : public Algo {
 private:
   GFMatsubara gf;
   const int sign;
   gf_type gt;
   using CM = ChainMatsubara;
   std::unique_ptr<CM> cm;
 public:
   explicit Algo_FDMmats(GFMatsubara gf, gf_type gt, const Params &P) : Algo(P), gf(gf), sign(gf_sign(gt)), gt(gt) {}
   void begin(const Step &) override { cm = std::make_unique<CM>(P, gt); }
   void calc(const Step &step, const Eigen &diagIi, const Eigen &diagIj, const Matrix &op1, const Matrix &op2, t_coef factor,
             const Invar &Ii, const Invar &Ij, const DensMatElements &rhoFDM, const Stats &stats) override
   {
     const size_t cutoff = P.mats;
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
         const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T); // a[ij] b[ji] exp(-beta e[i])
         const auto weightB = factor * conj_me(op1(j, i)) * op2(j, i) * (-sign) * wnf * exp(-Ej / P.T); // a[ij] b[ji] sign exp(-beta e[j])
         #pragma omp parallel for schedule(static)
         for (size_t n = 1; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
         if (gt == gf_type::fermionic || abs(energy) > WEIGHT_TOL)
           cm->add(size_t(0), (weightA + weightB) / (cmpl(0, ww(0, gt, P.T)) - energy));
         else // bosonic w=0 && Ei=Ej case
           cm->add(size_t(0), (-weightA / t_weight(P.T)));
       }
     }
     if (retj > 0 && alli > 0) {
       // rho [retj, retj] x B [retj, alli]
       const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, retj), ublas::range(0, alli));
       Matrix rho_op2(retj, alli);
       atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
       LOOP_D(i)
         LOOP_K(j)
           const auto energy = Ej - Ei;
           const auto weightA = factor * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ei / P.T);
           const auto weightB = factor * conj_me(op1(j, i)) * rho_op2(j, i) * (-sign);
           #pragma omp parallel for schedule(static)
           for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
         }
       }
     }
     if (allj > 0 && reti > 0) {
        // B [allj,reti] x rho [reti,reti]
        const ublas::matrix_range<const Matrix> op2cut(op2, ublas::range(0, allj), ublas::range(0, reti));
        Matrix op2_rho(allj, reti);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
        LOOP_K(i)
          LOOP_D(j)
            const auto energy = Ej - Ei;
            const auto weightA = factor * conj_me(op1(j, i)) * op2_rho(j, i);
            const auto weightB = factor * (-sign) * conj_me(op1(j, i)) * op2(j, i) * wnf * exp(-Ej / P.T);
            #pragma omp parallel for schedule(static)
            for (size_t n = 0; n < cutoff; n++) cm->add(n, (weightA + weightB) / (cmpl(0, ww(n, gt, P.T)) - energy));
          }
        }
     }
   }
   void end(const Step &step) override {
     gf.merge(*cm.get());
     cm.reset();
   }
   ~Algo_FDMmats() { gf.save(); }
   std::string rho_type() override { return "rhoFDM"; }
};
