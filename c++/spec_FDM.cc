// Recall: II=(Ij,Ii) <i|A|j> <j|B|i>. B is d^dag. We conjugate A.

class SPEC_FDMls : virtual public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "FDMls"; }
   string merge() { return "CFS"; }
};

class SPEC_FDMgt : virtual public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "FDMgt"; }
   string merge() { return "CFS"; }
};

class SPEC_FDM : public SPEC_FDMls, public SPEC_FDMgt
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &) { return new ChainSpectrumBinning; }
   void calc(const Eigen & a1, const Eigen & a2, const Matrix & a3, const Matrix & a4,
	     const BaseSpectrum & a5, t_factor a6, ChainSpectrum *a7, const Invar & a8, const Invar & a9)
   {
      SPEC_FDMgt::calc(a1,a2,a3,a4,a5,a6,a7,a8,a9);
      SPEC_FDMls::calc(a1,a2,a3,a4,a5,a6,a7,a8,a9);
   }
   string name() { return "FDM"; }
   string merge() { return "CFS"; }
};

#define LOOP_D(n) for (size_t n = ret##n; n < all##n; n++) { const t_eigen E##n = diagI##n.absenergy(n);
#define LOOP_K(n) for (size_t n = 0; n < ret##n; n++) { const t_eigen E##n = diagI##n.absenergy(n);

// *********** Greater correlation function ***********
void SPEC_FDMgt::calc(const Eigen & diagIi, const Eigen & diagIj,
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs, t_factor spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ii, const Invar &Ij)
{
   const double wnf = STAT::wnfactor[STAT::N];
   const Matrix & rhoi = rhoFDM[Ii];
   const Matrix & rhoj = rhoFDM[Ij];
   const size_t reti = (LAST_ITERATION() ? 0 : rhoi.size1());
   const size_t retj = (LAST_ITERATION() ? 0 : rhoj.size1());
   const size_t alli = diagIi.getnr();
   const size_t allj = diagIj.getnr();
   LOOP_D(i) LOOP_D(j) // A3
      DELTA d;
      d.energy = Ej-Ei;
      d.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei/P::T);
      cs->add(d.energy, d.weight);
   }}
   LOOP_D(i) LOOP_K(j) // A2
      DELTA d;
      d.energy = Ej-Ei;
      d.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei/P::T);
      cs->add(d.energy, d.weight);
   }}
   if (allj>0 && reti>0) {
      // B [D=allj,K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op2cut(op2II, range(0, allj), range(0, reti));
      Matrix op2_rho(allj, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
      LOOP_K(i) LOOP_D(j) // A1
	 DELTA d;
         d.energy = Ej-Ei;
         d.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2_rho(j, i);
         cs->add(d.energy, d.weight);
      }}
   }
}

// ************ Lesser correlation functions ***************
void SPEC_FDMls::calc(const Eigen & diagIi, const Eigen & diagIj,
		      const Matrix & op1II, const Matrix & op2II,
		      const BaseSpectrum &bs, t_factor spinfactor,
		      ChainSpectrum * cs,
		      const Invar &Ii, const Invar &Ij)
{
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   const double wnf = STAT::wnfactor[STAT::N];
   const Matrix & rhoi = rhoFDM[Ii];
   const Matrix & rhoj = rhoFDM[Ij];
   const size_t reti = (LAST_ITERATION() ? 0 : rhoi.size1());
   const size_t retj = (LAST_ITERATION() ? 0 : rhoj.size1());
   const size_t alli = diagIi.getnr();
   const size_t allj = diagIj.getnr();
   LOOP_D(i) LOOP_D(j) // B3
      DELTA d;
      d.energy = Ej-Ei;
      d.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej/P::T);
      cs->add(d.energy, d.weight);
   }}
   if (retj>0 && alli>0) {
      // rho [retj, retj] x B [K=retj, D=alli]
      const matrix_range<const Matrix> op2cut(op2II, range(0, retj), range(0, alli));
      Matrix rho_op2(retj, alli);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
      LOOP_D(i) LOOP_K(j) // B2
        DELTA d;
        d.energy = Ej-Ei;
        d.weight = spinfactor * CONJ_ME(op1II(j, i)) * rho_op2(j, i) * (-sign);
        cs->add(d.energy, d.weight);
      }}
   }
   LOOP_K(i) LOOP_D(j) // B1
     DELTA d;
     d.energy = Ej-Ei;
     d.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej/P::T);
     cs->add(d.energy, d.weight);
   }}
}

class SPEC_FDMmats : public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &bs) { return new ChainSpectrumMatsubara(bs.mt); }
   void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &,
	     const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &);
   string name() { return "FDMmats"; }
};

// *********** Matsubara axis version  ***********

void SPEC_FDMmats::calc(const Eigen & diagIi, const Eigen & diagIj,
			const Matrix & op1II, const Matrix & op2II,
			const BaseSpectrum &bs, t_factor spinfactor,
			ChainSpectrum * cs,
			const Invar &Ii, const Invar &Ij)
{
   const size_t cutoff = P::mats;
   // (-sign)=1 for fermionic case, (-sign)=-1 for bosonic case
   double sign = (bs.mt == matstype::bosonic ? S_BOSONIC : S_FERMIONIC);
   auto *csm = dynamic_cast<ChainSpectrumMatsubara*>(cs);
   const double wnf = STAT::wnfactor[STAT::N];
   const Matrix & rhoi = rhoFDM[Ii];
   const Matrix & rhoj = rhoFDM[Ij];
   const size_t reti = (LAST_ITERATION() ? 0 : rhoi.size1());
   const size_t retj = (LAST_ITERATION() ? 0 : rhoj.size1());
   const size_t alli = diagIi.getnr();
   const size_t allj = diagIj.getnr();
   LOOP_D(i) LOOP_D(j)
      DELTA dA; // A3
      dA.energy = Ej-Ei;
      dA.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei/P::T); // a[ij] b[ji] exp(-beta e[i])
      DELTA dB; // B3
      dB.energy = dA.energy;
      dB.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * (-sign) * wnf * exp(-Ej/P::T); // a[ij] b[ji] sign exp(-beta e[j])
#pragma omp parallel for schedule(static)
      for (size_t n = 1; n < cutoff; n++)
        csm->add(n, (dA.weight+dB.weight)/(cmpl(0,w(n,bs.mt))-dA.energy));
      if (bs.mt == matstype::fermionic || abs(dA.energy) > WEIGHT_TOL)
        csm->add(size_t(0), (dA.weight+dB.weight)/(cmpl(0,w(0,bs.mt))-dA.energy));
       else // bosonic w=0 && Ei=Ej case
        csm->add(size_t(0), (-dA.weight/t_weight(P::T)));
   }}
   if (retj>0 && alli>0) {
      // rho [retj, retj] x B [retj, alli]
      const matrix_range<const Matrix> op2cut(op2II, range(0, retj), range(0, alli));
      Matrix rho_op2(retj, alli);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
      LOOP_D(i) LOOP_K(j)
	 DELTA dA; // A2
         dA.energy = Ej-Ei;
	 dA.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ei/P::T);
	 DELTA dB; // B2
	 dB.energy = dA.energy;
	 dB.weight = spinfactor * CONJ_ME(op1II(j, i)) * rho_op2(j, i) * (-sign);
#pragma omp parallel for schedule(static)
	 for (size_t n = 0; n < cutoff; n++)
	    csm->add(n, (dA.weight+dB.weight)/(cmpl(0,w(n,bs.mt))-dB.energy));
      }}
   }
   if (allj>0 && reti>0) {
      // B [allj,reti] x rho [reti,reti]
      const matrix_range<const Matrix> op2cut(op2II, range(0, allj), range(0, reti));
      Matrix op2_rho(allj, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op2cut, rhoi, 0.0, op2_rho);
      LOOP_K(i) LOOP_D(j)
	 DELTA dA; // A1
         dA.energy = Ej-Ei;
         dA.weight = spinfactor * CONJ_ME(op1II(j, i)) * op2_rho(j, i);
         DELTA dB; // B1
         dB.energy = dA.energy;
         dB.weight = spinfactor * (-sign) * CONJ_ME(op1II(j, i)) * op2II(j, i) * wnf * exp(-Ej/P::T);
#pragma omp parallel for schedule(static)
         for (size_t n = 0; n < cutoff; n++)
	   csm->add(n, (dA.weight+dB.weight)/(cmpl(0,w(n,bs.mt))-dA.energy));
      }}
   }
}

// ************************************************************************************************************

class SPEC_FDM_v3mm : public SPEC
{
public:
   ChainSpectrum * make_cs(const BaseSpectrum &bs) { return new ChainSpectrumMatsubara2(matstype::fb); }
   string name() { return "FDM_v3mm"; }
   void calc_A(const Eigen &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const Matrix &,
	       const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &, const Invar &);
   void calc_B(const Eigen &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const Matrix &,
	       const BaseSpectrum &, t_factor, ChainSpectrum *, const Invar &, const Invar &, const Invar &);
};

template <typename T>
inline std::complex<T> boltz_fnc_noscale (T E1, T E2, T bzE1, T bzE2, T wn, short n)
{
   const T WEIGHT_TOL = 1e-8;
   if (n != 0 || abs(E1-E2) > WEIGHT_TOL) // integer test first!!
      return (bzE1-bzE2)/cmpl(E1-E2,wn); // bzE1 = exp(-E1/T), bzE2 = exp(-E2/T)
   else
      return -bzE1/P::T;
}

// Check list:
// 1. wn vs. 1, OK
// 2. factor /P::T vs *scT, OK
// 3. absolute energies (absenergy), no scale* factors, OK
// 4. rho*exp(Ei/P::T) vs. 1, OK
// 5. KKK correctness, OK
// 6. numerical prefactors
// 7. differences A vs. B: conjugation, index reversal

typedef std::vector<matrix<t_weight>> res_t;

#undef LOOP_D
#undef LOOP_K

#define LOOP_D(n) for (size_t n = ret##n; n < all##n; n++)
#define LOOP_K(n) for (size_t n = 0; n < ret##n; n++)
#define LOOP_X(n) for (size_t n = 0; n < all##n; n++)

#define LOOPVARS \
   const t_eigen Ei = diagi.absenergy(i); \
   const t_eigen Ej = diagj.absenergy(j); \
   const t_eigen El = diagl.absenergy(l); \
   const t_eigen bzEi = diagi.boltzmann(i); \
   const t_eigen bzEj = diagj.boltzmann(j); \
   const t_eigen bzEl = diagl.boltzmann(l);

// n bosonic, m fermionic
#define SUM(factor) \
   if (abs(mat) < v3mmcutoff) continue; \
   for (short n = 0; n < maxn; n++) { \
     for (short m = 0; m < maxm; m++) { \
	const t_weight w = spinfactor * mat * (factor); \
        res[omp_get_thread_num()](n,m) += w; \
     } \
   }


void SPEC_FDM_v3mm::calc_A(const Eigen & diagi, const Eigen & diagj, const Eigen & diagl,
			   const Matrix & op1, const Matrix & op2, const Matrix & op3,
			   const BaseSpectrum &bs, t_factor spinfactor,
			   ChainSpectrum *cs, const Invar & Ii, const Invar & Ij, const Invar & Il)
{
   const double v3mmcutoff = P::v3mmcutoff * sqr(P::T); // order of contributions is prop to 1/T^2
   size_t nr = omp_get_max_threads();
   nrglog('G', "nr=" << nr);
   const short maxn = P::mats;
   const short maxm = P::mats;
   res_t res(nr);
   for (size_t j = 0; j < nr; j++) {
      res[j].resize(maxn, maxm); // (n,m) order
      res[j].clear();
   }
   auto *csm = dynamic_cast<ChainSpectrumMatsubara2*>(cs);
   const double wnf = STAT::wnfactor[STAT::N];
   const Matrix & rhoi = rhoFDM[Ii];
   const Matrix & rhoj = rhoFDM[Ij];
   const Matrix & rhol = rhoFDM[Il];
   const size_t reti = (LAST_ITERATION() ? 0 : rhoi.size1()); // retained states
   const size_t retj = (LAST_ITERATION() ? 0 : rhoj.size1());
   const size_t retl = (LAST_ITERATION() ? 0 : rhol.size1());
   const size_t alli = diagi.getnr(); // all states
   const size_t allj = diagj.getnr();
   const size_t alll = diagl.getnr();
   my_assert(op1.size1() == allj);
   my_assert(op1.size2() == alli);
   my_assert(op2.size1() == allj);
   my_assert(op2.size2() == alll);
   my_assert(op3.size1() == alll);
   my_assert(op3.size2() == alli);
   my_assert(diagi.absenergy.size() == alli);
   my_assert(diagj.absenergy.size() == allj);
   my_assert(diagl.absenergy.size() == alll);
   my_assert(diagi.boltzmann.size() == alli);
   my_assert(diagj.boltzmann.size() == allj);
   my_assert(diagl.boltzmann.size() == alll);
#define SUMA13 SUM(1.0/(cmpl(Ei-Ej,wf(m))*cmpl(El-Ej,wf(m)+wb(n))));
   if (true) { // A1+A3, 1
//      cout << "A13, 1" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_X(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j,l) * op3(l, i) * wnf * bzEi;
	 SUMA13
      }}}
   }
   if (true) { // A1+A3, 2
//      cout << "A13, 2" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_X(i) { LOOP_D(j) { LOOP_X(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j,l) * op3(l, i) * wnf * bzEj;
         SUMA13
      }}}
   }
   if (retj>0 && alll>0) { // A1+A3, 3
//      cout << "A13, 3" << endl;
      // rho [retj,retj] x B [K=retj,D=alll]
      const matrix_range<const Matrix> op2cut(op2, range(0, retj), range(0, alll));
      Matrix rho_op2(retj, alll);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_X(i) { LOOP_K(j) { LOOP_D(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * rho_op2(j, l) * op3(l, i);
         SUMA13
      }}}
   }
   if (retj>0 && retl>0) { // A1+A3, 4
//      cout << "A13, 4" << endl;
      // rho [retj,retj] x B [K=retj,K=retl]
      const matrix_range<const Matrix> op2cut(op2, range(0, retj), range(0, retl));
      Matrix rho_op2(retj, retl);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhoj, op2cut, 0.0, rho_op2);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_K(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * rho_op2(j, l) * op3(l, i);
         SUMA13
      }}}
   }
   if (alll>0 && reti>0) { // A1+A3, 5
//      cout << "A13, 5" << endl;
      // C [D=alll, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, alll), range(0, reti));
      Matrix op3_rho(alll, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j, l) * op3_rho(l, i);
         SUMA13
      }}}
   }
   if (retl>0 && reti>0) { // A1+A3, 6
//      cout << "A13, 6" << endl;
      // C [K=retl, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix op3_rho(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j, l) * op3_rho(l, i);
         SUMA13
      }}}
   }
#undef SUMA13
#define SUMA24 SUM(1.0/(cmpl(El-Ei,wb(n))*cmpl(El-Ej,wf(m)+wb(n))));
   if (true) { // A2+A4, 1
//      cout << "A24, 1" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS // not LOOP_X(i), D case excluded
	 t_matel mat = CONJ_ME(op1(j,i)) * op2(j,l) * op3(l, i) * wnf * bzEl;
         SUMA24
      }}}
   }
   if (true) { // A2+A4, 2
//      cout << "A24, 2" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_K(l) { LOOPVARS // not LOOP_X(l), D case excluded
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j,l) * op3(l, i) * wnf * bzEi;
         SUMA24
      }}}
   }
   if (true) { // A2+A4, DD case
//      cout << "A24, DD" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS // DD case
	 t_matel mat = CONJ_ME(op1(j,i)) * op2(j,l) * op3(l, i) * wnf;
	 if (abs(mat*bzEl) < v3mmcutoff && abs(mat*bzEi) < v3mmcutoff) continue; // both!
	 for (short n = 0; n < maxn; n++) { // bosonic
	    for (short m = 0; m < maxm; m++) {
	       t_weight factor = boltz_fnc_noscale(El,Ei,bzEl,bzEi,wb(n),n)/cmpl(El-Ej,wb(n)+wf(m));
	       const t_weight w = spinfactor * mat * factor;
	       assert_isfinite(w);
	       res[omp_get_thread_num()](n,m) += w;
	    }
	 }
      }}}
   }
   if (retl>0 && alli>0) { // A2+A4, 3
//      cout << "A24, 3" << endl;
      // rho [retl,retl] x C [K=retl,D=alli]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, alli));
      Matrix rho_op3(retl, alli);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, rho_op3);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = CONJ_ME(op1(j,i)) * op2(j, l) * rho_op3(l, i);
         SUMA24
      }}}
   }
   if (retl>0 && reti>0) { // A2+A4, 4
//      cout << "A24, 4" << endl;
      // rho [retl,retl] x C [K=retl,K=reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix rho_op3(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, rho_op3);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = CONJ_ME(op1(j,i)) * op2(j, l) * rho_op3(l, i);
         SUMA24
      }}}
   }
   if (alll>0 && reti>0) { // A2+A4, 5
//      cout << "A24, 5" << endl;
      // C [D=alll, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, alll), range(0, reti));
      Matrix op3_rho(alll, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j, l) * op3_rho(l, i);
         SUMA24
      }}}
   }
   if (retl>0 && reti>0) { // A2+A4, 6
//      cout << "A24, 6" << endl;
      // C [K=retl, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix op3_rho(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(j,i)) * op2(j, l) * op3_rho(l, i);
         SUMA24
      }}}
   }
#undef SUMA24
   for (size_t j = 0; j < nr; j++)
      for (short n = 0; n < maxn; n++)
	 for (short m = 0; m < maxm; m++)
	    csm->add(m,n, res[j](n,m)); // note the order
}

void SPEC_FDM_v3mm::calc_B(const Eigen & diagi, const Eigen & diagj, const Eigen & diagl,
			   const Matrix & op1, const Matrix & op2, const Matrix & op3,
			   const BaseSpectrum &bs, t_factor spinfactor,
			   ChainSpectrum *cs, const Invar & Ii, const Invar & Ij, const Invar & Il)
{
//   cout << "B" << endl;
   const double v3mmcutoff = P::v3mmcutoff * sqr(P::T); // order of contributions is prop to 1/T^2
   size_t nr = omp_get_max_threads();
   nrglog('G', "nr=" << nr);
   const short maxn = P::mats;
   const short maxm = P::mats;
   res_t res(nr);
   for (size_t j = 0; j < nr; j++) {
      res[j].resize(maxn, maxm); // (n,m) order
      res[j].clear();
   }
   auto *csm = dynamic_cast<ChainSpectrumMatsubara2*>(cs);
   const double wnf = STAT::wnfactor[STAT::N];
   const Matrix & rhoi = rhoFDM[Ii];
   const Matrix & rhoj = rhoFDM[Ij];
   const Matrix & rhol = rhoFDM[Il];
   const size_t reti = (LAST_ITERATION() ? 0 : rhoi.size1()); // retained states
   const size_t retj = (LAST_ITERATION() ? 0 : rhoj.size1());
   const size_t retl = (LAST_ITERATION() ? 0 : rhol.size1());
   const size_t alli = diagi.getnr(); // all states
   const size_t allj = diagj.getnr();
   const size_t alll = diagl.getnr();
   my_assert(op1.size1() == alll);
   my_assert(op1.size2() == allj);
   my_assert(op2.size1() == alli);
   my_assert(op2.size2() == allj);
   my_assert(op3.size1() == alll);
   my_assert(op3.size2() == alli);
   my_assert(diagi.absenergy.size() == alli);
   my_assert(diagj.absenergy.size() == allj);
   my_assert(diagl.absenergy.size() == alll);
   my_assert(diagi.boltzmann.size() == alli);
   my_assert(diagj.boltzmann.size() == allj);
   my_assert(diagl.boltzmann.size() == alll);
#define SUMB12 SUM(1.0/(cmpl(Ej-Ei,wb(n)+wf(m))*cmpl(Ej-El,wf(m))));
   if (true) { // B1+B2, 1
//      cout << "B12, 1" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_X(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS
	 t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3(l,i) * wnf * bzEl;
	 SUMB12
      }}}
   }
   if (true) { // B1+B2, 2
//      cout << "B12, 2" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_X(i) { LOOP_D(j) { LOOP_X(l) { LOOPVARS
	 t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3(l,i) * wnf * bzEj;
	 SUMB12
      }}}
   }
   if (retj>0 && alll>0) { // B1+B2, 3
//      cout << "B12, 3" << endl;
      // rho [retj,retj] x A [K=retj,D=alll] (ConjTrans)
      const matrix_range<const Matrix> op1cut(op1, range(0, alll), range(0, retj)); // order
      Matrix rho_op1(retj, alll);
      atlas::gemm(CblasNoTrans, CblasConjTrans, 1.0, rhoj, op1cut, 0.0, rho_op1); // ConjTrans
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_X(i) { LOOP_K(j) { LOOP_D(l) { LOOPVARS
	    t_matel mat = rho_op1(j,l) * op2(i,j) * op3(l,i);
	    SUMB12
      }}}
   }
   if (retj>0 && retl>0) { // B1+B2, 4
//      cout << "B12, 4" << endl;
      // rho [retj,retj] x A [K=retj,K=retl] (ConjTrans)
      const matrix_range<const Matrix> op1cut(op1, range(0, retl), range(0, retj)); // order
      Matrix rho_op1(retj, retl);
      atlas::gemm(CblasNoTrans, CblasConjTrans, 1.0, rhoj, op1cut, 0.0, rho_op1); // ConjTrans
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_K(j) { LOOP_K(l) { LOOPVARS
	    t_matel mat = rho_op1(j,l) * op2(i,j) * op3(l,i);
	    SUMB12
      }}}
   }
   if (retl>0 && alli>0) { // B1+B2, 5
//      cout << "B12, 5" << endl;
      // rho [retl,retl] x C [K=retl,D=alli]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, alli));
      Matrix op3_rho(retl, alli);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_K(l) { LOOPVARS
	    t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3_rho(l,i);
	    SUMB12
      }}}
   }
   if (retl>0 && reti>0) { // B1+B2, 6
//      cout << "B12, 6" << endl;
      // rho [retl,retl] x C [K=retl,K=reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix op3_rho(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	    t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3_rho(l,i);
	    SUMB12
      }}}
   }
#undef SUMB12
#define SUMB34 SUM(1.0/(cmpl(Ej-Ei,wb(n)+wf(m))*cmpl(El-Ei,wb(n))));
   if (true) { // B3+B4, 1
//      cout << "B34, 1" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS // not LOOP_X(i), D case excluded
	 t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3(l,i) * wnf * bzEl;
	 SUMB34
      }}}
   }
   if (true) { // B3+B4, 2
//      cout << "B34, 2" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_K(l) { LOOPVARS // not LOOP_X(l), D case excluded
	 t_matel mat = - CONJ_ME(op1(l,j)) * op2(i,j) * op3(l,i) * wnf * bzEi;
	 SUMB34
      }}}
   }
   if (true) { // B3+B4, DD case
//      cout << "B34, 0/0" << endl;
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS  // DD case
	 t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * op3(l,i) * wnf;
	 if (abs(mat*bzEl) < v3mmcutoff && abs(mat*bzEi) < v3mmcutoff) continue;
	 for (short n = 0; n < maxn; n++) { // bosonic
	    for (short m = 0; m < maxm; m++) {
	       t_weight factor = boltz_fnc_noscale(El,Ei,bzEl,bzEi,wb(n),n)/cmpl(Ej-Ei,wb(n)+wf(m));
	       const t_weight w = spinfactor * mat * factor;
	       assert_isfinite(w);	
	       res[omp_get_thread_num()](n,m) += w;
	    }
	 }
      }}}
   }
   if (retl>0 && alli>0) { // B3+B4, 3
//      cout << "B34, 3" << endl;
      // rho [retl,retl] x C [K=retl,D=alli]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, alli));
      Matrix rho_op3(retl, alli);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, rho_op3);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_D(i) { LOOP_X(j) { LOOP_K(l) { LOOPVARS
	    t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * rho_op3(l,i);
	    SUMB34
      }}}
   }
   if (retl>0 && reti>0) { // B3+B4, 4
//      cout << "B34, 4" << endl;
      // rho [retl,retl] x C [K=retl,K=reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix rho_op3(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, rhol, op3cut, 0.0, rho_op3);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	    t_matel mat = CONJ_ME(op1(l,j)) * op2(i,j) * rho_op3(l,i);
	    SUMB34
      }}}
   }
   if (alll>0 && reti>0) { // B3+B4, 5
//      cout << "B34, 5" << endl;
      // C [D=alll, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, alll), range(0, reti));
      Matrix op3_rho(alll, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_X(j) { LOOP_D(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(l,j)) * op2(i,j) * op3_rho(l,i);
         SUMB34
      }}}
   }
   if (retl>0 && reti>0) { // B3+B4, 6
//      cout << "B34, 6" << endl;
      // C [K=retl, K=reti] x rho [reti,reti]
      const matrix_range<const Matrix> op3cut(op3, range(0, retl), range(0, reti));
      Matrix op3_rho(retl, reti);
      atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, op3cut, rhoi, 0.0, op3_rho);
#pragma omp parallel for schedule(dynamic) collapse(3)
      LOOP_K(i) { LOOP_D(j) { LOOP_K(l) { LOOPVARS
	 t_matel mat = -CONJ_ME(op1(l,j)) * op2(i,j) * op3_rho(l,i);
         SUMB34
      }}}
   }
#undef SUMB34      
   for (size_t j = 0; j < nr; j++)
      for (short n = 0; n < maxn; n++)
	 for (short m = 0; m < maxm; m++)
	    csm->add(m,n, res[j](n,m)); // note the order
}
