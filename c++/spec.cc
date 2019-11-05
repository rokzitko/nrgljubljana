#ifndef _spec_cc_
# define _spec_cc_

// GENERAL NOTE: CONJ_ME() and reversal of indexes is appropriate for
// spectral functions of fermionic operators, since we need to
// distinguish between d and d^\dag. For generic G_{AB}, this is an
// unnecessary complication and confusing. For real symmetric
// operators, the result is not affected. For complex Hermitian, the
// same is true. But it is ugly.
 
// For argument 'sign' in calc_generic_* functions. This is the only
// difference between the GFs for bosonic and fermionic operators.
const double S_FERMIONIC = -1.0;
const double S_BOSONIC = 1.0;

#include "spec_FT.cc"
#include "spec_DMNRG.cc"
#include "spec_FDM.cc"
#include "spec_CFS.cc"

// Calculate (finite temperature) spectral function 1/Pi Im << op1^\dag(t)
// op2(0) >>. Required spin direction is determined by 'SPIN'. For SPIN=0
// both spin direction are equivalent. For QSZ, we need to differentiate
// the two.

template <typename FactorFnc, typename CheckSpinFnc>
void calc_generic(const BaseSpectrum &bs,
		  const DiagInfo &diag, 
		  FactorFnc &factorfnc,
		  CheckSpinFnc &checkspinfnc)
{
   nrglog('g', "calc_generic() " << bs.fullname());
   // Strategy: we loop through all subspace pairs and check whether
   // they have non-zero irreducible matrix elements.
   std::vector<Twoinvar> tasks;
   LOOP_const(diag, i)
      LOOP_const(diag, j) {
	 const Twoinvar II = make_pair(INVAR(j), INVAR(i));
	 if (bs.op1.count(II) && bs.op2.count(II))
	    if (checkspinfnc(INVAR(j), INVAR(i), bs.spin))
	       tasks.push_back(II);
      }
   ChainSpectrum *cs = bs.spectype->make_cs(bs);
   for (const auto & II : tasks) {
      Invar Ij, Ii;
      tie(Ij, Ii) = II;
      const Eigen & diagi = diag.find(Ii)->second;
      const Eigen & diagj = diag.find(Ij)->second;
      // Note that rho is always "diagonal in subspace indexes", i.e.
      // rho(I1,I2)=delta_(I1,I2) rho_I1. There is thus only one
      // multiplicity factor in play here.
      const t_factor spinfactor = factorfnc(Ii, Ij);
      my_assert(!num_equal(spinfactor, 0.0)); // bug trap
      const Matrix & op1II = bs.op1.find(II)->second;
      const Matrix & op2II = bs.op2.find(II)->second;
      if (logletter('G')) nrgdump2(Ij, Ii) << endl;
      bs.spectype->calc(diagi, diagj, op1II, op2II, bs, spinfactor, cs, Ii, Ij);
   }
   nrglog('*', "Merging " << bs.fullname());
   bs.spec->merge(cs);
   delete cs;
}

template <typename FactorFnc>
void calc_generic3(const BaseSpectrum &bs,
		   const DiagInfo &diag,
		   FactorFnc &factorfnc)
{
   nrglog('g', "calc_generic3() " << bs.fullname());
   ChainSpectrum *cs = bs.spectype->make_cs(bs);
   LOOP_const(diag, i) 
      LOOP_const(diag, j) 
	 LOOP_const(diag, l) {
	    const Invar Ii = INVAR(i);
	    const Invar Ij = INVAR(j);
	    const Invar Il = INVAR(l);
	    const Eigen &diagi = diag.find(Ii)->second;
	    const Eigen &diagj = diag.find(Ij)->second;
	    const Eigen &diagl = diag.find(Il)->second;
	    const auto cji = make_pair(Ij, Ii); // conj
	    const auto jl = make_pair(Ij, Il);
	    const auto li = make_pair(Il, Ii);
	    // A_ij x B_jl x C_li, A=op1, B=op2, C=op3
	    if (bs.op1.count(cji) && bs.op2.count(jl) && bs.op3.count(li)) {
	       if (logletter('G')) nrgdump3(Ii, Ij, Il) << endl;
	       const t_factor spinfactor = factorfnc(Ii, Ij);
	       my_assert(!num_equal(spinfactor, 0.0)); // bug trap
	       const Matrix &op1 = bs.op1.find(cji)->second; // conj : A_ij=(a)_ij=(a+)_ji*
	       const Matrix &op2 = bs.op2.find(jl)->second; // B_jl=(b+)_jl
	       const Matrix &op3 = bs.op3.find(li)->second; // C_li=n_li
	       bs.spectype->calc_A(diagi, diagj, diagl, op1, op2, op3, bs, spinfactor, cs, Ii, Ij, Il);
	    }
	    const auto ij = make_pair(Ii, Ij);
	    const auto clj = make_pair(Il, Ij); // conj
	    if (bs.op1.count(clj) && bs.op2.count(ij) && bs.op3.count(li)) {
	       if (logletter('G')) nrgdump3(Ii, Ij, Il) << endl;
	       const t_factor spinfactor = factorfnc(Ii, Ij); // XXX
	       my_assert(!num_equal(spinfactor, 0.0)); // bug trap
	       const Matrix &op1 = bs.op1.find(clj)->second; // conj : A_jl=(a)_jl=(a+)_lj*
	       const Matrix &op2 = bs.op2.find(ij)->second; // B_ij=(b+)_ij
	       const Matrix &op3 = bs.op3.find(li)->second; // C_li=n_li
	       bs.spectype->calc_B(diagi, diagj, diagl, op1, op2, op3, bs, spinfactor, cs, Ii, Ij, Il);
	    }
	 }
   nrglog('*', "Merging " << bs.fullname());
   bs.spec->merge(cs);
   delete cs;
}

#endif // _spec_cc_
