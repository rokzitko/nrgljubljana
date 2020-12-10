// dmnrg.h - Density-matrix NRG
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _dmnrg_hpp_
#define _dmnrg_hpp_

#include <memory>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cmath>

#include "operators.hpp"
#include "symmetry.hpp"
#include "step.hpp"
#include "eigen.hpp"
#include "params.hpp"
#include "invar.hpp"
#include "traits.hpp"
#include "store.hpp"
#include "time_mem.hpp"
#include "stats.hpp"
#include "numerics.hpp"

namespace NRG {

// Check if the trace of the density matrix equals 'ref_value'.
template<scalar S, typename MF>
void check_trace_rho(const DensMatElements<S> &m, MF mult, const double ref_value = 1.0) {
  if (!num_equal(m.trace(mult), ref_value))
    throw std::runtime_error("check_trace_rho() failed");
}

// Calculate rho_N, the density matrix at the last NRG iteration. It is
// normalized to 1. Note: in CFS approach, we consider all states in the
// last iteration to be "discarded".
// For the details on the full Fock space approach see:
// F. B. Anders, A. Schiller, Phys. Rev. Lett. 95, 196801 (2005).
// F. B. Anders, A. Schiller, Phys. Rev. B 74, 245113 (2006).
// R. Peters, Th. Pruschke, F. B. Anders, Phys. Rev. B 74, 245114 (2006).
template<scalar S, typename MF>
auto init_rho(const Step &step, const DiagInfo<S> &diag, MF mult) {
  DensMatElements<S> rho;
  for (const auto &[I, eig]: diag)
    rho[I] = eig.diagonal_exp(step.scT()) / grand_canonical_Z(step.scT(), diag, mult);
    // NOTE: diagonal_exp() and grand_canonical_Z() both use the round-off-error corrected eigenvalues.
  check_trace_rho(rho, mult);
  return rho;
}

// Calculation of the contribution from subspace I1 of rhoN (density matrix at iteration N) to rhoNEW (density matrix
// at iteration N-1)
template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void cdmI(const size_t i,          // Subspace index
          const Matrix &rhoN,      // rho^N
          const Eigen<S> &diagI1,  // contains U_{I1}
          Matrix &rhoNEW,          // rho^{N-1}
          const t_coef factor,     // multiplicative factor that accounts for multiplicity
          const Params &P)
{
  nrglog('D', "cdmI i=" << i << " factor=" << factor);
  const auto U = diagI1.U.get(i);
  rotate<S>(rhoNEW, factor, U, rhoN);
}

// Calculation of the shell-N REDUCED DENSITY MATRICES: Calculate rho at previous iteration (N-1) from rho
// at the current iteration (N, rho)
template<scalar S>
auto calc_densitymatrix_iterN(const DiagInfo<S> &diag, const DensMatElements<S> &rho,
                              const size_t N, const Store<S> &store, const Symmetry<S> *Sym, const Params &P) {
  nrglog('D', "calc_densitymatrix_iterN N=" << N);
  DensMatElements<S> rhoPrev;
  for (const auto &[I, dimsub] : store[N-1]) { // loop over all subspaces at *previous* iteration
    const auto dim  = dimsub.kept();
    rhoPrev[I]      = Zero_matrix<S>(dim);
    if (dim == 0) continue;
    const auto ns = Sym->new_subspaces(I);
    for (const auto &[i, sub] : ns | ranges::views::enumerate)
      if (const auto x = rho.find(sub); x != rho.end())
        if (const auto y = diag.find(sub); y != diag.end())
          cdmI(i, x->second, y->second, rhoPrev[I], double(Sym->mult(sub)) / double(Sym->mult(I)), P);
  }
  return rhoPrev;
}

// Returns true if all the required density matrices are already saved on the disk.
inline bool already_computed(const std::string &prefix, const Params &P) {
  for (auto N = P.Nmax - 1; N > P.Ninit; N--) {
    const std::string fn = P.workdir->rhofn(N-1, prefix); // note the minus 1
    if (!file_exists(fn)) {
      std::cout << fn << " not found. Computing." << std::endl;
      return false;
    }
  }
  return true;
}

// calc_densitymatrix() is called prior to starting the NRG procedure for the second time. Here we calculate the
// shell-N density matrices for all iteration steps.
template<scalar S>
void calc_densitymatrix(DensMatElements<S> &rho, const Store<S> &store, const Symmetry<S> *Sym,
                        MemTime &mt, const Params &P, const std::string filename = fn_rho) {
  if (P.resume && already_computed(filename, P)) return;
  check_trace_rho(rho, Sym->multfnc()); // Must be 1.
  if (P.ZBW) return;
  const auto section_timing = mt.time_it("DM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[DM] " << N << std::endl;
    DiagInfo<S> diag_loaded(N, P);
    auto rhoPrev = calc_densitymatrix_iterN(diag_loaded, rho, N, store, Sym, P);
    check_trace_rho(rhoPrev, Sym->multfnc()); // Make sure rho is normalized to 1.
    rhoPrev.save(N-1, P, filename);
    rho.swap(rhoPrev);
  }
}

// ****************** Calculation of the FULL REDUCED DENSITY MATIRICES

// Calculate rho(N), the shell-N density matrix, computed using
// the discarded states at shell N.
// Must be called AFTER calc_ZnD().
// Called fron nrg.cc immediately after the first NRG run (with N=Nmax-1),
// and also from calc_fulldensity_iterN (with lower N).
// A. Weichselbaum, J. von Delft, Phys. Rev. Lett. 99, 076402 (2007)
// T. A. Costi, V. Zlatic, Phys. Rev. B 81, 235127 (2010)
// H. Zhang, X. C. Xie, Q. Sun, Phys. Rev. B 82, 075111 (2010)
template<scalar S, typename MF>
DensMatElements<S> init_rho_FDM(const size_t N, const Store<S> &store, const Stats<S> &stats, 
                                MF mult, const double T) {
  DensMatElements<S> rhoFDM;
  for (const auto &[I, ds] : store[N]) {
    rhoFDM[I] = Zero_matrix<S>(ds.max());
    if (stats.ZnDNd[N] != 0.0)
      for (const auto i: ds.all())
        rhoFDM[I](i, i) = exp(-ds.eig.values.abs_zero(i) / T) * stats.wn[N] / stats.ZnDNd[N];
  }
  if (stats.wn[N] != 0.0) { // note: wn \propto ZnDNd, so this is the same condition as above
    // Trace should be equal to the total weight of the shell-N contribution to the FDM.
    const auto tr = rhoFDM.trace(mult);
    const auto diff = (tr - stats.wn[N]) / stats.wn[N]; // relative error
    if (!num_equal(diff, 0.0, 1e-8)) my_assert(stats.wn[N] < 1e-12); // OK if small enough overall
  }
  return rhoFDM;
}

template<scalar S>
auto calc_fulldensitymatrix_iterN(const Step &step, // only required for step::last()
                                  const DiagInfo<S> &diag,
                                  const DensMatElements<S> &rhoFDM, // input
                                  const size_t N, const Store<S> &store, const Stats<S> &stats,
                                  const Symmetry<S> *Sym, const Params &P) {
  nrglog('D', "calc_fulldensitymatrix_iterN N=" << N);
  DensMatElements<S> rhoDD;
  DensMatElements<S> rhoFDMPrev;
  if (!step.last(N))
    rhoDD = init_rho_FDM(N, store, stats, Sym->multfnc(), P.T);
  for (const auto &[I, ds] : store[N-1]) { // loop over all subspaces at *previous* iteration
    const auto subs = Sym->new_subspaces(I);
    const auto dim  = ds.kept();
    rhoFDMPrev[I]   = Zero_matrix<S>(dim);
    if (!dim) continue;
    for (const auto i : Sym->combs()) {
      const auto sub = subs[i];
      // DM construction for non-Abelian symmetries: must include the ratio of multiplicities as a coefficient.
      const auto coef = double(Sym->mult(sub)) / double(Sym->mult(I));
      if (const auto y = diag.find(sub); y != diag.end()) {
        // Contribution from the KK sector.
        if (const auto x1 = rhoFDM.find(sub); x1 != rhoFDM.end())
          cdmI(i, x1->second, y->second, rhoFDMPrev[I], coef, P);
        // Contribution from the DD sector. rhoDD -> rhoFDMPrev
        if (!step.last(N))
          if (const auto x2 = rhoDD.find(sub); x2 != rhoDD.end())
            cdmI(i, x2->second, y->second, rhoFDMPrev[I], coef, P);
      }
    } // over combinations
  } // over subspaces
  return rhoFDMPrev;
}

template<scalar S>
void calc_fulldensitymatrix(const Step &step, DensMatElements<S> &rhoFDM, const Store<S> &store, const Stats<S> &stats,
                            const Symmetry<S> *Sym, MemTime &mt, const Params &P, const std::string &filename = fn_rhoFDM) {
  if (P.resume && already_computed(filename, P)) return;
  if (P.ZBW) return;
  const auto section_timing = mt.time_it("FDM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[FDM] " << N << std::endl;
    DiagInfo<S> diag_loaded(N, P);
    auto rhoFDMPrev = calc_fulldensitymatrix_iterN(step, diag_loaded, rhoFDM, N, store, stats, Sym, P);
    const auto tr       = rhoFDMPrev.trace(Sym->multfnc());
    const auto expected = std::accumulate(stats.wn.begin() + N, stats.wn.begin() + P.Nmax, 0.0);
    const auto diff     = (tr - expected) / expected;
    nrglog('w', "tr[rhoFDM(" << N << ")]=" << tr << " sum(wn)=" << expected << " diff=" << diff);
    my_assert(num_equal(diff, 0.0));
    rhoFDMPrev.save(N-1, P, filename);
    rhoFDM.swap(rhoFDMPrev);
  }
}

} // namespace

#endif
