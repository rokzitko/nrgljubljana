// dmnrg.h - Density-matrix NRG
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _dmnrg_h_
#define _dmnrg_h_

#include <fstream>
#include <iomanip>

#include "operators.h"
#include "symmetry.h"
#include "step.h"
#include "eigen.h"
#include "params.h"
#include "invar.h"
#include "traits.h"
#include "store.h"
#include "time_mem.h"

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>
namespace atlas = boost::numeric::bindings::atlas;

// Check if the trace of the density matrix equals 'ref_value'.
template<typename S>
void check_trace_rho(const DensMatElements<S> &m, std::shared_ptr<Symmetry<S>> Sym, const double ref_value = 1.0) {
  if (!num_equal(m.trace(Sym->multfnc()), ref_value))
    throw std::runtime_error("check_trace_rho() failed");
}

// Calculate rho_N, the density matrix at the last NRG iteration. It is
// normalized to 1. Note: in CFS approach, we consider all states in the
// last iteration to be "discarded".
// For the details on the full Fock space approach see:
// F. B. Anders, A. Schiller, Phys. Rev. Lett. 95, 196801 (2005).
// F. B. Anders, A. Schiller, Phys. Rev. B 74, 245113 (2006).
// R. Peters, Th. Pruschke, F. B. Anders, Phys. Rev. B 74, 245114 (2006).
template<typename S>
auto init_rho(const Step &step, const DiagInfo<S> &diag, std::shared_ptr<Symmetry<S>> Sym) {
  DensMatElements<S> rho;
  for (const auto &[I, eig]: diag)
    rho[I] = eig.diagonal_exp(step.scT()) / grand_canonical_Z(step, diag, Sym);
  check_trace_rho(rho, Sym);
  return rho;
}

// Calculation of the contribution from subspace I1 of rhoN (density matrix at iteration N) to rhoNEW (density matrix
// at iteration N-1)
template<typename S>
void cdmI(const size_t i,        // Subspace index
          const Invar &I1,       // Quantum numbers corresponding to subspace i
          const typename traits<S>::Matrix &rhoN,    // rho^N
          const Eigen<S> &diagI1,   // contains U_{I1}
          typename traits<S>::Matrix &rhoNEW,        // rho^{N-1}
          const size_t N,
          const typename traits<S>::t_coef factor, // multiplicative factor that accounts for multiplicity
          const AllSteps<S> &dm,
          const Params &P)
{
  my_assert(i < P.combs);
  nrglog('D', "cdmI i=" << i << " I1=" << I1 << " factor=" << factor);
  // Range of indexes r and r' in matrix C^{QS,N}_{r,r'}, cf. Eq. (3.55) in my dissertation.
  const auto dim = rhoNEW.size2();
  // number of states taken into account in the density-matrix at *current* (Nth) stage (in subspace I1)
  const auto nromega = rhoN.size2();
  if (nromega == 0 || dim == 0) return;   // continue only if connection exists
  // rmax (info[I1].rmax[i]) is the range of r in U^N_I1(omega|ri), only those states that we actually kept..
  const auto rmax = dm[N].at(I1).rmax.rmax(i);
  if (rmax == 0) return;    // rmax can be zero in the case a subspace has been completely truncated
  my_assert(rmax == dim);   // Otherwise, rmax must equal dim
  // Check range of omega: do the dimensions of C^N_I1(omega omega') and U^N_I1(omega|r1) match?
  my_assert(nromega <= diagI1.getnrstored());
  const ublas::matrix_range<const typename traits<S>::Matrix> U(diagI1.matrix, ublas::range(0, nromega), dm[N].at(I1).rmax.ubview(i));
  typename traits<S>::Matrix T(dim, nromega);
  using t_coef = typename traits<S>::t_coef;
  atlas::gemm(CblasConjTrans, CblasNoTrans, t_coef(1.0), U, rhoN, t_coef(0.0), T);    // T <- U^dag rhoN
  atlas::gemm(CblasNoTrans, CblasNoTrans, factor, T, U, t_coef(1.0), rhoNEW); // rhoNEW <- rhoNEW + factor T U
}

// Calculation of the shell-N REDUCED DENSITY MATRICES: Calculate rho at previous iteration (N-1, rhoPrev) from rho
// at the current iteration (N, rho)
template<typename S>
void calc_densitymatrix_iterN(const DiagInfo<S> &diag,
                              const DensMatElements<S> &rho, // input
                              DensMatElements<S> &rhoPrev,   // output
                              const size_t N, const AllSteps<S> &dm, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  nrglog('D', "calc_densitymatrix_iterN N=" << N);
  for (const auto &[I, dimsub] : dm[N - 1]) { // loop over all subspaces at *previous* iteration
    const auto dim  = dimsub.kept;
    rhoPrev[I]      = Zero_matrix<S>(dim);
    if (dim == 0) continue;
    const auto ns = Sym->new_subspaces(I);
    for (const auto &[i, sub] : ns | ranges::views::enumerate) {
      const auto x = rho.find(sub);
      const auto y = diag.find(sub);
      if (x != rho.end() && y != diag.end())
        cdmI(i, sub, x->second, y->second, rhoPrev[I], N, double(Sym->mult(sub)) / double(Sym->mult(I)), dm, P);
    }
  }
}

inline bool file_exists(const std::string &fn)
{
  std::ofstream F(fn, std::ios::binary | std::ios::out);
  return bool(F);
}

// Returns true if all the required density matrices are already saved on the disk.
inline bool already_computed(const std::string &prefix, const Params &P) {
  for (auto N = P.Nmax - 1; N > P.Ninit; N--) {
    const std::string fn = P.workdir.rhofn(N-1, prefix); // note the minus 1
    if (!file_exists(fn)) {
      std::cout << fn << " not found. Computing." << std::endl;
      return false;
    }
  }
  return true;
}

// calc_densitymatrix() is called prior to starting the NRG procedure for the second time. Here we calculate the
// shell-N density matrices for all iteration steps.
template<typename S>
void calc_densitymatrix(DensMatElements<S> &rho, const AllSteps<S> &dm, std::shared_ptr<Symmetry<S>> Sym,
                        MemTime &mt, const Params &P, const std::string filename = fn_rho) {
  if (P.resume && already_computed(filename, P)) return;
  check_trace_rho(rho, Sym); // Must be 1.
  if (P.ZBW) return;
  mt.time_it("DM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[DM] " << N << std::endl;
    DiagInfo<S> diag_loaded(N, P);
    DensMatElements<S> rhoPrev;
    calc_densitymatrix_iterN(diag_loaded, rho, rhoPrev, N, dm, Sym, P);
    check_trace_rho(rhoPrev, Sym); // Make sure rho is normalized to 1.
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
template<typename S>
DensMatElements<S> init_rho_FDM(const size_t N, const AllSteps<S> &dm, const Stats<S> &stats, 
                                     std::shared_ptr<Symmetry<S>> Sym, const double T) {
  DensMatElements<S> rhoFDM;
  for (const auto &[I, ds] : dm[N]) {
    rhoFDM[I] = Zero_matrix<S>(ds.max());
    if (stats.ZnDNd[N] != 0.0)
      for (const auto i: ds.all())
        rhoFDM[I](i, i) = exp(-ds.eig.absenergyN[i] / T) * stats.wn[N] / stats.ZnDNd[N];
  }
  if (stats.wn[N] != 0.0) { // note: wn \propto ZnDNd, so this is the same condition as above
    // Trace should be equal to the total weight of the shell-N contribution to the FDM.
    const auto tr = rhoFDM.trace(Sym->multfnc());
    const auto diff = (tr - stats.wn[N]) / stats.wn[N]; // relative error
    if (!num_equal(diff, 0.0, 1e-8)) my_assert(stats.wn[N] < 1e-12); // OK if small enough overall
  }
  return rhoFDM;
}

template<typename S>
void calc_fulldensitymatrix_iterN(const Step &step, // only required for step::last()
                                  const DiagInfo<S> &diag,
                                  const DensMatElements<S> &rhoFDM, // input
                                  DensMatElements<S> &rhoFDMPrev,   // output
                                  const size_t N, const AllSteps<S> &dm, const Stats<S> &stats,
                                  std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  nrglog('D', "calc_fulldensitymatrix_iterN N=" << N);
  DensMatElements<S> rhoDD;
  if (!step.last(N))
    rhoDD = init_rho_FDM(N, dm, stats, Sym, P.T);
  for (const auto &[I, ds] : dm[N - 1]) { // loop over all subspaces at *previous* iteration
    const auto subs = Sym->new_subspaces(I);
    const auto dim  = ds.kept;
    rhoFDMPrev[I]   = Zero_matrix<S>(dim);
    if (!dim) continue;
    for (const auto i : Sym->combs()) {
      const auto sub = subs[i];
      // DM construction for non-Abelian symmetries: must include the ratio of multiplicities as a coefficient.
      const auto coef = double(Sym->mult(sub)) / double(Sym->mult(I));
      // Contribution from the KK sector.
      const auto x1 = rhoFDM.find(sub);
      const auto y = diag.find(sub);
      if (x1 != rhoFDM.end() && y != diag.end())
        cdmI(i, sub, x1->second, y->second, rhoFDMPrev[I], N, coef, dm, P);
      // Contribution from the DD sector. rhoDD -> rhoFDMPrev
      if (!step.last(N))
        if (const auto x2 = rhoDD.find(sub); x2 !=rhoDD.end() && y != diag.end())
          cdmI(i, sub, x2->second, y->second, rhoFDMPrev[I], N, coef, dm, P);
      // (Exception: for the N-1 iteration, the rhoPrev is already initialized with the DD sector of the last iteration.) }
    } // over combinations
  } // over subspaces
}

// Sum of statistical weights from site N to the end of the Wilson chain.
template<typename S>
auto sum_wn(const size_t N, const Stats<S> &stats, const Params &P) {
  auto sum = 0.0;
  for (size_t n = P.Nmax - 1; n >= N; n--) sum += stats.wn[n]; // XXX: accumulate
  return sum;
}

template<typename S>
void calc_fulldensitymatrix(const Step &step, DensMatElements<S> &rhoFDM, const AllSteps<S> &dm, const Stats<S> &stats,
                            std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P, const std::string filename = fn_rhoFDM) {
  if (P.resume && already_computed(filename, P)) return;
  if (P.ZBW) return;
  mt.time_it("FDM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[FDM] " << N << std::endl;
    DiagInfo<S> diag_loaded(N, P);
    DensMatElements<S> rhoFDMPrev;
    calc_fulldensitymatrix_iterN(step, diag_loaded, rhoFDM, rhoFDMPrev, N, dm, stats, Sym, P);
    const auto tr       = rhoFDMPrev.trace(Sym->multfnc());
    const auto expected = sum_wn(N, stats, P);
    const auto diff     = (tr - expected) / expected;
    nrglog('w', "tr[rhoFDM(" << N << ")]=" << tr << " sum(wn)=" << expected << " diff=" << diff);
    my_assert(num_equal(diff, 0.0));
    rhoFDMPrev.save(N-1, P, filename);
    rhoFDM.swap(rhoFDMPrev);
  }
}

#endif // _dmnrg_h_
