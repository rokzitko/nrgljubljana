// dmnrg.h - Density-matrix NRG
// Copyright (C) 2009-2021 Rok Zitko

#ifndef _dmnrg_hpp_
#define _dmnrg_hpp_

#include <memory>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cmath>
#include <chrono>
#include <limits>
#include <utility>

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

#include <fmt/format.h>

namespace NRG {

// Check if the trace of the density matrix equals 'ref_value'.
template<scalar S, typename MF>
void check_trace_rho(const DensMatElements<S> &m, MF mult, const double ref_value = 1.0) {
  const auto tr = m.trace(mult);
  if (!num_equal(tr, ref_value))
    throw std::runtime_error(fmt::format("check_trace_rho() failed, tr={}, ref_value={}", tr, ref_value));
}

class DmnrgProfile {
 public:
  using clock = std::chrono::steady_clock;
  using time_point = clock::time_point;
  enum class Contribution { generic, kk, dd };

  DmnrgProfile(std::string label, const size_t N, const bool enabled) : label_(std::move(label)), N_(N), enabled_(enabled) {}

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }
  [[nodiscard]] static time_point now() { return clock::now(); }
  [[nodiscard]] static double elapsed(const time_point start) { return std::chrono::duration<double>(clock::now() - start).count(); }

  void add_diag_load(const double seconds) { diag_load_seconds_ += seconds; }
  void add_init(const double seconds) { init_seconds_ += seconds; }
  void add_trace(const double seconds) { trace_seconds_ += seconds; }
  void add_save(const double seconds) { save_seconds_ += seconds; }
  void add_missing_rho() { ++missing_rho_; }
  void add_missing_diag() { ++missing_diag_; }
  void add_zero_dim() { ++zero_dim_; }
  void add_zero_omega() { ++zero_omega_; }
  void add_zero_rmax() { ++zero_rmax_; }

  void add_rotate(const size_t q, const size_t r, const size_t element_size, const double seconds,
                  const Contribution contribution = Contribution::generic, const bool diagonal = false) {
    ++contributions_;
    rotate_seconds_ += seconds;
    if (contribution == Contribution::kk) {
      ++kk_contributions_;
      kk_rotate_seconds_ += seconds;
    } else if (contribution == Contribution::dd) {
      ++dd_contributions_;
      dd_rotate_seconds_ += seconds;
    }
    q_min_ = std::min(q_min_, q);
    q_max_ = std::max(q_max_, q);
    r_min_ = std::min(r_min_, r);
    r_max_ = std::max(r_max_, r);
    if (diagonal)
      multiply_adds_ += static_cast<long double>(q) * static_cast<long double>(r)
                     + static_cast<long double>(q) * static_cast<long double>(r) * static_cast<long double>(r);
    else
      multiply_adds_ += static_cast<long double>(q) * static_cast<long double>(q) * static_cast<long double>(r)
                     + static_cast<long double>(q) * static_cast<long double>(r) * static_cast<long double>(r);
    u_bytes_ += static_cast<long double>(q) * static_cast<long double>(r) * static_cast<long double>(element_size);
    rho_bytes_ += static_cast<long double>(diagonal ? q : q * q) * static_cast<long double>(element_size);
    out_bytes_ += static_cast<long double>(r) * static_cast<long double>(r) * static_cast<long double>(element_size);
  }

  void report() const {
    if (!enabled_) return;
    const auto q_min = contributions_ ? q_min_ : 0;
    const auto r_min = contributions_ ? r_min_ : 0;
    std::cout << fmt::format(
        "[{} profile] N={} contrib={} rotate={:.6g}s kk={} kk_rotate={:.6g}s dd={} dd_rotate={:.6g}s diag_load={:.6g}s init={:.6g}s trace={:.6g}s save={:.6g}s "
        "q={}..{} r={}..{} madd={:.6g} U={:.6g}MB rho={:.6g}MB out={:.6g}MB missing_rho={} missing_diag={} zero_dim={} zero_omega={} zero_rmax={}",
        label_, N_, contributions_, rotate_seconds_, kk_contributions_, kk_rotate_seconds_, dd_contributions_, dd_rotate_seconds_,
        diag_load_seconds_, init_seconds_, trace_seconds_, save_seconds_,
        q_min, q_max_, r_min, r_max_, static_cast<double>(multiply_adds_), megabytes(u_bytes_), megabytes(rho_bytes_), megabytes(out_bytes_),
        missing_rho_, missing_diag_, zero_dim_, zero_omega_, zero_rmax_) << std::endl;
  }

 private:
  [[nodiscard]] static double megabytes(const long double bytes) { return static_cast<double>(bytes / (1024.0L * 1024.0L)); }

  std::string label_;
  size_t N_{};
  bool enabled_{};
  size_t contributions_{};
  size_t kk_contributions_{};
  size_t dd_contributions_{};
  size_t missing_rho_{};
  size_t missing_diag_{};
  size_t zero_dim_{};
  size_t zero_omega_{};
  size_t zero_rmax_{};
  size_t q_min_ = std::numeric_limits<size_t>::max();
  size_t q_max_{};
  size_t r_min_ = std::numeric_limits<size_t>::max();
  size_t r_max_{};
  long double multiply_adds_{};
  long double u_bytes_{};
  long double rho_bytes_{};
  long double out_bytes_{};
  double rotate_seconds_{};
  double kk_rotate_seconds_{};
  double dd_rotate_seconds_{};
  double diag_load_seconds_{};
  double init_seconds_{};
  double trace_seconds_{};
  double save_seconds_{};
};

// Calculate rho_N, the density matrix at the last NRG iteration. It is
// normalized to 1. Note: in CFS approach, we consider all states in the
// last iteration to be "discarded".
// For the details on the full Fock space approach see:
// F. B. Anders, A. Schiller, Phys. Rev. Lett. 95, 196801 (2005).
// F. B. Anders, A. Schiller, Phys. Rev. B 74, 245113 (2006).
// R. Peters, Th. Pruschke, F. B. Anders, Phys. Rev. B 74, 245114 (2006).
template<scalar S>
auto init_rho(const Step &step, const DiagInfo<S> &diag_in, const Symmetry<S> *Sym, const Params &P) {
  const auto mult = Sym->multfnc();
  if (P.project == ""s) {
    DensMatElements<S> rho;
    const auto Z = grand_canonical_Z(step.scT(), diag_in, mult);
    for (const auto &[I, eig] : diag_in)
      rho[I] = eig.diagonal_exp(step.scT()) / Z;
    if (P.checkrho) check_trace_rho(rho, mult);
    return rho;
  }
  const auto diag = Sym->project(diag_in, P.project);
  DensMatElements<S> rho;
  const auto Z = grand_canonical_Z(step.scT(), diag, mult);
  for (const auto &[I, eig] : diag)
    rho[I] = eig.diagonal_exp(step.scT()) / Z;
  if (P.checkrho) check_trace_rho(rho, mult);
  return rho;
}

// Calculation of the contribution from subspace I1 of rhoN (density matrix at iteration N) to rhoNEW (density matrix
// at iteration N-1)
template<scalar S, typename Matrix = Matrix_traits<S>, typename t_coef = coef_traits<S>>
void cdmI(const size_t i,        // Subspace index
          const Invar &I1,       // Quantum numbers corresponding to subspace i
          const Matrix &rhoN,    // rho^N
          const Eigen<S> &diagI1,   // contains U_{I1}
          Matrix &rhoNEW,        // rho^{N-1}
          const size_t N,
          const t_coef factor, // multiplicative factor that accounts for multiplicity
          const BackiterStore &store_all,
          const Params &P,
          DmnrgProfile *profile = nullptr,
          const DmnrgProfile::Contribution contribution = DmnrgProfile::Contribution::generic,
          const bool diagonal_rho = false)
{
  my_assert(i < P.combs);
  nrglog('D', "cdmI i=" << i << " I1=" << I1 << " factor=" << factor);
  // Range of indexes r and r' in matrix C^{QS,N}_{r,r'}, cf. Eq. (3.55) in my dissertation.
  const auto dim = size2(rhoNEW);
  // number of states taken into account in the density-matrix at *current* (Nth) stage (in subspace I1)
  const auto nromega = size2(rhoN);
  if (nromega == 0) { if (profile != nullptr) profile->add_zero_omega(); return; }
  if (dim == 0) { if (profile != nullptr) profile->add_zero_dim(); return; } // continue only if connection exists
  // rmax (info[I1].rmax[i]) is the range of r in U^N_I1(omega|ri), only those states that we actually kept..
  const auto rmax = store_all[N].at(I1).rmax.rmax(i);
  if (rmax == 0) { if (profile != nullptr) profile->add_zero_rmax(); return; } // rmax can be zero in the case a subspace has been completely truncated
  my_assert(rmax == dim);   // Otherwise, rmax must equal dim
  // Check range of omega: do the dimensions of C^N_I1(omega omega') and U^N_I1(omega|r1) match?
  my_assert(nromega <= diagI1.getnrstored());
  const auto &U0 = diagI1.U.get(i);
  const auto U = NRG::submatrix_const(U0, {0, nromega}, {0, size2(U0)});
  const auto rotate_start = profile != nullptr && profile->enabled() ? DmnrgProfile::now() : DmnrgProfile::time_point{};
  if (diagonal_rho)
    rotate_diagonal<S>(rhoNEW, factor, U, rhoN);
  else
    rotate<S>(rhoNEW, factor, U, rhoN);
  if (profile != nullptr && profile->enabled())
    profile->add_rotate(nromega, dim, sizeof(S), DmnrgProfile::elapsed(rotate_start), contribution, diagonal_rho);
}

// Calculation of the shell-N REDUCED DENSITY MATRICES: Calculate rho at previous iteration (N-1) from rho
// at the current iteration (N, rho)
template<scalar S>
auto calc_densitymatrix_iterN(const DiagInfo<S> &diag, const DensMatElements<S> &rho,
                              const size_t N, const BackiterStore &store_all, const Symmetry<S> *Sym, const Params &P,
                              DmnrgProfile *profile = nullptr) {
  nrglog('D', "calc_densitymatrix_iterN N=" << N);
  DensMatElements<S> rhoPrev;
  for (const auto &[I, dimsub] : store_all[N - 1]) { // loop over all subspaces at *previous* iteration
    const auto dim  = dimsub.kept();
    rhoPrev[I]      = zero_matrix<S>(dim);
    if (dim == 0) { if (profile != nullptr) profile->add_zero_dim(); continue; }
    const auto ns = Sym->new_subspaces(I);
    for (const auto &[i, sub] : ns | ranges::views::enumerate) {
      const auto x = rho.find(sub);
      const auto y = diag.find(sub);
      if (x == rho.end()) { if (profile != nullptr) profile->add_missing_rho(); continue; }
      if (y == diag.end()) { if (profile != nullptr) profile->add_missing_diag(); continue; }
      cdmI(i, sub, x->second, y->second, rhoPrev[I], N, double(Sym->mult(sub)) / double(Sym->mult(I)), store_all, P, profile);
    }
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
void calc_densitymatrix(DensMatElements<S> &rho, const BackiterStore &store_all, const Symmetry<S> *Sym,
                        MemTime &mt, const Params &P, const std::string filename = fn_rho) {
  if (P.resume && already_computed(filename, P)) return;
  if (P.checkrho) check_trace_rho(rho, Sym->multfnc()); // Must be 1.
  const auto section_timing = mt.time_it("DM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[DM] " << N << std::endl;
    DmnrgProfile profile("DM", N, P.logletter('Y'));
    auto timer = DmnrgProfile::now();
    const DiagInfo<S> diag_loaded(N, P);
    profile.add_diag_load(DmnrgProfile::elapsed(timer));
    auto rhoPrev = calc_densitymatrix_iterN(diag_loaded, rho, N, store_all, Sym, P, &profile); // need store_all for backiteration!
    if (P.checkrho) {
      timer = DmnrgProfile::now();
      check_trace_rho(rhoPrev, Sym->multfnc()); // Make sure rho is normalized to 1.
      profile.add_trace(DmnrgProfile::elapsed(timer));
    }
    timer = DmnrgProfile::now();
    rhoPrev.save(N-1, P, filename);
    profile.add_save(DmnrgProfile::elapsed(timer));
    profile.report();
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
DensMatElements<S> init_rho_FDM(const size_t N, const ThermoStore<S> &store, const Stats<S> &stats,
                                MF mult, const double T, const bool checkrho) {
  DensMatElements<S> rhoFDM;
  for (const auto &[I, ds] : store[N]) {
    rhoFDM[I] = zero_matrix<S>(ds.max());
    if (stats.ZnDNd[N] != 0.0)
      for (const auto i: ds.all())
        rhoFDM[I](i, i) = exp(-ds.eig.values.abs_zero(i) / T) * stats.wn[N] / stats.ZnDNd[N];
  }
  if (checkrho && stats.wn[N] != 0.0) { // note: wn \propto ZnDNd, so this is the same condition as above
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
                                  const size_t N, const ThermoStore<S> &store, const BackiterStore &store_all, const Stats<S> &stats,
                                  const Symmetry<S> *Sym, const Params &P,
                                  DmnrgProfile *profile = nullptr) {
  nrglog('D', "calc_fulldensitymatrix_iterN N=" << N);
  DensMatElements<S> rhoDD;
  DensMatElements<S> rhoFDMPrev;
  if (!step.last(N)) {
    const auto timer = profile != nullptr && profile->enabled() ? DmnrgProfile::now() : DmnrgProfile::time_point{};
    rhoDD = init_rho_FDM(N, store, stats, Sym->multfnc(), P.T, P.checkrho); // store here!
    if (profile != nullptr && profile->enabled()) profile->add_init(DmnrgProfile::elapsed(timer));
  }
  for (const auto &[I, ds] : store_all[N - 1]) { // loop over all subspaces at *previous* iteration, hence store_all here
    const auto subs = Sym->new_subspaces(I);
    const auto dim  = ds.kept();
    rhoFDMPrev[I]   = zero_matrix<S>(dim);
    if (!dim) { if (profile != nullptr) profile->add_zero_dim(); continue; }
    for (const auto i : Sym->combs()) {
      const auto sub = subs[i];
      // DM construction for non-Abelian symmetries: must include the ratio of multiplicities as a coefficient.
      const auto coef = double(Sym->mult(sub)) / double(Sym->mult(I));
      // Contribution from the KK sector.
      const auto x1 = rhoFDM.find(sub);
      const auto y = diag.find(sub);
      if (x1 == rhoFDM.end()) {
        if (profile != nullptr) profile->add_missing_rho();
      } else if (y == diag.end()) {
        if (profile != nullptr) profile->add_missing_diag();
      } else {
        cdmI(i, sub, x1->second, y->second, rhoFDMPrev[I], N, coef, store_all, P, profile, DmnrgProfile::Contribution::kk);
      }
      // Contribution from the DD sector. rhoDD -> rhoFDMPrev
      if (!step.last(N)) {
        if (const auto x2 = rhoDD.find(sub); x2 == rhoDD.end()) {
          if (profile != nullptr) profile->add_missing_rho();
        } else if (y == diag.end()) {
          if (profile != nullptr) profile->add_missing_diag();
        } else {
          cdmI(i, sub, x2->second, y->second, rhoFDMPrev[I], N, coef, store_all, P, profile, DmnrgProfile::Contribution::dd, true);
        }
      }
      // (Exception: for the N-1 iteration, the rhoPrev is already initialized with the DD sector of the last iteration.) }
    } // over combinations
  } // over subspaces
  return rhoFDMPrev;
}

template<scalar S>
void calc_fulldensitymatrix(const Step &step, DensMatElements<S> &rhoFDM, const ThermoStore<S> &store, const BackiterStore &store_all, const Stats<S> &stats,
                            const Symmetry<S> *Sym, MemTime &mt, const Params &P, const std::string &filename = fn_rhoFDM) {
  if (P.resume && already_computed(filename, P)) return;
  const auto section_timing = mt.time_it("FDM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    std::cout << "[FDM] " << N << std::endl;
    DmnrgProfile profile("FDM", N, P.logletter('Y'));
    auto timer = DmnrgProfile::now();
    const DiagInfo<S> diag_loaded(N, P); // = load_and_project(N, Sym, P);
    profile.add_diag_load(DmnrgProfile::elapsed(timer));
    auto rhoFDMPrev        = calc_fulldensitymatrix_iterN(step, diag_loaded, rhoFDM, N, store, store_all, stats, Sym, P, &profile);
    if (P.checkrho) {
      timer = DmnrgProfile::now();
      const auto tr          = rhoFDMPrev.trace(Sym->multfnc());
      const auto expected    = std::accumulate(stats.wn.begin() + N, stats.wn.begin() + P.Nmax, 0.0);
      const auto diff        = (tr - expected) / expected;
      nrglog('w', "tr[rhoFDM(" << N << ")]=" << tr << " sum(wn)=" << expected << " diff=" << diff);
      my_assert(num_equal(diff, 0.0));
      profile.add_trace(DmnrgProfile::elapsed(timer));
    }
    timer = DmnrgProfile::now();
    rhoFDMPrev.save(N-1, P, filename);
    profile.add_save(DmnrgProfile::elapsed(timer));
    profile.report();
    rhoFDM.swap(rhoFDMPrev);
  }
}

} // namespace

#endif
