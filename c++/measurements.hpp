#ifndef _measurements_hpp_
#define _measurements_hpp_

#include "step.hpp"
#include "eigen.hpp"
#include "operators.hpp"
#include "symmetry.hpp"
#include "traits.hpp"
#include "output.hpp"
#include "store.hpp"
#include "time_mem.hpp"
#include "oprecalc.hpp"

namespace NRG {

template<typename S>
CONSTFNC auto calc_trace_singlet(const Step &step, const DiagInfo<S> &diag,
                                 const MatrixElements<S> &n, std::shared_ptr<Symmetry<S>> Sym) {
  typename traits<S>::t_matel tr{};
  for (const auto &[I, eig] : diag) {
    const auto & nI = n.at({I,I});
    const auto dim = eig.getnrstored();
    my_assert(dim == nI.size2());
    typename traits<S>::t_matel sum{};
    for (const auto r : range0(dim)) sum += exp(-step.TD_factor() * eig.value_zero(r)) * nI(r, r);
    tr += Sym->mult(I) * sum;
  }
  return tr; // note: t_expv = t_matel
}

// Measure thermodynamic expectation values of singlet operators
template<typename S>
void measure_singlet(const Step &step, Stats<S> &stats, const DiagInfo<S> &diag, const IterInfo<S> &a,
                            Output<S> &output, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  const auto Z = ranges::accumulate(diag, 0.0, [&Sym, &step](auto total, const auto &d) { const auto &[I, eig] = d;
    return total + Sym->mult(I) * ranges::accumulate(eig.value_zero, 0.0,
                                                     [f=step.TD_factor()](auto sum, const auto &x) { return sum + exp(-f*x); }); });
  for (const auto &[name, m] : a.ops)  stats.expv[name] = calc_trace_singlet(step, diag, m, Sym) / Z;
  for (const auto &[name, m] : a.opsg) stats.expv[name] = calc_trace_singlet(step, diag, m, Sym) / Z;
  output.custom->field_values(step.Teff());
}

template<typename T>
T trace_contract(const ublas::matrix<T> &A, const ublas::matrix<T> &B, const size_t range)
{
  T sum{};
  for (const auto i : range0(range))
       for (const auto j : range0(range))
      sum += A(i, j) * B(j, i);
  return sum;
}

template<typename S>
CONSTFNC auto calc_trace_fdm_kept(const size_t ndx, const MatrixElements<S> &n, const DensMatElements<S> &rhoFDM,
                                  const Store<S> &store, std::shared_ptr<Symmetry<S>> Sym) {
  typename traits<S>::t_matel tr{};
  for (const auto &[I, rhoI] : rhoFDM)
    tr += Sym->mult(I) * trace_contract(rhoI, n.at({I,I}), store[ndx].at(I).kept()); // over kept states ONLY
  return tr;
}

template<typename S>
void measure_singlet_fdm(const Step &step, Stats<S> &stats, const DiagInfo<S> &diag, const IterInfo<S> &a,
                         Output<S> &output,  const DensMatElements<S> &rhoFDM,
                         const Store<S> &store, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  for (const auto &[name, m] : a.ops)  stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, store, Sym);
  for (const auto &[name, m] : a.opsg) stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, store, Sym);
  output.customfdm->field_values(P.T);
}

// Calculate grand canonical partition function at current NRG energy shell. This is not the same as the true
// partition function of the full problem! Instead this is the Z_N that is used to initialize the density matrix,
// i.e. rho = 1/Z_N \sum_{l} exp{-beta E_l} |l;N> <l;N|.  grand_canonical_Z() is also used to calculate stats.Zft,
// that is used to compute the spectral function with the conventional approach, as well as stats.Zgt for G(T)
// calculations, stats.Zchit for chi(T) calculations.
template<typename S>
auto grand_canonical_Z(const Step &step, const DiagInfo<S> &diag, std::shared_ptr<Symmetry<S>> Sym, const double factor = 1.0) {
  double ZN{};
  for (const auto &[I, eig]: diag) 
    for (const auto &i : eig.kept()) // sum over all kept states
      ZN += Sym->mult(I) * exp(-eig.value_zero(i) * step.scT() * factor);
  my_assert(ZN >= 1.0);
  return ZN;
}

// Calculate partial statistical sums, ZnD*, and the grand canonical Z (stats.ZZG), computed with respect to absolute
// energies. calc_ZnD() must be called before the second NRG run.
template<typename S>
void calc_ZnD(const Store<S> &store, Stats<S> &stats, std::shared_ptr<Symmetry<S>> Sym, const double T) {
  mpf_set_default_prec(400); // this is the number of bits, not decimal digits!
  for (const auto N : store.Nall()) {
    my_mpf ZnDG, ZnDN; // arbitrary-precision accumulators to avoid precision loss
    mpf_set_d(ZnDG, 0.0);
    mpf_set_d(ZnDN, 0.0);
    for (const auto &[I, ds] : store[N])
      for (const auto i : ds.all()) {
        my_mpf g, n;
        mpf_set_d(g, Sym->mult(I) * exp(-ds.eig.absenergyG[i]/T)); // absenergyG >= 0.0
        mpf_set_d(n, Sym->mult(I) * exp(-ds.eig.absenergyN[i]/T)); // absenergyN >= 0.0
        mpf_add(ZnDG, ZnDG, g);
        mpf_add(ZnDN, ZnDN, n);
      }
    mpf_set(stats.ZnDG[N], ZnDG);
    mpf_set(stats.ZnDN[N], ZnDN);
    stats.ZnDNd[N] = mpf_get_d(stats.ZnDN[N]);
  }
  // Note: for ZBW, Nlen=Nmax+1. For Ninit=Nmax=0, index 0 will thus be included here.
  my_mpf ZZG;
  mpf_set_d(ZZG, 0.0);
  for (const auto N : store.Nall()) {
    my_mpf a;
    mpf_set(a, stats.ZnDG[N]);
    my_mpf b;
    mpf_set_d(b, Sym->nr_combs());
    mpf_pow_ui(b, b, store.Nend - N - 1);
    my_mpf c;
    mpf_mul(c, a, b);
    mpf_add(ZZG, ZZG, c);
  }
  stats.ZZG = mpf_get_d(ZZG);
  std::cout << "ZZG=" << HIGHPREC(stats.ZZG) << std::endl;
  for (const auto N : store.Nall()) {
    const double w  = std::pow(Sym->nr_combs(), store.Nend - N - 1) / stats.ZZG; // ZZZ
    stats.wnfactor[N] = w; // These ratios enter the terms for the spectral function.
    stats.wn[N] = w * mpf_get_d(stats.ZnDG[N]); // This is w_n defined after Eq. (8) in the WvD paper.
  }
  const auto sumwn = ranges::accumulate(stats.wn, 0.0);
  std::cout << "sumwn=" << sumwn << " sumwn-1=" << sumwn - 1.0 << std::endl;
  my_assert(num_equal(sumwn, 1.0));  // Check the sum-rule.
}

template<typename S>
void report_ZnD(Stats<S> &stats, const Params &P) {
  for (const auto N : P.Nall())
    std::cout << "ZG[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDG[N])) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "ZN[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDN[N])) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "w[" << N << "]=" << HIGHPREC(stats.wn[N]) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "wfactor[" << N << "]=" << HIGHPREC(stats.wnfactor[N]) << std::endl;
}

// TO DO: use Boost.Multiprecision instead of low-level GMP calls
// https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html
template<typename S>
void fdm_thermodynamics(const Store<S> &store, Stats<S> &stats, std::shared_ptr<Symmetry<S>> Sym, const double T)
{
  stats.td_fdm.T = T;
  stats.Z_fdm = stats.ZZG*exp(-stats.GS_energy/T); // this is the true partition function
  stats.td_fdm.F = stats.F_fdm = -log(stats.ZZG)*T+stats.GS_energy; // F = -k_B*T*log(Z)
  // We use multiple precision arithmetics to ensure sufficient accuracy in the calculation of
  // the variance of energy and thus the heat capacity.
  my_mpf E, E2;
  mpf_set_d(E, 0.0);
  mpf_set_d(E2, 0.0);
  for (const auto N : store.Nall())
    if (stats.wn[N] > 1e-16)
      for (const auto &[I, ds] : store[N])
        for (const auto i : ds.all()) {
          my_mpf weight;
          mpf_set_d(weight, stats.wn[N] * Sym->mult(I) * exp(-ds.eig.absenergyN[i]/T));
          mpf_div(weight, weight, stats.ZnDN[N]);
          my_mpf e;
          mpf_set_d(e, ds.eig.absenergy[i]);
          my_mpf e2;
          mpf_mul(e2, e, e);
          mpf_mul(e, e, weight);
          mpf_mul(e2, e2, weight);
          mpf_add(E, E, e);
          mpf_add(E2, E2, e2);
        }
  stats.td_fdm.E = stats.E_fdm = mpf_get_d(E);
  my_mpf sqrE;
  mpf_mul(sqrE, E, E);
  my_mpf varE;
  mpf_sub(varE, E2, sqrE);
  stats.td_fdm.C = stats.C_fdm = mpf_get_d(varE)/pow(T,2);
  stats.td_fdm.S = stats.S_fdm = (stats.E_fdm-stats.F_fdm)/T;
  std::cout << std::endl;
  std::cout << "Z_fdm=" << HIGHPREC(stats.Z_fdm) << std::endl;
  std::cout << "F_fdm=" << HIGHPREC(stats.F_fdm) << std::endl;
  std::cout << "E_fdm=" << HIGHPREC(stats.E_fdm) << std::endl;
  std::cout << "C_fdm=" << HIGHPREC(stats.C_fdm) << std::endl;
  std::cout << "S_fdm=" << HIGHPREC(stats.S_fdm) << std::endl;
  std::cout << std::endl;
  stats.td_fdm.save_values();
}

// We calculate thermodynamic quantities before truncation to make better use of the available states. Here we
// compute quantities which are defined for all symmetry types. Other calculations are performed by calculate_TD
// member functions defined in symmetry.h
template<typename S>
void calculate_TD(const Step &step, const DiagInfo<S> &diag, Stats<S> &stats, Output<S> &output, 
                  std::shared_ptr<Symmetry<S>> Sym, const double additional_factor = 1.0) {
  // Rescale factor for energies. The energies are expressed in units of omega_N, thus we need to appropriately
  // rescale them to calculate the Boltzmann weights at the temperature scale Teff (Teff=scale/betabar).
  const auto rescale_factor = step.TD_factor() * additional_factor;
  auto mult = [Sym](const auto &I) { return Sym->mult(I); };
  const auto Z  = diag.trace([](double x) { return 1; },        rescale_factor, mult); // partition function
  const auto E  = diag.trace([](double x) { return x; },        rescale_factor, mult); // Tr[beta H]
  const auto E2 = diag.trace([](double x) { return pow(x,2); }, rescale_factor, mult); // Tr[(beta H)^2]
  stats.Z = Z;
  stats.td.T  = step.Teff();
  stats.td.E  = E/Z;               // beta <H>
  stats.td.E2 = E2/Z;              // beta^2 <H^2>
  stats.td.C  = E2/Z - pow(E/Z,2); // C/k_B=beta^2(<H^2>-<H>^2)
  stats.td.F  = -log(Z);           // F/(k_B T)=-ln(Z)
  stats.td.S  = E/Z+log(Z);        // S/k_B=beta<H>+ln(Z)
  Sym->calculate_TD(step, diag, stats, rescale_factor);  // symmetry-specific calculation routine
  stats.td.save_values();
}

template<typename S>
void calculate_spectral_and_expv(const Step &step, Stats<S> &stats, Output<S> &output, Oprecalc<S> &oprecalc,
                                 const DiagInfo<S> &diag, const IterInfo<S> &iterinfo, const Store<S> &store,
                                 std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P) {
  // Zft is used in the spectral function calculations using the conventional approach. We calculate it here, in
  // order to avoid recalculations later on.
  stats.Zft = grand_canonical_Z(step, diag, Sym);
  if (std::string(P.specgt) != "" || std::string(P.speci1t) != "" || std::string(P.speci2t) != "")
    stats.Zgt = grand_canonical_Z(step, diag, Sym, 1.0/(P.gtp*step.scT()) ); // exp(-x*gtp)
  if (std::string(P.specchit) != "") 
    stats.Zchit = grand_canonical_Z(step, diag, Sym, 1.0/(P.chitp*step.scT()) ); // exp(-x*chitp)
  DensMatElements<S> rho, rhoFDM;
  if (step.dmnrg()) {
    if (P.need_rho()) {
      rho.load(step.ndx(), P, fn_rho, P.removefiles);
      check_trace_rho(rho, Sym); // Check if Tr[rho]=1, i.e. the normalization
    }
    if (P.need_rhoFDM())
      rhoFDM.load(step.ndx(), P, fn_rhoFDM, P.removefiles);
  }
  oprecalc.sl.calc(step, diag, rho, rhoFDM, stats, Sym, mt, P);
  if (step.nrg()) {
    measure_singlet(step, stats, diag, iterinfo, output, Sym, P);
    iterinfo.dump_diagonal(P.dumpdiagonal);
  }
  if (step.dmnrg() && P.fdmexpv && step.N() == P.fdmexpvn) measure_singlet_fdm(step, stats, diag, iterinfo, output, rhoFDM, store, Sym, P);
}

// Perform calculations of physical quantities. Called prior to NRG iteration (if calc0=true) and after each NRG
// step.
template<typename S>
void perform_basic_measurements(const Step &step, const DiagInfo<S> &diag, std::shared_ptr<Symmetry<S>> Sym,
                                Stats<S> &stats, Output<S> &output) {
  output.dump_all_energies(diag, step.ndx());
  calculate_TD(step, diag, stats, output, Sym);
  output.annotated.dump(step, diag, stats, Sym);
}

} // namespace

#endif
