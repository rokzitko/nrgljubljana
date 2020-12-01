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

template<typename T, typename S>
auto trace_exp(const ublas::vector<T> &e, const ublas::matrix<S> &m, const double factor) { // Tr[exp(-factor*e) m]
  my_assert(e.size() == m.size1() && e.size() == m.size2());
  return ranges::accumulate(range0(e.size()), S{}, {}, [&e, &m, factor](const auto r){ return exp(-factor * e(r)) * m(r, r); });
}
   
// note: t_expv = t_matel, thus the return type is OK
template<typename S, typename MF, typename t_matel = matel_traits<S>>
CONSTFNC auto calc_trace_singlet(const DiagInfo<S> &diag, const MatrixElements<S> &m, MF mult, const double factor) {
  return ranges::accumulate(diag, S{}, {}, [&m, &mult, factor](const auto &x){
    const auto [I, eig] = x; return mult(I) * trace_exp(eig.value_zero, m.at({I,I}), factor); });
}

template<typename T>
auto sum_of_exp(T values, const double factor)
{
  return ranges::accumulate(values, 0.0, {}, [factor](const auto &x){ return exp(-factor*x); });
}      

// Measure thermodynamic expectation values of singlet operators
template<typename S, typename MF>
void measure_singlet(const double factor, Stats<S> &stats, const Operators<S> &a, MF mult, const DiagInfo<S> &diag) {
  const auto Z = ranges::accumulate(diag, 0.0, {}, [mult, factor](const auto &d) { const auto &[I, eig] = d;
                                                   return mult(I) * sum_of_exp(eig.value_zero, factor); });
  for (const auto &[name, m] : a.ops)  stats.expv[name] = calc_trace_singlet(diag, m, mult, factor) / Z;
  for (const auto &[name, m] : a.opsg) stats.expv[name] = calc_trace_singlet(diag, m, mult, factor) / Z;
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

template<typename S, typename MF, typename t_matel = matel_traits<S>>
CONSTFNC auto calc_trace_fdm_kept(const size_t ndx, const MatrixElements<S> &n, const DensMatElements<S> &rhoFDM,
                                  const Store<S> &store, MF mult) {
  t_matel tr{};
  for (const auto &[I, rhoI] : rhoFDM)
    tr += mult(I) * trace_contract(rhoI, n.at({I,I}), store[ndx].at(I).kept()); // over kept states ONLY
  return tr;
}

template<typename S, typename MF>
void measure_singlet_fdm(const size_t ndx, Stats<S> &stats, const Operators<S> &a, MF mult,
                         const DensMatElements<S> &rhoFDM, const Store<S> &store) {
  for (const auto &[name, m] : a.ops)  stats.fdmexpv[name] = calc_trace_fdm_kept(ndx, m, rhoFDM, store, mult);
  for (const auto &[name, m] : a.opsg) stats.fdmexpv[name] = calc_trace_fdm_kept(ndx, m, rhoFDM, store, mult);
}

// Calculate grand canonical partition function at current NRG energy shell. This is not the same as the true
// partition function of the full problem! Instead this is the Z_N that is used to initialize the density matrix,
// i.e. rho = 1/Z_N \sum_{l} exp{-beta E_l} |l;N> <l;N|.  grand_canonical_Z() is also used to calculate stats.Zft,
// that is used to compute the spectral function with the conventional approach, as well as stats.Zgt for G(T)
// calculations, stats.Zchit for chi(T) calculations.
template<typename S, typename MF>
auto grand_canonical_Z(const double factor, const DiagInfo<S> &diag, MF mult) {
  double ZN{};
  for (const auto &[I, eig]: diag) 
    for (const auto &i : eig.kept()) // sum over all kept states
      ZN += mult(I) * exp(-factor * eig.value_zero(i));
  my_assert(ZN >= 1.0);
  return ZN;
}

// Calculate partial statistical sums, ZnD*, and the grand canonical Z (stats.ZZG), computed with respect to absolute
// energies. calc_ZnD() must be called before the second NRG run.
template<typename S>
void calc_ZnD(const Store<S> &store, Stats<S> &stats, const Symmetry<S> *Sym, const double T) {
  mpf_set_default_prec(400); // this is the number of bits, not decimal digits!
  for (const auto N : store.Nall()) {
    my_mpf ZnDG, ZnDN; // arbitrary-precision accumulators to avoid precision loss
    mpf_set_d(ZnDG, 0.0);
    mpf_set_d(ZnDN, 0.0);
    for (const auto &[I, ds] : store[N])
      for (const auto i : ds.all()) {
        my_mpf g, n;
        mpf_set_d(g, Sym->mult(I) * exp(-ds.eig.absenergyG[i]/T));     // absenergyG >= 0.0
        mpf_set_d(n, Sym->mult(I) * exp(-ds.eig.absenergy_zero[i]/T)); // absenergy_zero >= 0.0
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
void fdm_thermodynamics(const Store<S> &store, Stats<S> &stats, const Symmetry<S> *Sym, const double T)
{
  stats.Z_fdm = stats.ZZG*exp(-stats.GS_energy/T); // this is the true partition function
  stats.F_fdm = -log(stats.ZZG)*T+stats.GS_energy; // F = -k_B*T*log(Z)
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
          mpf_set_d(weight, stats.wn[N] * Sym->mult(I) * exp(-ds.eig.absenergy_zero[i]/T));
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
  stats.E_fdm = mpf_get_d(E);
  my_mpf sqrE;
  mpf_mul(sqrE, E, E);
  my_mpf varE;
  mpf_sub(varE, E2, sqrE);
  stats.C_fdm = mpf_get_d(varE)/pow(T,2);
  stats.S_fdm = (stats.E_fdm-stats.F_fdm)/T;
  std::cout << std::endl;
  std::cout << "Z_fdm=" << HIGHPREC(stats.Z_fdm) << std::endl;
  std::cout << "F_fdm=" << HIGHPREC(stats.F_fdm) << std::endl;
  std::cout << "E_fdm=" << HIGHPREC(stats.E_fdm) << std::endl;
  std::cout << "C_fdm=" << HIGHPREC(stats.C_fdm) << std::endl;
  std::cout << "S_fdm=" << HIGHPREC(stats.S_fdm) << std::endl;
  std::cout << std::endl;
  stats.td_fdm.set("T", T);
  stats.td_fdm.set("F_fdm", stats.F_fdm); 
  stats.td_fdm.set("E_fdm", stats.E_fdm);
  stats.td_fdm.set("C_fdm", stats.C_fdm);
  stats.td_fdm.set("S_fdm", stats.S_fdm);
  stats.td_fdm.save_values();
}

// We calculate thermodynamic quantities before truncation to make better use of the available states. Here we
// compute quantities which are defined for all symmetry types. Other calculations are performed by calculate_TD
// member functions defined in symmetry.h
template<typename S>
void calculate_TD(const Step &step, const DiagInfo<S> &diag, Stats<S> &stats,
                  const Symmetry<S> *Sym, const double additional_factor = 1.0) {
  // Rescale factor for energies. The energies are expressed in units of omega_N, thus we need to appropriately
  // rescale them to calculate the Boltzmann weights at the temperature scale Teff (Teff=scale/betabar).
  const auto rescale_factor = step.TD_factor() * additional_factor;
  const auto Z  = diag.trace([](double x) { return 1; },        rescale_factor, Sym->multfnc()); // partition function
  const auto E  = diag.trace([](double x) { return x; },        rescale_factor, Sym->multfnc()); // Tr[beta H]
  const auto E2 = diag.trace([](double x) { return pow(x,2); }, rescale_factor, Sym->multfnc()); // Tr[(beta H)^2]
  stats.Z = Z;
  stats.td.set("T",     step.Teff());
  stats.td.set("<E>",   E/Z);               // beta <H>
  stats.td.set("<E^2>", E2/Z);              // beta^2 <H^2>
  stats.td.set("C",     E2/Z - pow(E/Z,2)); // C/k_B=beta^2(<H^2>-<H>^2)
  stats.td.set("F",     -log(Z));           // F/(k_B T)=-ln(Z)
  stats.td.set("S",     E/Z+log(Z));        // S/k_B=beta<H>+ln(Z)
  Sym->calculate_TD(step, diag, stats, rescale_factor);  // symmetry-specific calculation routine
  stats.td.save_values();
}

template<typename S, typename MF>
void calc_Z(const Step &step, Stats<S> &stats, const DiagInfo<S> &diag, MF mult, const Params &P) {
  stats.Zft = grand_canonical_Z(step.scT(), diag, mult);
  if (std::string(P.specgt) != "" || std::string(P.speci1t) != "" || std::string(P.speci2t) != "")
    stats.Zgt = grand_canonical_Z(1.0/P.gtp, diag, mult); // exp(-x*gtp)
  if (std::string(P.specchit) != "") 
    stats.Zchit = grand_canonical_Z(1.0/P.chitp, diag, mult); // exp(-x*chitp)
}

template<typename S, typename MF>
void calculate_spectral_and_expv(const Step &step, Stats<S> &stats, Output<S> &output, Oprecalc<S> &oprecalc,
                                 const DiagInfo<S> &diag, const Operators<S> &operators, const Store<S> &store,
                                 MF mult, MemTime &mt, const Params &P) {
  // Calculate the density matrices
  DensMatElements<S> rho, rhoFDM;
  if (step.dmnrg()) {
    if (P.need_rho()) {
      rho.load(step.ndx(), P, fn_rho, P.removefiles);
      check_trace_rho(rho, mult); // Check if Tr[rho]=1, i.e. the normalization
    }
    if (P.need_rhoFDM())
      rhoFDM.load(step.ndx(), P, fn_rhoFDM, P.removefiles);
  }
  // Calculate all spectral functions
  calc_Z(step, stats, diag, mult, P); // required for FT and CFS approaches
  oprecalc.sl.calc(step, diag, rho, rhoFDM, stats, mt, P);
  // Calculate all expectation values
  if (step.nrg()) {
    measure_singlet(step.TD_factor(), stats, operators, mult, diag);
    output.custom->field_values(step.Teff());
    operators.dump_diagonal(P.dumpdiagonal);
  }
  if (step.dmnrg() && P.fdmexpv && step.N() == P.fdmexpvn) {
    measure_singlet_fdm(step.N(), stats, operators, mult, rhoFDM, store);
    output.customfdm->field_values(P.T);
  }
}

// Perform calculations of physical quantities. Called prior to NRG iteration (if calc0=true) and after each NRG
// step.
template<typename S>
void perform_basic_measurements(const Step &step, const DiagInfo<S> &diag, const Symmetry<S> *Sym,
                                Stats<S> &stats, Output<S> &output) {
  output.dump_all_energies(step.ndx(), diag);
  output.annotated.dump(step, diag, stats, Sym->multfnc());
  calculate_TD(step, diag, stats, Sym);
}

} // namespace

#endif
