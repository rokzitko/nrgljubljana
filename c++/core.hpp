#ifndef _core_hpp_
#define _core_hpp_

#include <cstddef>
#include <set>
#include <algorithm>
#include <range/v3/all.hpp>

#include "constants.hpp"
#include "traits.hpp"
#include "invar.hpp"
#include "eigen.hpp"
#include "operators.hpp"
#include "subspaces.hpp"
#include "store.hpp"
#include "step.hpp"
#include "stats.hpp"
#include "spectral.hpp"
#include "coef.hpp"
#include "tridiag.hpp"
#include "diag.hpp"
#include "symmetry.hpp"
#include "matrix.hpp"
#include "recalc.hpp"
#include "read-input.hpp"
#include "spectrum.hpp"
#include "algo.hpp"
#include "dmnrg.hpp"
#include "splitting.hpp"
#include "output.hpp"
#include "oprecalc.hpp"
#include "measurements.hpp"
#include "truncation.hpp"
#include "diag_mpi.hpp"
#include "diag_openmp.hpp"
#include "diag_serial.hpp"
#include "h5.hpp"
#include "io.hpp"

#include <fmt/format.h>
#include <fmt/color.h>

namespace NRG {

// Determine the ranges of index r
template<scalar S>
SubspaceDimensions::SubspaceDimensions(const Invar &I, const InvarVec &ancestors, const DiagInfo<S> &diagprev,
                                       const Symmetry<S> *Sym, const bool ignore_inequality) : ancestors(ancestors) {
  for (const auto &[i, anc] : ancestors | ranges::views::enumerate) {
    const bool coupled = Sym->triangle_inequality(I, anc, Sym->QN_subspace(i));
    dims.push_back(coupled || ignore_inequality? diagprev.size_subspace(anc) : 0);
  }
  // The triangle inequality test here is *required*. There are cases where a candidate subspace exists (as generated
  // from the In vector as one of the "combinations"), but it is actually decoupled from space I, because the
  // triangle inequality is not satisfied. [Set ignore_inequality=true to disable the check for testing purposes.]
}

// Determine the structure of matrices in the new NRG shell
template<scalar S>
SubspaceStructure::SubspaceStructure(const DiagInfo<S> &diagprev, const Symmetry<S> *Sym) {
  for (const auto &I : new_subspaces(diagprev, Sym))
    (*this)[I] = SubspaceDimensions{I, Sym->ancestors(I), diagprev, Sym};
}

// Subspaces for the new iteration
template<scalar S>
auto new_subspaces(const DiagInfo<S> &diagprev, const Symmetry<S> *Sym) {
  std::set<Invar> subspaces;
  for (const auto &I : diagprev.subspaces()) {
    const auto all = Sym->new_subspaces(I);
    const auto non_empty = all | ranges::views::filter([Sym](const auto &In) { return Sym->Invar_allowed(In); }) | ranges::to<std::vector>();
    std::copy(non_empty.begin(), non_empty.end(), std::inserter(subspaces, subspaces.end()));
  }
  return subspaces;
}

template<scalar S>
auto hamiltonian(const Step &step, const Invar &I, const Opch<S> &opch, const Coef<S> &coef,
                 const DiagInfo<S> &diagprev, const Output<S> &output, const Symmetry<S> *Sym, const Params &P) {
  const auto anc = Sym->ancestors(I);
  const SubspaceDimensions rm{I, anc, diagprev, Sym};
  const auto dim = rm.total();
  nrglog('i', std::endl << "Subspace (" << I << ") dim=" << dim); // skip a line
  if (P.logletter('s')) {
    std::cout << "Ancestors of (" << I << "): ";
    rm.info();
  }
  auto h = zero_matrix<S>(dim);
  for (const auto i : Sym->combs()) {
    const auto range = rm.view(i);
    for (const auto & [n, r] : range | ranges::views::enumerate)
      h(r,r) = P.nrg_step_scale_factor() * diagprev.at(anc[i]).values.corr(n); // H_{N+1}=\lambda^{1/2} H_N+\xi_N (hopping terms)
  }
  Sym->make_matrix(h, step, rm, I, anc, opch, coef);  // Symmetry-type-specific matrix initialization steps
  if (P.logletter('m')) dump_matrix(h);
  if (P.h5raw && (P.h5all || (P.h5last && step.last())) && P.h5ham)
    h5_dump_matrix(*output.h5raw, std::to_string(step.ndx()+1) + "/hamiltonian/" + I.name() + "/matrix", h);
  return h;
}

template<scalar S>
auto load_or_compute_diag(const Step &step, const Operators<S> &operators, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diagprev,
                          Output<S> &output, const TaskList &tasklist, const Symmetry<S> *Sym, DiagEngine<S> *eng, MemTime &mt,
                          const Params &P, const double diagratio) {
  if (step.nrg()) {
    if (!(P.resume && P.laststored.has_value() && step.ndx() <= P.laststored.value())) {
      const auto section_timing = mt.time_it("diag");
      return eng->diagonalisations(step, operators.opch, coef, diagprev, output, tasklist.get(), DiagParams(P, diagratio), Sym, P);
    }
    return DiagInfo<S>(step.ndx(), P, false);
  }
  auto diag = DiagInfo<S>(step.ndx(), P, P.removefiles);
  diag.subtract_GS_energy(stats.GS_energy);
  return diag;
}

template<scalar S>
void initialize_diag_energy_reference(const Params &P, Stats<S> &stats, DiagInfo<S> &diag) {
  if (P.floquet)
    stats.Egs = 0.0;
  else
    stats.Egs = diag.find_Egs();
  diag.set_shift_Egs(stats.Egs);
}

template<scalar S>
void prepare_diag_for_truncation(const Step &step, DiagInfo<S> &diag, const Symmetry<S> *Sym, const Params &P) {
  Clusters<S> clusters(diag, P.fixeps, P);
  truncate_prepare(step, diag, Sym->multfnc(), P);
}

template<scalar S>
auto do_diag(const Step &step, const Operators<S> &operators, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diagprev,
             Output<S> &output, const TaskList &tasklist, const Symmetry<S> *Sym, DiagEngine<S> *eng, MemTime &mt, const Params &P) {
  step.infostring();
  Sym->show_coefficients(step, coef);
  double diagratio = P.diagratio; // non-const
  DiagInfo<S> diag;
  while (true) {
    try {
      diag = load_or_compute_diag(step, operators, coef, stats, diagprev, output, tasklist, Sym, eng, mt, P, diagratio);
      initialize_diag_energy_reference(P, stats, diag);
      prepare_diag_for_truncation(step, diag, Sym, P);
      break;
    }
    catch (NotEnough &e) {
      color_print(P.pretty_out, fmt::emphasis::bold | fg(fmt::color::yellow), "Insufficient number of states computed.\n");
      if (!(step.nrg() && P.restart)) break;
      diagratio = std::min(diagratio * P.restartfactor, 1.0);
      color_print(P.pretty_out, fmt::emphasis::bold | fg(fmt::color::yellow), "\nRestarting this iteration step. diagratio={}\n\n", diagratio);
    }
  }
  return diag;
}

// Operator sumrules
template<scalar S, typename F> 
auto norm(const MatrixElements<S> &m, const Symmetry<S> *Sym, F factor_fnc, const int SPIN) {
  weight_traits<S> sum{};
  for (const auto &[II, mat] : m) {
    const auto & [I1, Ip] = II;
    if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
    sum += factor_fnc(Ip, I1) * frobenius_norm(mat);
  }
  return 2.0 * sum.real(); // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
}

template<scalar S>
void operator_sumrules(const Operators<S> &a, const Symmetry<S> *Sym) {
  // We check sum rules wrt some given spin (+1/2, by convention). For non-spin-polarized calculations, this is
  // irrelevant (0).
  const int SPIN = Sym->isfield() ? 1 : 0;
  for (const auto &[name, m] : a.opd)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensFactorFnc(), SPIN) << std::endl;
  for (const auto &[name, m] : a.opq)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensquadFactorFnc(), 0) << std::endl;
  for (const auto && [i, ch] : a.opch | ranges::views::enumerate)
    for (const auto && [j, m] : ch | ranges::views::enumerate)
      std::cout << "norm[f," << i << "," << j << "]=" << norm(m, Sym, Sym->SpecdensFactorFnc(), SPIN) << std::endl;
}

// Perform processing after a successful NRG step.
template<scalar S>
void handle_floquet_postprocess(const Step &step, Operators<S> &operators, Stats<S> &stats, DiagInfo<S> &diag,
                                const SubspaceStructure &substruct, Oprecalc<S> &oprecalc, const Params &P) {
  split_in_blocks(diag, substruct, false); // false = don't discard
  my_assert(P.extra_params.count("Omega") > 0);
  const auto scale = step.scale();
  const auto _Omega = std::stod(P.extra_params.at("Omega"));
  const auto Omega = _Omega / scale;
  nrglog('0' , "Omega=" << _Omega << " rescaled=" << Omega);
  auto mnew = oprecalc.recalculate_operator_m(operators, step, diag, substruct, P);
  if (P.logletter('1'))
    dump_diagonal_op("m", mnew, 0);
  double emin = std::numeric_limits<double>::max();
  for (auto &[I, eig] : diag) {
    if (P.logletter('2') || P.logletter('3'))
      std::cout << "Floquet: subspace " << I << std::endl;
    for (size_t i = 0; i < eig.values.size(); i++) {
      const auto e = eig.values.raw(i);
      emin = std::min(emin, e);
      const auto m = mnew[Twoinvar(I, I)](i, i);
      const auto mOmega = std::real(m) * Omega;
      const auto e0 = e - mOmega;
      const auto x = e0 + abs(mOmega); // second term: penalize high-m states
      nrglog('2', "i=" << i << " e=" << e << " m=" << m << " e0=e-m*Omega=" << e0 << " x=" << x);
      nrglog('3', "rs i=" << i << " e=" << e*scale << " m=" << m << " e0=e-m*Omega=" << e0*scale << " x=" << x*scale);
      eig.values.set_crit(i, x);
    }
  }
  stats.Egs = emin;
  std::cout << "Egs=" << stats.Egs << std::endl;
  const auto Clw = diag.find_Clw();
  std::cout << "Clw=" << Clw << std::endl;
  diag.shift(stats.Egs, Clw);
  if (P.logletter('4'))
    diag.report(true);
  diag.sort_by_c();
  split_in_blocks(diag, substruct, true);
}

template<scalar S>
void finalize_nrg_iteration_metadata(const Step &step, Stats<S> &stats, DiagInfo<S> &diag, const Params &P) {
  stats.update(step);
  for (auto &eig : diag.eigs()) {
    eig.values.set_scale(step.scale());
    eig.values.set_T_shift(stats.total_energy);
  }
}

template<scalar S>
void persist_nrg_iteration_outputs(const Step &step, Stats<S> &stats, DiagInfo<S> &diag, Output<S> &output,
                                   const SubspaceStructure &substruct, const Symmetry<S> *Sym, const Params &P) {
  if (P.dm && !(P.resume && P.laststored.has_value() && step.ndx() <= P.laststored.value()))
    diag.save(step.ndx(), P);
  perform_basic_measurements(step, diag, Sym, stats, output, P);
  if (P.h5raw && (P.h5all || (P.h5last && step.last()))) {
    diag.h5save(*output.h5raw, std::to_string(step.ndx()+1) + "/eigen/", P.h5vectors);
    if (P.h5U)
      h5save_blocks(*output.h5raw, std::to_string(step.ndx()+1) + "/U/", diag, substruct);
  }
}

template<scalar S>
void archive_iteration_state(const Step &step, const DiagInfo<S> &diag, const SubspaceStructure &substruct,
                             ThermoStore<S> &store, BackiterStore &store_all, const Symmetry<S> *Sym, const Params &P) {
  store_all[step.ndx()] = BackiterSubs(diag, substruct);
  if (P.project == ""s) {
    store[step.ndx()] = ThermoSubs<S>(diag, step.last());
  } else {
    const auto projected = Sym->project(diag, P.project);
    store[step.ndx()] = ThermoSubs<S>(projected, step.last());
  }
}

template<scalar S>
void recalculate_and_measure(const Step &step, Operators<S> &operators, Stats<S> &stats, DiagInfo<S> &diag,
                             Oprecalc<S> &oprecalc, Output<S> &output, const SubspaceStructure &substruct,
                             const BackiterStore &store_all, const Symmetry<S> *Sym, MemTime &mt, const Params &P) {
  oprecalc.recalculate_operators(operators, step, diag, substruct, P);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag, operators, store_all, mt, Sym, P);
}

template<scalar S>
void after_diag(const Step &step, Operators<S> &operators, Stats<S> &stats, DiagInfo<S> &diag, Output<S> &output,
                const SubspaceStructure &substruct, ThermoStore<S> &store, BackiterStore &store_all, Oprecalc<S> &oprecalc, const Symmetry<S> *Sym,
                MemTime &mt, const Params &P) {
  nrglog('@', "after_diag()");
  if (step.nrg()) {
    if (P.floquet) {
      handle_floquet_postprocess(step, operators, stats, diag, substruct, oprecalc, P);
    } else {
      split_in_blocks(diag, substruct, true); // true = discard
    }
    finalize_nrg_iteration_metadata(step, stats, diag, P);
    persist_nrg_iteration_outputs(step, stats, diag, output, substruct, Sym, P);
  }
  if (P.do_recalc_all(step.get_runtype())) { // Either ...
    recalculate_and_measure(step, operators, stats, diag, oprecalc, output, substruct, store_all, Sym, mt, P);
  }
  nrglog('@', "truncate_perform()");
  diag.truncate_perform();                               // Actual truncation occurs at this point
  archive_iteration_state(step, diag, substruct, store, store_all, Sym, P);
  if (!step.last()) {
    nrglog('@', "recalc_irreducible()");
    recalc_irreducible(step, diag, substruct, operators.opch, Sym, mt, P);
    if (P.dump_f) operators.opch.dump();
  }
  if (P.do_recalc_kept(step.get_runtype())) { // ... or ...
    recalculate_and_measure(step, operators, stats, diag, oprecalc, output, substruct, store_all, Sym, mt, P);
  }
  if (P.checksumrules) operator_sumrules(operators, Sym);
  if (P.h5raw && (P.h5all || (P.h5last && step.last())) && P.h5ops)
    operators.h5save(*output.h5raw, std::to_string(step.ndx()+1));
}

// Perform one iteration step
template<scalar S>
auto iterate(const Step &step, Operators<S> &operators, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diagprev,
             Output<S> &output, ThermoStore<S> &store, BackiterStore &store_all, Oprecalc<S> &oprecalc, const Symmetry<S> *Sym, DiagEngine<S> *eng, MemTime &mt, const Params &P) {
  SubspaceStructure substruct{diagprev, Sym};
  TaskList tasklist{substruct, !P.silent}; // verbose = !P.silent
  if (P.h5raw && (P.h5all || (P.h5last && step.last())) && P.h5struct)
    substruct.h5save(*output.h5raw, std::to_string(step.ndx()+1) + "/structure");
  auto diag = do_diag(step, operators, coef, stats, diagprev, output, tasklist, Sym, eng, mt, P);
  after_diag(step, operators, stats, diag, output, substruct, store, store_all, oprecalc, Sym, mt, P);
  operators.trim_matrices(diag);
  diag.clear_eigenvectors();
  mt.brief_report();
  return diag;
}

// Perform calculations with quantities from 'data' file
template<scalar S>
void docalc0(Step &step, const Operators<S> &operators, const DiagInfo<S> &diag0, Stats<S> &stats, Output<S> &output,
             Oprecalc<S> &oprecalc, const Symmetry<S> *Sym, MemTime &mt, const Params &P) {
  step.set(P.Ninit - 1); // in the usual case with Ninit=0, this will result in N=-1
  std::cout << std::endl << "Before NRG iteration";
  std::cout << " (N=" << step.N() << ")" << std::endl;
  perform_basic_measurements(step, diag0, Sym, stats, output, P);
  BackiterStore empty_st(0, 0);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag0, operators, empty_st, mt, Sym, P);
  if (P.checksumrules) operator_sumrules(operators, Sym);
}

template<scalar S>
auto nrg_loop(Step &step, Operators<S> &operators, const Coef<S> &coef, Stats<S> &stats, DiagInfo<S> diag,
              Output<S> &output, ThermoStore<S> &store, BackiterStore &store_all, Oprecalc<S> &oprecalc, const Symmetry<S> *Sym, DiagEngine<S> *eng, MemTime &mt, const Params &P) {
  for (step.init(); !step.end(); step.next())
    diag = iterate(step, operators, coef, stats, diag, output, store, store_all, oprecalc, Sym, eng, mt, P);
  step.set(step.lastndx());
  return diag;
}

} // namespace

#endif
