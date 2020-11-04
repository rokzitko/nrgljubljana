#ifndef _core_hpp_
#define _core_hpp_

#include <cstddef>
#include <set>
#include <algorithm>
#include <omp.h>
#include <h5cpp/all>
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
#include "mpi_diag.hpp"

namespace NRG {

// Determine the ranges of index r
template<typename S>
SubspaceDimensions::SubspaceDimensions(const Invar &I, const InvarVec &ancestors, const DiagInfo<S> &diagprev, 
                                       std::shared_ptr<Symmetry<S>> Sym) : ancestors(ancestors) {
  for (const auto &[i, anc] : ancestors | ranges::views::enumerate)
    dims.push_back(Sym->triangle_inequality(I, anc, Sym->QN_subspace(i)) ? diagprev.size_subspace(anc) : 0);
  // The triangle inequality test here is *required*. There are cases where a candidate subspace exists (as generated
  // from the In vector as one of the "combinations"), but it is actually decoupled, because the triangle inequality
  // is not satisfied.
}

// Determine the structure of matrices in the new NRG shell
template<typename S>
SubspaceStructure::SubspaceStructure(const DiagInfo<S> &diagprev, std::shared_ptr<Symmetry<S>> Sym) {
  for (const auto &I : new_subspaces(diagprev, Sym))
    (*this)[I] = SubspaceDimensions{I, Sym->ancestors(I), diagprev, Sym};
}

// Subspaces for the new iteration
template<typename S>
auto new_subspaces(const DiagInfo<S> &diagprev, std::shared_ptr<Symmetry<S>> Sym) {
  std::set<Invar> subspaces;
  for (const auto &I : diagprev.subspaces()) {
    const auto all = Sym->new_subspaces(I);
    const auto non_empty = all | ranges::views::filter([&Sym](const auto &In) { return Sym->Invar_allowed(In); }) | ranges::to<std::vector>();
    std::copy(non_empty.begin(), non_empty.end(), std::inserter(subspaces, subspaces.end()));
  }
  return subspaces;
}

template<typename S>
typename traits<S>::Matrix hamiltonian(const Step &step, const Invar &I, const Opch<S> &opch, const Coef<S> &coef, 
                                       const DiagInfo<S> &diagprev, const Output<S> &output, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  const auto anc = Sym->ancestors(I);
  const SubspaceDimensions rm{I, anc, diagprev, Sym};
  auto h = Zero_matrix<S>(rm.total());
  for (const auto i : Sym->combs()) {
    const auto range = rm.view(i);
    for (const auto & [n, r] : range | ranges::views::enumerate)
      h(r,r) = P.nrg_step_scale_factor() * diagprev.at(anc[i]).value_zero(n); // H_{N+1}=\lambda^{1/2} H_N+\xi_N (hopping terms)
  }
  Sym->make_matrix(h, step, rm, I, anc, opch, coef);  // Symmetry-type-specific matrix initialization steps
  if (P.logletter('m')) dump_matrix(h);
  if (output.h5raw)
    h5::write(output.h5raw.value(), std::to_string(step.ndx()+1) + "/hamiltonian/" + I.name() + "/matrix", h);
  return h;
}

template<typename S>
auto diagonalisations_OpenMP(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, const Output<S> &output,
                             const std::vector<Invar> &tasks, const DiagParams &DP, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  DiagInfo<S> diagnew;
  const auto nr = tasks.size();
  size_t itask = 0;
  // cppcheck-suppress unreadVariable symbolName=nth
  const int nth = P.diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
  for (itask = 0; itask < nr; itask++) {
    const Invar I  = tasks[itask];
    auto h = hamiltonian(step, I, opch, coef, diagprev, output, Sym, P); // non-const, consumed by diagonalise()
    const int thid = omp_get_thread_num();
#pragma omp critical
    { nrglog('(', "Diagonalizing " << I << " size=" << h.size1() << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
    auto e = diagonalise<S>(h, DP, -1); // -1 = not using MPI
#pragma omp critical
    { diagnew[I] = e; }
  }
  return diagnew;
}

// Build matrix H(ri;r'i') in each subspace and diagonalize it
template<typename S>
auto diagonalisations(const Step &step, const Opch<S> &opch, const Coef<S> &coef, const DiagInfo<S> &diagprev, 
                      const Output<S> &output, const std::vector<Invar> &tasks, const double diagratio, 
                      std::shared_ptr<Symmetry<S>> Sym, MPI_diag &mpi, MemTime &mt, const Params &P) {
  const auto section_timing = mt.time_it("diag");
  return P.diag_mode == "MPI" ? mpi.diagonalisations_MPI<S>(step, opch, coef, diagprev, output, tasks, DiagParams(P, diagratio), Sym, P) 
                              : diagonalisations_OpenMP(step, opch, coef, diagprev, output, tasks, DiagParams(P, diagratio), Sym, P);
}

template<typename S>
auto do_diag(const Step &step, IterInfo<S> &iterinfo, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diagprev, const Output<S> &output,
             SubspaceStructure &substruct, std::shared_ptr<Symmetry<S>> Sym, MPI_diag &mpi, MemTime &mt, const Params &P) {
  step.infostring();
  Sym->show_coefficients(step, coef);
  auto tasks = substruct.task_list();
  double diagratio = P.diagratio; // non-const
  DiagInfo<S> diag;
  while (true) {
    try {
      if (step.nrg()) {
        if (!(P.resume && int(step.ndx()) <= P.laststored))
          diag = diagonalisations(step, iterinfo.opch, coef, diagprev, output, tasks, diagratio, Sym, mpi, mt, P); // compute in first run
        else
          diag = DiagInfo<S>(step.ndx(), P, false); // or read from disk
      }
      if (step.dmnrg()) {
        diag = DiagInfo<S>(step.ndx(), P, P.removefiles); // read from disk in second run
        diag.subtract_GS_energy(stats.GS_energy);
      }
      stats.Egs = diag.find_groundstate();
      if (step.nrg()) // should be done only once!
        diag.subtract_Egs(stats.Egs);
      Clusters<S> clusters(diag, P.fixeps);
      truncate_prepare(step, diag, Sym, P);
      break;
    }
    catch (NotEnough &e) {
      fmt::print(fmt::emphasis::bold | fg(fmt::color::yellow), "Insufficient number of states computed.\n");
      if (!(step.nrg() && P.restart)) break;
      diagratio = std::min(diagratio * P.restartfactor, 1.0);
      fmt::print(fmt::emphasis::bold | fg(fmt::color::yellow), "\nRestarting this iteration step. diagratio={}\n\n", diagratio);
    }
  }
  return diag;
}

// Absolute energies. Must be called in the first NRG run after stats.total_energy has been updated, but before
// store_transformations(). absenergyG is updated to its correct values (referrenced to absolute 0) in
// shift_abs_energies().
template<typename S>
void calc_abs_energies(const Step &step, DiagInfo<S> &diag, const Stats<S> &stats) {
  for (auto &eig : diag.eigs()) {
    eig.absenergy_zero = eig.value_zero * step.scale();    // referenced to the lowest energy in current NRG step (not modified later on)
    eig.absenergy = eig.absenergy_zero;                    // absolute energies (not modified later on)
    std::transform(eig.absenergy.begin(), eig.absenergy.end(), eig.absenergy.begin(), [v = stats.total_energy](auto x) { return x + v; });
    eig.absenergyG = eig.absenergy;                        // referenced to the absolute 0 (updated by shft_abs_energies())
  }
}

// Operator sumrules
template<typename S, typename F> 
auto norm(const MatrixElements<S> &m, std::shared_ptr<Symmetry<S>> Sym, F factor_fnc, const int SPIN) {
  typename traits<S>::t_weight sum{};
  for (const auto &[II, mat] : m) {
    const auto & [I1, Ip] = II;
    if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
    sum += factor_fnc(Ip, I1) * frobenius_norm(mat);
  }
  return 2.0 * sum.real(); // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
}

template<typename S>
void operator_sumrules(const IterInfo<S> &a, std::shared_ptr<Symmetry<S>> Sym) {
  // We check sum rules wrt some given spin (+1/2, by convention). For non-spin-polarized calculations, this is
  // irrelevant (0).
  const int SPIN = Sym->isfield() ? 1 : 0;
  for (const auto &[name, m] : a.opd)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensFactorFnc(), SPIN) << std::endl;
  for (const auto &[name, m] : a.opq)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensquadFactorFnc(), 0) << std::endl;
}

// Perform processing after a successful NRG step. Also called from doZBW() as a final step.
template<typename S>
void after_diag(const Step &step, IterInfo<S> &iterinfo, Stats<S> &stats, DiagInfo<S> &diag, Output<S> &output,
                SubspaceStructure &substruct, Store<S> &store, Oprecalc<S> &oprecalc, std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P) {
  stats.total_energy += stats.Egs * step.scale(); // stats.Egs has already been initialized
  std::cout << "Total energy=" << HIGHPREC(stats.total_energy) << "  Egs=" << HIGHPREC(stats.Egs) << std::endl;
  stats.rel_Egs[step.ndx()] = stats.Egs;
  stats.abs_Egs[step.ndx()] = stats.Egs * step.scale();
  stats.energy_offsets[step.ndx()] = stats.total_energy;
  if (step.nrg()) {
    calc_abs_energies(step, diag, stats);  // only in the first run, in the second one the data is loaded from file!
    if (P.dm && !(P.resume && int(step.ndx()) <= P.laststored))
      diag.save(step.ndx(), P);
    perform_basic_measurements(step, diag, Sym, stats, output); // Measurements are performed before the truncation!
  }
  if (output.h5raw)
    diag.h5save(output.h5raw.value(), std::to_string(step.ndx()+1) + "/eigen/", step.dmnrg());
  if (!P.ZBW) {
    split_in_blocks(diag, substruct);
    if (output.h5raw)
      h5save_blocks(output.h5raw.value(), std::to_string(step.ndx()+1) + "/U/", diag, substruct);
  }
  if (P.do_recalc_all(step.get_runtype())) { // Either ...
    oprecalc.recalculate_operators(iterinfo, step, diag, substruct);
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, store, Sym, mt, P);
  }
  if (!P.ZBW)
    diag.truncate_perform();                        // Actual truncation occurs at this point
  store[step.ndx()] = Subs(diag, substruct, step.last());  // Store information about subspaces and states for DM algorithms
  if (!step.last()) {
    recalc_irreducible(step, diag, substruct, iterinfo.opch, Sym, mt, P);
    if (P.dump_f) iterinfo.opch.dump();
  }
  if (P.do_recalc_kept(step.get_runtype())) { // ... or ...
    oprecalc.recalculate_operators(iterinfo, step, diag, substruct);
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, store, Sym, mt, P);
  }
  if (P.do_recalc_none())  // ... or this
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, store, Sym, mt, P);
  if (P.checksumrules) operator_sumrules(iterinfo, Sym);
  if (output.h5raw)
    iterinfo.h5save(output.h5raw.value(), std::to_string(step.ndx()+1));
}

// Perform one iteration step
template<typename S>
auto iterate(const Step &step, IterInfo<S> &iterinfo, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diagprev,
             Output<S> &output, Store<S> &store, Oprecalc<S> &oprecalc, std::shared_ptr<Symmetry<S>> Sym, MPI_diag &mpi, MemTime &mt, const Params &P) {
  SubspaceStructure substruct{diagprev, Sym};
  if (output.h5raw)
    substruct.h5save(output.h5raw.value(), std::to_string(step.ndx()+1) + "/structure");
  auto diag = do_diag(step, iterinfo, coef, stats, diagprev, output, substruct, Sym, mpi, mt, P);
  after_diag(step, iterinfo, stats, diag, output, substruct, store, oprecalc, Sym, mt, P);
  iterinfo.trim_matrices(diag);
  diag.clear_eigenvectors();
  mt.brief_report();
  return diag;
}

// Perform calculations with quantities from 'data' file
template<typename S>
void docalc0(Step &step, const IterInfo<S> &iterinfo, const DiagInfo<S> &diag0, Stats<S> &stats, Output<S> &output, 
             Oprecalc<S> &oprecalc, std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P) {
  step.set(P.Ninit - 1); // in the usual case with Ninit=0, this will result in N=-1
  std::cout << std::endl << "Before NRG iteration";
  std::cout << " (N=" << step.N() << ")" << std::endl;
  perform_basic_measurements(step, diag0, Sym, stats, output);
  Store<S> empty_st(0, 0);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag0, iterinfo, empty_st, Sym, mt, P);
  if (P.checksumrules) operator_sumrules(iterinfo, Sym);
}

// doZBW() takes the place of iterate() called from main_loop() in the case of zero-bandwidth calculation.
// It replaces do_diag() and calls after_diag() as the last step.
template<typename S>
auto nrg_ZBW(Step &step, IterInfo<S> &iterinfo, Stats<S> &stats, const DiagInfo<S> &diag0, Output<S> &output, 
             Store<S> &store, Oprecalc<S> &oprecalc, std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P) {
  std::cout << std::endl << "Zero bandwidth calculation" << std::endl;
  step.set_ZBW();
  // --- begin do_diag() equivalent
  DiagInfo<S> diag;
  if (step.nrg())
    diag = diag0;
  if (step.dmnrg()) {
    diag = DiagInfo<S>(step.ndx(), P, P.removefiles);
    diag.subtract_GS_energy(stats.GS_energy);
  }
  stats.Egs = diag.find_groundstate();
  if (step.nrg())      
    diag.subtract_Egs(stats.Egs);
  truncate_prepare(step, diag, Sym, P); // determine # of kept and discarded states
  // --- end do_diag() equivalent
  SubspaceStructure substruct{};
  after_diag(step, iterinfo, stats, diag, output, substruct, store, oprecalc, Sym, mt, P);
  return diag;
}

template<typename S>
auto nrg_loop(Step &step, IterInfo<S> &iterinfo, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diag0,
              Output<S> &output, Store<S> &store, Oprecalc<S> &oprecalc, std::shared_ptr<Symmetry<S>> Sym, MPI_diag &mpi, MemTime &mt, const Params &P) {
  auto diag = diag0;
  for (step.init(); !step.end(); step.next())
    diag = iterate(step, iterinfo, coef, stats, diag, output, store, oprecalc, Sym, mpi, mt, P);
  step.set(step.lastndx());
  return diag;
}

} // namespace

#endif
