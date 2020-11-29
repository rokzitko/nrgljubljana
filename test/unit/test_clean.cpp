#include <string>
#include <sstream>
#include <gtest/gtest.h>

// Reproduce test/c++/test0_clean calculation

#include "test_common.hpp"
#include "test_clean.hpp"
#include <basicio.hpp>
#include <core.hpp>

using namespace NRG;

TEST(Clean, H) { // NOLINT
  Params P;
  setup_P_clean(P);

  auto SymSP = setup_Sym<double>(P); // get the shared pointer
  auto Sym = SymSP.get(); // get the raw pointer
  Stats<double> stats(P, Sym->get_td_fields(), 0.0);
  Step step{P, RUNTYPE::NRG};
  EXPECT_EQ(step.ndx(), 0);
  Store<double> store(P.Ninit,P.Nlen);
  auto diagprev = setup_diag_clean(P, Sym);
  auto operators = setup_operators_clean<double>(operators, diagprev);
  operators.opch = setup_opch_clean(P, Sym, diagprev);
  auto coef = setup_coef_clean<double>(P);
  EXPECT_EQ(coef.xi.max(0), P.Nmax);
  auto output = Output(step.get_runtype(), operators, stats, P);

  SubspaceStructure substruct{diagprev, Sym};
  TaskList tasklist{substruct};
  DiagInfo<double> diag;
  for(const auto &I : tasklist.get()) {
    auto h = hamiltonian<double>(step, I, operators.opch, coef, diagprev, output, Sym, P);
    dump_matrix(h);
    auto e = diagonalise(h, DiagParams(P, 1.0), -1); // -1 = not using MPI
    std::cout << "eig(" << I << ")=" << e.value_orig << std::endl;
    diag[I] = e;
  }

  stats.Egs = diag.Egs_subtraction();
  stats.update(step);
  truncate_prepare(step, diag, Sym->multfnc(), P);
  calc_abs_energies(step, diag, stats);
  calculate_TD(step, diag, stats, Sym);
  split_in_blocks(diag, substruct);
  MemTime mt;
  auto oprecalc = Oprecalc<double>(step.get_runtype(), operators, SymSP, mt, P);
  oprecalc.recalculate_operators(operators, step, diag, substruct);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag, operators, store, Sym->multfnc(), mt, P);
  diag.truncate_perform();
  EXPECT_EQ(step.last(), true);
  store[step.ndx()] = Subs(diag, substruct, step.last());
  recalc_irreducible(step, diag, substruct, operators.opch, Sym, mt, P);
  operators.opch.dump();
//X  operators.trim_matrices(diag);
//X  diag.clear_eigenvectors();

  mt.brief_report();
  EXPECT_EQ(step.ndx(), step.lastndx());
  stats.GS_energy = stats.total_energy;
  store.shift_abs_energies(stats.GS_energy); // A

  auto rho = init_rho(step, diag, Sym->multfnc());
  rho.save(step.lastndx(), P, fn_rho);
  calc_densitymatrix(rho, store, Sym, mt, P);

  calc_ZnD(store, stats, Sym, P.T);
  fdm_thermodynamics(store, stats, Sym, P.T);
  auto rhoFDM = init_rho_FDM(step.lastndx(), store, stats, Sym->multfnc(), P.T);
  rhoFDM.save(step.lastndx(), P, fn_rhoFDM);
  calc_fulldensitymatrix(step, rhoFDM, store, stats, Sym, mt, P);

  // single-step calculation: no need to recalculate diag & operators, but we need to run
  // subtract_GS_energy()
  diag.subtract_GS_energy(stats.GS_energy); // B

  Step step_dmnrg{P, RUNTYPE::DMNRG};
  auto oprecalc_dmnrg = Oprecalc<double>(step_dmnrg.get_runtype(), operators, SymSP, mt, P);
  auto output_dmnrg = Output(step_dmnrg.get_runtype(), operators, stats, P);
  calculate_spectral_and_expv(step_dmnrg, stats, output_dmnrg, oprecalc_dmnrg, diag, operators, store, Sym->multfnc(), mt, P);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
