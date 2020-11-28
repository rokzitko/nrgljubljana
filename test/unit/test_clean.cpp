#include <string>
#include <sstream>
#include <gtest/gtest.h>

// Reproduce test0_clean calculation

#include "test_common.hpp"
#include <basicio.hpp>
#include <core.hpp>

using namespace NRG;

TEST(Clean, H) { // NOLINT
  MemTime mt;
  Params P;
  EXPECT_EQ(P.Lambda.value(), 2.0);
  EXPECT_EQ(P.discretization.value(), "Z"s);
  EXPECT_EQ(P.Ninit.value(), 0);
  P.Nmax = P.Nlen = 1;
  P.keep.setvalue(100);
  P.ops.setvalue("n_f n_f^2");
  P.specs.setvalue("n_f-n_f");
  P.finite.setvalue(true);
  P.finitemats.setvalue(true);
  P.mats.setvalue(10);
  P.T.setvalue(0.1); // temperature

  auto SymSP = setup_Sym<double>(P); // get the shared pointer
  auto Sym = SymSP.get(); // get the raw pointer
  Stats<double> stats(P);
  stats.td.allfields.add(Sym->get_td_fields(), 1);
  stats.total_energy = 0.0; // cf. data file for this case
  Step step{P, RUNTYPE::NRG};
  Store<double> store(P.Ninit,P.Nlen);

  auto diagprev = setup_diag_clean(P, Sym);
  auto n = new_subspaces(diagprev, Sym);
  SubspaceStructure substruct{diagprev, Sym};
  auto tasks = substruct.task_list();

  auto operators = Operators<double>();
  setup_operators_clean(operators, diagprev);
  operators.opch = setup_opch_clean(P, Sym, diagprev);
  auto & opch = operators.opch;
  auto coef = setup_coef_clean<double>(P);
  EXPECT_EQ(coef.xi.max(0), P.Nmax);
  auto output = Output(step.get_runtype(), operators, stats, P);

  DiagInfo<double> diag;
  for(const auto &I : tasks) {
    auto h = hamiltonian<double>(step, I, opch, coef, diagprev, output, Sym, P);
    dump_matrix(h);
    auto e = diagonalise(h, DiagParams(P, 1.0), -1); // -1 = not using MPI
    std::cout << "eig(" << I << ")=" << e.value_orig << std::endl;
    diag[I] = e;
  }

  stats.Egs = diag.find_groundstate();
  diag.subtract_Egs(stats.Egs);
  truncate_prepare(step, diag, Sym->multfnc(), P);
  stats.total_energy += stats.Egs * step.scale();
  std::cout << "Total energy=" << HIGHPREC(stats.total_energy) << "  Egs=" << HIGHPREC(stats.Egs) << std::endl;
  calc_abs_energies(step, diag, stats);
  calculate_TD(step, diag, stats, Sym);
  split_in_blocks(diag, substruct);
  auto oprecalc = Oprecalc<double>(step.get_runtype(), operators, SymSP, mt, P);
  oprecalc.recalculate_operators(operators, step, diag, substruct);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag, operators, store, Sym->multfnc(), mt, P);
  diag.truncate_perform();
  EXPECT_EQ(step.last(), true);
  store[step.ndx()] = Subs(diag, substruct, step.last());
  recalc_irreducible(step, diag, substruct, operators.opch, Sym, mt, P);
  operators.opch.dump(); 
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
