#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <core.hpp>

using namespace NRG;

TEST(Core, H1) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  EXPECT_EQ(P.combs, 4);
  EXPECT_EQ(SymSP->input_subspaces().size(), 4); // In
  EXPECT_EQ(SymSP->nr_combs(), 4);
  auto Sym = SymSP.get(); // get the raw pointer

  Stats<double> stats(P, Sym->get_td_fields(), 0.0);
  Step step{P, RUNTYPE::NRG};
  Store<double> store(0,1);

  auto diag = setup_diag(P, Sym);
  auto n = new_subspaces(diag, Sym);
  std::cout << n << std::endl;
  SubspaceStructure substruct{diag, Sym};
  auto tasks = substruct.task_list();
  auto I = tasks.front();

  auto anc = Sym->ancestors(I);
  SubspaceDimensions rm{I, anc, diag, Sym}; // equal to substruct[I]
  auto h = Zero_matrix<double>(rm.total());
}


TEST(Core, H2) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get(); // get the raw pointer
  Stats<double> stats(P, Sym->get_td_fields(), 0.0);
  Step step{P, RUNTYPE::NRG};
  Store<double> store(0,1);

  auto diag = setup_diag(P, Sym);
  auto n = new_subspaces(diag, Sym);
  std::cout << n << std::endl;
  SubspaceStructure substruct{diag, Sym};
  auto tasks = substruct.task_list();
  auto I = tasks.front();

//  auto h = hamiltonian<double>(step, I, opch, coef, diag, output, Sym, P);
//  dump_matrix(h);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
