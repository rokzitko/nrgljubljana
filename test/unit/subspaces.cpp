#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <invar.hpp>
#include <subspaces.hpp>
#include <core.hpp>

using namespace NRG;

TEST(Core, H) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  EXPECT_EQ(P.combs, 4);
  EXPECT_EQ(SymSP->input_subspaces().size(), 4); // In
  EXPECT_EQ(SymSP->nr_combs(), 4);
  auto Sym = SymSP.get(); // get the raw pointer

  Stats<double> stats(P);
  stats.td.allfields.add(Sym->get_td_fields(), 1);
  stats.total_energy = 0.0;
  Step step{P, RUNTYPE::NRG};
  Store<double> store(0,1);

  auto diag = setup_diag(P, Sym);
  auto n = new_subspaces(diag, Sym);
  std::cout << n << std::endl;
  SubspaceStructure substruct{diag, Sym};
  auto tasks = substruct.task_list();
  auto I = tasks.front();
}

TEST(Subspaces, SubspaceStructure) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag(P, Sym);
  SubspaceStructure substruct{diag, Sym};
  auto tasks = substruct.task_list();
  std::cout << tasks << std::endl;
  substruct.dump();
  EXPECT_EQ(substruct.at_or_null(Invar(0,2)).total(), 5);
  EXPECT_EQ(substruct.at_or_null(Invar(0,0)).total(), 0);
}

TEST(Subspaces, SubspaceDimensions) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag3(P, Sym);
  Invar I(2,1);
  InvarVec ancestors = { Invar(0,1), Invar(1,2), Invar(2,1) }; // dims 2,3,4
  SubspaceDimensions sd(I, ancestors, diag, Sym, true); // true = ignore triangle inequality
  sd.dump();
  EXPECT_EQ(sd.combs(), 3);
  EXPECT_EQ(sd.rmax(0), 2);
  EXPECT_EQ(sd.rmax(1), 3);
  EXPECT_EQ(sd.rmax(2), 4);
  EXPECT_EQ(sd[0], 2);
  EXPECT_EQ(sd[1], 3);
  EXPECT_EQ(sd[2], 4);
  EXPECT_EQ(sd.exists(0), true);
  EXPECT_EQ(sd.exists(1), true);
  EXPECT_EQ(sd.exists(2), true);
  EXPECT_EQ(sd.offset(0), 0);
  EXPECT_EQ(sd.offset(1), 2);
  EXPECT_EQ(sd.offset(2), 5);
  EXPECT_EQ(sd.chunk(1), std::make_pair(0ul,2ul)); // first=first element, second=length
  EXPECT_EQ(sd.chunk(2), std::make_pair(2ul,3ul));
  EXPECT_EQ(sd.chunk(3), std::make_pair(5ul,4ul));
  EXPECT_EQ(sd.total(), 9);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
