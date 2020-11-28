#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <invar.hpp>
#include <subspaces.hpp>
#include <core.hpp>

using namespace NRG;

TEST(Core, H) { // NOLINT
  auto [P, SymSP] = test_setup_basic<double>();
  EXPECT_EQ(P.combs, 4);
  EXPECT_EQ(SymSP->input_subspaces().size(), 4); // In
//  EXPECT_EQ(SymSP->nr_combs(), 4);
/*  
  auto Sym = SymSP.get();
  EXPECT_EQ(Sym->input_subspaces().size(), 4); // In
  EXPECT_EQ(Sym->nr_combs(), 4);
  
  Stats<double> stats(P);
  stats.td.allfields.add(Sym->get_td_fields(), 1);
  stats.total_energy = 0.0;
  Step step{P, RUNTYPE::NRG};
  Store<double> store(0,1);
  
  std::set<Invar> subspaces;
  for (const auto &I : diag.subspaces()) {
    const auto all = Sym->new_subspaces(I);
    std::cout << I << " -> " << all << std::endl;
  }
*/
//  auto n = new_subspaces(diag, Sym);
//  std::cout << n << std::endl;
//  SubspaceStructure substruct{diag, Sym};
//  auto tasks = substruct.task_list();
//  auto I = tasks.front();
}
      

/*
TEST(Subspaces, SubspaceStructure) { // NOLINT
  auto [P, Sym, diag] = test_setup_diag();
  SubspaceStructure substruct{diag, Sym.get()};
  auto tasks = substruct.task_list();
  std::cout << tasks << std::endl;
}
*/

TEST(Subspaces, SubspaceDimensions) { // NOLINT
/*  auto [P, Sym] = test_setup_basic();
  Invar I(2,1);
  const auto anc = Sym->ancestors(I);
  SubspaceDimensions rm{I, anc, diagprev, Sym.get()};
 *./
  
  /*
  auto [P, Sym, diag] = test_setup_diag3();
  Invar I(2,1);
  InvarVec ancestors = { Invar(0,1), Invar(1,2), Invar(2,0) }; // dims 2,3,4
  SubspaceDimensions sd(I, ancestors, diag, Sym.get()); // use get()!
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
  EXPECT_EQ(sd.chunk(1), std::make_pair(0ul,2ul));
  EXPECT_EQ(sd.chunk(2), std::make_pair(2ul,5ul));
  EXPECT_EQ(sd.chunk(3), std::make_pair(5ul,9ul));
*/
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
