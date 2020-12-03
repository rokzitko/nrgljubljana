#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <invar.hpp>
#include <symmetry.hpp>

using namespace NRG;

TEST(Symmetry, Symmetry) { // NOLINT
  Params P;
  P.symtype.setvalue("QS");
  P.set_channels(1);
  auto Sym = set_symmetry<double>(P);
  EXPECT_EQ(P.combs, 4);
  EXPECT_EQ(Sym->nr_combs(), 4);
  const auto In = Sym->input_subspaces();
  // Differences wrt QNs in the previous shell to get *new* (Q,Sz)
  EXPECT_EQ(In[0], Invar(1,0));  // empty: (Q+1,Sz) -> (Q,Sz)
  EXPECT_EQ(In[1], Invar(0,-1)); // spin up: (Q,Sz-1/2) -> (Q,Sz)
  EXPECT_EQ(In[2], Invar(0,1));  // spin down: (Q,Sz+1/2) -> (Q,Sz)
  EXPECT_EQ(In[3], Invar(-1,0)); // double occupancy: (Q-1,Sz) -> (Q,Sz)
  // (Q,S) quantum numbers of of the new site
  EXPECT_EQ(Sym->QN_subspace(0), Invar(-1,1));
  EXPECT_EQ(Sym->QN_subspace(1), Invar(0,2));
  EXPECT_EQ(Sym->QN_subspace(2), Invar(0,2));
  EXPECT_EQ(Sym->QN_subspace(3), Invar(1,1));
  EXPECT_EQ(Sym->ancestor(Invar(0,2), 0), Invar(1,2));
  EXPECT_EQ(Sym->ancestor(Invar(0,2), 1), Invar(0,1));
  EXPECT_EQ(Sym->ancestor(Invar(0,2), 2), Invar(0,3));
  EXPECT_EQ(Sym->ancestor(Invar(0,2), 3), Invar(-1,2));
  auto anc = Sym->ancestors(Invar(0,2));
  EXPECT_EQ(anc[0], Invar(1,2));
  EXPECT_EQ(anc[1], Invar(0,1));
  EXPECT_EQ(anc[2], Invar(0,3));
  EXPECT_EQ(anc[3], Invar(-1,2));
  auto subs = Sym->new_subspaces(Invar(0,0));
  EXPECT_EQ(subs[0], Invar(-1,0));
  EXPECT_EQ(subs[1], Invar(0,1));
  EXPECT_EQ(subs[2], Invar(0,-1));
  EXPECT_EQ(subs[3], Invar(1,0));
  EXPECT_EQ(Sym->mult(Invar(1,2)), 2);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
