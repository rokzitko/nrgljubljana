#include <gtest/gtest.h>

#include <initializer_list>
#include <vector>

#include <core.hpp>
#include <recalc.hpp>
#include "test_common.hpp"

using namespace NRG;

#ifdef NRG_SYM_ALL
namespace {

DiagInfo<double> one_state_diag(const std::initializer_list<Invar> sectors) {
  DiagInfo<double> diag;
  for (const auto &I : sectors)
    diag[I] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  return diag;
}

SubspaceStructure empty_ancestor_subspaces(const DiagInfo<double> &diag, const Symmetry<double> *Sym) {
  const DiagInfo<double> empty_previous_shell;
  SubspaceStructure substruct;
  for (const auto &I : diag.subspaces())
    substruct[I] = SubspaceDimensions{I, Sym->ancestors(I), empty_previous_shell, Sym};
  return substruct;
}

void expect_block(const MatrixElements<double> &blocks, const Invar &I1, const Invar &Ip) {
  EXPECT_TRUE(blocks.contains(Twoinvar(I1, Ip))) << "Missing block " << I1 << " -> " << Ip;
}

} // namespace
#endif

TEST(recalc, split_in_blocks_Eigen) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag(P, Sym);
#ifdef DISABLED_DUE_TO_ISSUES
  SubspaceStructure substruct{diag, Sym};
  substruct.dump();
#endif
}

#ifdef NRG_SYM_ALL
TEST(recalc, SPSU2LR_operator_ranks) { // NOLINT
  Params P;
  const auto Sym = set_symmetry<double>(P, "SPSU2LR", 2);
  const Invar center(5, 1);
  const auto diag = one_state_diag({center, Invar(4, 1), Invar(6, 1), Invar(3, 1), Invar(7, 1)});
  const auto substruct = empty_ancestor_subspaces(diag, Sym.get());
  const MatrixElements<double> old;

  const auto singlet = Sym->recalc_singlet(diag, substruct, old, 1);
  expect_block(singlet, center, center);

  const auto doublet = Sym->recalc_doublet(diag, substruct, old);
  expect_block(doublet, center, Invar(4, 1));
  expect_block(doublet, center, Invar(6, 1));

  const auto triplet = Sym->recalc_triplet(diag, substruct, old);
  expect_block(triplet, center, center);
  expect_block(triplet, center, Invar(3, 1));
  expect_block(triplet, center, Invar(7, 1));
}

TEST(recalc, SPU1LR_doublet_preserves_parity) { // NOLINT
  Params P;
  const auto Sym = set_symmetry<double>(P, "SPU1LR", 2);
  const Invar center(5, 1);
  const auto diag = one_state_diag({center, Invar(4, 1), Invar(6, 1)});
  const auto substruct = empty_ancestor_subspaces(diag, Sym.get());

  const auto doublet = Sym->recalc_doublet(diag, substruct, {});
  expect_block(doublet, center, Invar(4, 1));
  expect_block(doublet, center, Invar(6, 1));
}

TEST(recalc, SPSU2T_doublet_rank) { // NOLINT
  Params P;
  const auto Sym = set_symmetry<double>(P, "SPSU2T", 3);
  const Invar center(5, 3);
  const auto diag = one_state_diag({center, Invar(4, 2), Invar(6, 2), Invar(4, 3), Invar(6, 3), Invar(4, 4), Invar(6, 4)});
  const auto substruct = empty_ancestor_subspaces(diag, Sym.get());

  const auto doublet = Sym->recalc_doublet(diag, substruct, {});
  for (const int ss : {4, 6})
    for (const int t : {2, 3, 4})
      expect_block(doublet, center, Invar(ss, t));
}

TEST(recalc, DBLSU2_doublet_ranks) { // NOLINT
  Params P;
  const auto Sym = set_symmetry<double>(P, "DBLSU2", 2);
  const Invar center(5, 5);
  const auto diag = one_state_diag({center, Invar(4, 5), Invar(6, 5), Invar(5, 4), Invar(5, 6)});
  const auto substruct = empty_ancestor_subspaces(diag, Sym.get());

  const auto doublet = Sym->recalc_doublet(diag, substruct, {});
  expect_block(doublet, center, Invar(4, 5));
  expect_block(doublet, center, Invar(6, 5));
  expect_block(doublet, center, Invar(5, 4));
  expect_block(doublet, center, Invar(5, 6));
}

TEST(recalc, QSZTZ_triplet_components) { // NOLINT
  Params P;
  const auto Sym = set_symmetry<double>(P, "QSZTZ", 3);
  const Invar center(0, 5, 0);
  const auto diag = one_state_diag({center, Invar(0, 3, 0), Invar(0, 7, 0)});
  const auto substruct = empty_ancestor_subspaces(diag, Sym.get());

  const auto triplet = Sym->recalc_triplet(diag, substruct, {});
  expect_block(triplet, center, center);
  expect_block(triplet, center, Invar(0, 3, 0));
  expect_block(triplet, center, Invar(0, 7, 0));
}

TEST(recalc, corrected_spectral_factors) { // NOLINT
  Params P;
  auto Sym = set_symmetry<double>(P, "SPSU2LR", 2);
  EXPECT_DOUBLE_EQ(Sym->dynamicsusceptibility_factor(Invar(3, 1), Invar(5, 1)), 5.0 / 3.0);
  EXPECT_DOUBLE_EQ(Sym->specdens_factor(Invar(3, 1), Invar(4, 1)), 2.0);

  Sym = set_symmetry<double>(P, "DBLSU2", 2);
  EXPECT_DOUBLE_EQ(Sym->specdens_factor(Invar(3, 5), Invar(4, 5)), 10.0);
  EXPECT_DOUBLE_EQ(Sym->specdens_factor(Invar(3, 5), Invar(3, 6)), 9.0);

  Sym = set_symmetry<double>(P, "SPSU2T", 3);
  EXPECT_DOUBLE_EQ(Sym->specdens_factor(Invar(3, 0), Invar(4, 0)), 0.0);
}
#endif

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
