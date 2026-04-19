#include <vector>

#include <gtest/gtest.h>

#include "test_common.hpp"
#include "test_clean.hpp"
#include <store.hpp>
#include <subspaces.hpp>

using namespace NRG;

TEST(Store, StoredEigenPreservesValuesAndRanges) { // NOLINT
  Params P;
  RawEigen<double> raw(3, 3);
  raw.val = {1.0, 2.0, 4.0};
  raw.vec = id_matrix<double>(3);
  Step step(P);
  NRG::Eigen<double> eig(std::move(raw), step);
  eig.values.set_scale(2.0);
  eig.values.set_shift(1.0);
  eig.values.set_T_shift(3.0);
  eig.truncate_prepare(2);
  eig.truncate_perform();

  StoredEigen<double> stored(eig, false);

  EXPECT_EQ(stored.kept(), 2U);
  EXPECT_EQ(stored.total(), 3U);
  EXPECT_EQ(stored.min(), 2U);
  EXPECT_EQ(stored.max(), 3U);
  EXPECT_EQ(std::vector<size_t>(stored.all().begin(), stored.all().end()), (std::vector<size_t>{2U}));
  EXPECT_DOUBLE_EQ(stored.values.abs_zero(2), 6.0);
  EXPECT_DOUBLE_EQ(stored.values.abs_T(2), 9.0);

  stored.subtract_GS_energy(1.5);
  EXPECT_DOUBLE_EQ(stored.values.abs_G(2), 7.5);

  StoredEigen<double> stored_last(eig, true);
  EXPECT_EQ(std::vector<size_t>(stored_last.all().begin(), stored_last.all().end()), (std::vector<size_t>{0U, 1U, 2U}));
}

TEST(Store, ThermoAndBackiterBuildersKeepExpectedMetadata) { // NOLINT
  Params P;
  setup_P_clean(P);
  auto SymSP = setup_Sym<double>(P);
  auto *Sym = SymSP.get();
  auto diag = setup_diag_clean(P, Sym);
  SubspaceStructure substruct{diag, Sym};

  auto thermo = make_thermo_subs(diag, false);
  auto backiter = make_backiter_subs(diag, substruct, false);

  ASSERT_EQ(thermo.size(), diag.size());
  ASSERT_EQ(backiter.size(), diag.size());

  for (const auto &[I, eig] : diag) {
    ASSERT_TRUE(thermo.count(I));
    ASSERT_TRUE(backiter.count(I));
    EXPECT_EQ(thermo.at(I).kept(), eig.getnrkept());
    EXPECT_EQ(thermo.at(I).total(), eig.getdim());
    EXPECT_EQ(backiter.at(I).kept(), eig.getnrkept());
    EXPECT_EQ(backiter.at(I).total(), eig.getdim());
    EXPECT_EQ(backiter.at(I).rmax.total(), substruct.at_or_null(I).total());
    EXPECT_DOUBLE_EQ(thermo.at(I).eig.values.abs_zero(0), eig.values.abs_zero(0));
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
