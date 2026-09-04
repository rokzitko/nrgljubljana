#include <string>
#include <sstream>
#include <vector>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <core.hpp>
#include <floquet.hpp>
#include <truncation.hpp>

using namespace NRG;

TEST(Core, HighestRetainedUsesSafeguardWithoutFullSort) { // NOLINT
  Params P;
  P.keep = 2;
  P.keepmin = 1;
  P.safeguard = 0.1;
  P.safeguardmax = 2;

  DiagInfo<double> diag;
  diag[Invar()] = NRG::Eigen<double>(std::vector<double>{2.0, 0.0, 1.0, 1.05}, 1.0);

  Step step(P);
  EXPECT_DOUBLE_EQ(highest_retained(step, diag, P), 1.05);
}

TEST(Core, FloquetCriterionPreservesExtendedZoneOrdering) { // NOLINT
  const Invar sector;
  DiagInfo<double> diag;
  diag[sector] = NRG::Eigen<double>(std::vector<double>{10.0, 4.0, 7.0}, 1.0);
  MatrixElements<double> mode;
  mode[Twoinvar(sector, sector)] = zero_matrix<double>(3);
  mode.at(Twoinvar(sector, sector))(0, 0) = 2.0;
  mode.at(Twoinvar(sector, sector))(1, 1) = -1.0;
  mode.at(Twoinvar(sector, sector))(2, 2) = 0.5;

  const auto prepared = calculate_floquet_criteria(diag, mode, 2.0, {-1.0, 2.0});
  EXPECT_DOUBLE_EQ(prepared.minimum_energy, 4.0);
  EXPECT_EQ(prepared.values.at(sector), (std::vector<double>{10.0, 8.0, 7.0}));

  for (size_t i = 0; i < prepared.values.at(sector).size(); i++)
    diag.at(sector).values.set_crit(i, prepared.values.at(sector)[i]);
  diag.sort_by_c();
  const auto lowest_criterion = diag.find_Clw();
  diag.shift(prepared.minimum_energy, lowest_criterion);
  EXPECT_EQ(diag.at(sector).values.all_rel(), (std::vector<double>{3.0, 0.0, 6.0}));
  EXPECT_EQ(diag.at(sector).values.all_crit(), (std::vector<double>{0.0, 1.0, 3.0}));
}

TEST(Core, FloquetCriterionRejectsMissingBlockWithoutInsertion) { // NOLINT
  DiagInfo<double> diag;
  diag[Invar()] = NRG::Eigen<double>(std::vector<double>{1.0}, 1.0);
  MatrixElements<double> mode;
  const auto original_size = mode.size();

  EXPECT_THROW(calculate_floquet_criteria(diag, mode, 1.0, {0.0, 0.0}), std::invalid_argument);
  EXPECT_EQ(mode.size(), original_size);
}

TEST(Core, FloquetCriterionValidatesNumbersAndShapeBeforeMutation) { // NOLINT
  const Invar sector;
  DiagInfo<std::complex<double>> diag;
  diag[sector] = NRG::Eigen<std::complex<double>>(std::vector<double>{1.0, 2.0}, 1.0);
  MatrixElements<std::complex<double>> mode;
  mode[Twoinvar(sector, sector)] = zero_matrix<std::complex<double>>(2);
  mode.at(Twoinvar(sector, sector))(0, 0) = -1.0;
  mode.at(Twoinvar(sector, sector))(1, 1) = 1.0;
  const auto original_criteria = diag.at(sector).values.all_crit();

  auto wrong_shape = mode;
  wrong_shape.at(Twoinvar(sector, sector)).resize(1, 1);
  EXPECT_THROW(calculate_floquet_criteria(diag, wrong_shape, 1.0, {-1.0, 1.0}), std::invalid_argument);

  auto complex_diagonal = mode;
  complex_diagonal.at(Twoinvar(sector, sector))(0, 0) = std::complex<double>(-1.0, 1e-3);
  EXPECT_THROW(calculate_floquet_criteria(diag, complex_diagonal, 1.0, {-1.0, 1.0}), std::invalid_argument);

  auto out_of_range = mode;
  out_of_range.at(Twoinvar(sector, sector))(1, 1) = 2.0;
  EXPECT_THROW(calculate_floquet_criteria(diag, out_of_range, 1.0, {-1.0, 1.0}), std::invalid_argument);

  auto nonfinite = mode;
  nonfinite.at(Twoinvar(sector, sector))(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(calculate_floquet_criteria(diag, nonfinite, 1.0, {-1.0, 1.0}), std::invalid_argument);

  EXPECT_EQ(diag.at(sector).values.all_crit(), original_criteria);
}

TEST(Core, FloquetCriterionAcceptsRoundoffScaleImaginaryPart) { // NOLINT
  const Invar sector;
  DiagInfo<std::complex<double>> diag;
  diag[sector] = NRG::Eigen<std::complex<double>>(std::vector<double>{1.0}, 1.0);
  MatrixElements<std::complex<double>> mode;
  mode[Twoinvar(sector, sector)] = zero_matrix<std::complex<double>>(1);
  mode.at(Twoinvar(sector, sector))(0, 0) = std::complex<double>(0.5, 1e-15);

  EXPECT_NO_THROW(calculate_floquet_criteria(diag, mode, 1.0, {-1.0, 1.0}));
}

TEST(Core, FloquetModeOperatorIsRetainedForDmSweep) { // NOLINT
  Params P;
  P.floquet = true;
  Step dm_step(P, RUNTYPE::DMNRG);
  Oprecalc<double>::Ops selected;

  EXPECT_TRUE(selected.do_s("m", P, dm_step));
  EXPECT_FALSE(selected.do_s("other", P, dm_step));
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
