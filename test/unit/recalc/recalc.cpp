#include <gtest/gtest.h>
#include <gsl/gsl_sf_coupling.h>

#include <cmath>
#include <complex>
#include <initializer_list>
#include <filesystem>
#include <numeric>
#include <string>
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

std::vector<Recalc<double>> triplet_table(const bool qstz, const int ss1, const int branch) {
  const int q1 = 0, tz1 = 0;
  if (qstz) {
    if (branch == 0) return {
#include "qstz/qstz-triplets.dat"
    };
    if (branch == 1) return {
#include "qstz/qstz-tripletp.dat"
    };
    return {
#include "qstz/qstz-tripletm.dat"
    };
  }
  if (branch == 0) return {
#include "qs/qs-3ch-triplets.dat"
  };
  if (branch == 1) return {
#include "qs/qs-3ch-tripletp.dat"
  };
  return {
#include "qs/qs-3ch-tripletm.dat"
  };
}

bool spin_triangle(const int a, const int b, const int c) {
  return a >= 0 && b >= 0 && c >= 0 && std::abs(a-b) <= c && c <= a+b && (a+b+c) % 2 == 0;
}

// All arguments are twice the angular momentum; contract the q=+1 component at M=S.
double triplet_cg_sum(const int s, const int sp, const int j, const int jp, const int l) {
  const auto cg = [](int a, int ma, int b, int mb, int c, int mc) {
    if (!spin_triangle(a, b, c) || ma+mb != mc || std::abs(ma) > a || std::abs(mb) > b || std::abs(mc) > c ||
        (a-ma) % 2 || (b-mb) % 2 || (c-mc) % 2) return 0.0;
    return ((a-b+mc)/2 % 2 ? -1.0 : 1.0) * std::sqrt(c+1.0) * gsl_sf_coupling_3j(a, b, c, ma, mb, -mc);
  };
  double sum = 0;
  for (int mu = -l; mu <= l; mu += 2)
    sum += cg(j, s-mu, l, mu, s, s) * cg(jp, s-2-mu, l, mu, sp, s-2) * cg(jp, s-2-mu, 2, 2, j, s-mu);
  return sum / cg(sp, s-2, 2, 2, s, s);
}

void check_triplet_coefficients(const bool qstz) {
  Params P;
  const auto Sym = set_symmetry<double>(P, qstz ? "QSTZ" : "QS", 3);
  for (int ss1 = 1; ss1 <= 21; ++ss1)
    for (const int branch : {0, 1, -1}) {
      const int s = ss1-1, sp = s+2*branch;
      if (!spin_triangle(s, sp, 2)) continue;
      const auto I1 = qstz ? Invar(0, ss1, 0) : Invar(0, ss1);
      const auto Ip = qstz ? Invar(0, ss1+2*branch, 0) : Invar(0, ss1+2*branch);
      size_t checked = 0;
      for (const auto &row : triplet_table(qstz, ss1, branch)) {
        SCOPED_TRACE(::testing::Message() << "ss1=" << ss1 << " branch=" << branch << " row=" << row.i1 << ',' << row.ip);
        ASSERT_GE(row.i1, 1u);
        ASSERT_LE(row.i1, Sym->nr_combs());
        ASSERT_GE(row.ip, 1u);
        ASSERT_LE(row.ip, Sym->nr_combs());
        EXPECT_EQ(row.IN1, Sym->ancestor(I1, row.i1-1));
        EXPECT_EQ(row.INp, Sym->ancestor(Ip, row.ip-1));
        const auto local = Sym->QN_subspace(row.i1-1);
        EXPECT_EQ(local, Sym->QN_subspace(row.ip-1));
        const int j = row.IN1.get("SS")-1, jp = row.INp.get("SS")-1, l = local.get("SS")-1;
        if (!spin_triangle(j, l, s) || !spin_triangle(jp, l, sp) || !spin_triangle(j, jp, 2)) continue;
        ++checked;
        [&] { // A fatal row assertion must not hide failures in subsequent rows.
          ASSERT_TRUE(std::isfinite(row.factor));
          EXPECT_NEAR(row.factor, triplet_cg_sum(s, sp, j, jp, l), 1e-12);
        }();
      }
      EXPECT_GT(checked, 0u);
    }
}

struct TripletEndpoint {
  const char *name;
  bool qstz;
  int ss1, branch;
  size_t i1, ip;
  double expected;
};
class TripletEndpoints : public ::testing::TestWithParam<TripletEndpoint> {};

template<typename Scalar>
void check_triplet_recalc(const TripletEndpoint &e, const Scalar amplitude) {
  Params P;
  const auto Sym = set_symmetry<Scalar>(P, e.qstz ? "QSTZ" : "QS", 3);
  const auto I1 = e.qstz ? Invar(0, e.ss1, 0) : Invar(0, e.ss1);
  const auto Ip = e.qstz ? Invar(0, e.ss1+2*e.branch, 0) : Invar(0, e.ss1+2*e.branch);
  const auto IN1 = Sym->ancestor(I1, e.i1-1), INp = Sym->ancestor(Ip, e.ip-1);
  const double factor = triplet_cg_sum(e.ss1-1, e.ss1-1+2*e.branch, IN1.get("SS")-1, INp.get("SS")-1,
                                       Sym->QN_subspace(e.i1-1).get("SS")-1);
  ASSERT_NEAR(factor, e.expected, 1e-12);
  DiagInfo<Scalar> previous, diag;
  for (const auto &I : {IN1, INp}) previous[I] = NRG::Eigen<Scalar>(std::vector<double>{0.0}, 1.0);
  SubspaceStructure substruct;
  for (const auto &I : {I1, Ip}) {
    auto &sub = substruct[I] = SubspaceDimensions{I, Sym->ancestors(I), previous, Sym.get()};
    ASSERT_GT(sub.total(), 0);
    diag[I] = NRG::Eigen<Scalar>(std::vector<double>(sub.total(), 0.0), 1.0); // Identity eigenvectors.
    split_in_blocks_Eigen(diag.at(I), sub);
  }
  ASSERT_EQ(substruct.at(I1).rmax(e.i1-1), 1u);
  ASSERT_EQ(substruct.at(Ip).rmax(e.ip-1), 1u);
  MatrixElements<Scalar> cold;
  cold[Twoinvar(IN1, INp)] = Matrix_traits<Scalar>::Constant(1, 1, amplitude);
  const auto result = Sym->recalc_triplet(diag, substruct, cold);
  ASSERT_TRUE(result.contains(Twoinvar(I1, Ip)));
  const auto &block = result.at(Twoinvar(I1, Ip));
  ASSERT_EQ(size1(block), substruct.at(I1).total());
  ASSERT_EQ(size2(block), substruct.at(Ip).total());
  ASSERT_TRUE(block.allFinite());
  EXPECT_NEAR(std::abs(block(substruct.at(I1).offset(e.i1-1), substruct.at(Ip).offset(e.ip-1)) - factor*amplitude), 0.0, 1e-12);
}

TEST(recalc, QSTZ_triplet_coefficients) { check_triplet_coefficients(true); }
TEST(recalc, QS_three_channel_triplet_coefficients) { check_triplet_coefficients(false); }
TEST_P(TripletEndpoints, real_recalc) { check_triplet_recalc(GetParam(), 2.0); }
TEST_P(TripletEndpoints, complex_recalc) { check_triplet_recalc(GetParam(), std::complex<double>(0.75, -1.25)); }
INSTANTIATE_TEST_SUITE_P(recalc, TripletEndpoints, ::testing::Values(
  TripletEndpoint{"QSTZ_s_16_15_ss2", true, 2, 0, 16, 15, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_19_18_ss2", true, 2, 0, 19, 18, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_22_21_ss2", true, 2, 0, 22, 21, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_51_50_ss2", true, 2, 0, 51, 50, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_54_53_ss2", true, 2, 0, 54, 53, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_57_56_ss2", true, 2, 0, 57, 56, -2*std::sqrt(2.0)/3},
  TripletEndpoint{"QSTZ_s_42_41_ss2", true, 2, 0, 42, 41, -std::sqrt(5.0/6)},
  TripletEndpoint{"QSTZ_p_42_40_ss1", true, 1, 1, 42, 40, 1.0},
  TripletEndpoint{"QS_s_42_41_ss2", false, 2, 0, 42, 41, -std::sqrt(5.0/6)},
  TripletEndpoint{"QSTZ_s_control", true, 3, 0, 16, 15, -std::sqrt(5.0/3)/2},
  TripletEndpoint{"QSTZ_p_control", true, 3, 1, 42, 40, 1/(2*std::sqrt(5.0))},
  TripletEndpoint{"QSTZ_m_control", true, 3, -1, 40, 42, 1/std::sqrt(6.0)},
  TripletEndpoint{"QS_s_control", false, 2, 0, 16, 15, -2*std::sqrt(2.0)/3}),
  [](const auto &test_case) { return test_case.param.name; });

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

TEST(recalc, h5save_blocks_skips_zero_width_ancestors) { // NOLINT
  Params P;
  const auto SymSP = setup_Sym<double>(P);
  auto *Sym = SymSP.get();
  const auto previous = setup_diag<double>(P, Sym);
  const SubspaceStructure substruct(previous, Sym);
  DiagInfo<double> diag;

  for (const auto &[I, sub] : substruct) {
    if (sub.total() == 0) continue;
    std::vector<double> values(sub.total());
    std::iota(values.begin(), values.end(), 0.0);
    diag[I] = NRG::Eigen<double>(values, 1.0);
    split_in_blocks_Eigen(diag.at(I), sub, true);
  }

  {
    auto h5 = H5Easy::File("recalc-blocks.h5", H5Easy::File::Overwrite);
    h5save_blocks(h5, "U/", diag, substruct);
    for (const auto &[I, sub] : substruct)
      for (const auto i : range0(sub.combs())) {
        const auto path = "U/" + I.name() + "/" + sub.ancestor(i).name() + "|nr=" + std::to_string(i);
        EXPECT_EQ(h5.exist(path), sub.exists(i));
        if (sub.exists(i)) {
          const auto block = H5Easy::load<EigenMatrix<double>>(h5, path);
          EXPECT_EQ(size1(block), diag.at(I).getnrstored());
          EXPECT_EQ(size2(block), sub.rmax(i));
        }
      }
  }
  std::filesystem::remove("recalc-blocks.h5");
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
