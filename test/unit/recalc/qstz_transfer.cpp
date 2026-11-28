#include <gtest/gtest.h>

#include <bit>
#include <complex>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <core.hpp>

namespace {
using namespace NRG;

template<typename S>
void check_seed(S tiny) {
  initInvar({{"Q", additive}, {"SS", additive}, {"TZ", additive}});
  const Invar bra(1, 2, 0);
  MatrixElements<S> seed;
  EXPECT_EQ(qstz_transfer::infer("empty", seed), 0u);
  for (int tz : {-1, 0, 1, 2}) seed[{bra, Invar(0, 1, tz)}] = Matrix_traits<S>::Zero(2, 3);
  for (unsigned mask = 0; mask < 8; ++mask) {
    for (int d : {-1, 0, 1}) seed.at({bra, Invar(0, 1, -d)})(1, 2) = mask & (1u << (d+1)) ? tiny : S{-0.0};
    EXPECT_EQ(qstz_transfer::infer("mixed", seed), mask);
  }
  seed.at({bra, Invar(0, 1, 2)})(0, 0) = tiny;
  // The first three transfers already fill the mask; the rest must still be validated.
  EXPECT_THROW(qstz_transfer::infer("bad", seed), std::invalid_argument);
}

TEST(QstzTransfer, ExactReal) { check_seed(std::numeric_limits<double>::denorm_min()); }
TEST(QstzTransfer, ExactImaginary) { check_seed(std::complex<double>(0, std::numeric_limits<double>::denorm_min())); }

TEST(QstzTransfer, SparseInputAndStructuralFallbacks) {
  initInvar({{"Q", additive}, {"SS", additive}, {"TZ", additive}});
  Params P;
  const Invar bra(1, 2, 0);
  DiagInfo<double> diag;
  for (const auto &I : {bra, Invar(0, 1, -1), Invar(0, 1, 0), Invar(0, 1, 1)})
    diag[I] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  std::istringstream data("2\n1 2 0  0 1 -1  2.5\n1 2 0  0 1 0  0\n");
  CustomOp<double> named;
  named["Aep"] = MatrixElements<double>(data, diag);
  // Mask of Aep (7: all transfers, no restriction)
  const auto mask = [&](const bool enabled) {
    P.doublet_transfer_mask = enabled;
    const auto masks = doublet_transfer::seed_masks<double>(named, P, [](const auto &name, const auto &seed) -> std::optional<unsigned> {
      return qstz_transfer::infer(name, seed);
    });
    EXPECT_FALSE(masks.contains("unknown"));
    const auto it = masks.find("Aep");
    return it == masks.end() ? 7u : it->second;
  };
  EXPECT_EQ(mask(false), 7u);
  EXPECT_EQ(mask(true), 4u);
  P.h5ops = true;
  EXPECT_EQ(mask(true), 7u);
  P.h5ops = false;
  P.checksumrules = true;
  EXPECT_EQ(mask(true), 7u);
  P.checksumrules = false;
  for (const char *log : {"f", "r", "R"}) {
    P.logstr = std::string(log);
    EXPECT_EQ(mask(true), 7u) << log;
  }
  P.logstr = "";
  for (auto *request : {&P.specgt, &P.speci1t, &P.speci2t}) {
    *request = "Aep-Aep";
    EXPECT_EQ(mask(true), 7u);
    *request = "";
  }
  const auto &seed = named.at("Aep");
  ASSERT_EQ(seed.size(), 2u);
  EXPECT_EQ(seed.at({bra, Invar(0, 1, -1)})(0, 0), 2.5);
  EXPECT_EQ(seed.at({bra, Invar(0, 1, 0)})(0, 0), 0.0);
  EXPECT_FALSE(seed.contains({bra, Invar(0, 1, 1)}));
}

#ifdef NRG_SYM_ALL
TEST(QstzTransfer, NamedRecalcBothSpinSigns) {
  Params P;
  const auto sym = set_symmetry<double>(P, "QSTZ", 3);
  const Invar bra(1, 4, 0);
  const Twoinvar accidental_zero{bra, Invar(0, 5, -1)};
  std::vector<Invar> targets{bra};
  for (int tz : {-1, 0, 1}) for (int ss : {3, 5}) targets.emplace_back(0, ss, tz);

  // For both spin branches row (1,1) in the table copies the first ancestor's block with unit weight.
  DiagInfo<double> previous, diag;
  SubspaceStructure sub;
  for (const auto &I : targets) previous[sym->ancestor(I, 0)] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  for (const auto &I : {accidental_zero.first, accidental_zero.second}) previous[I] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  for (const auto &I : targets) {
    sub[I] = SubspaceDimensions(I, sym->ancestors(I), previous, sym.get());
    ASSERT_EQ(sub.at(I).rmax(0), 1u);
    diag[I] = NRG::Eigen<double>(std::vector<double>(sub.at(I).total(), 0.0), 1.0);
    split_in_blocks_Eigen(diag.at(I), sub.at(I)); // identity eigenvectors
  }

  P.finite = true;
  P.Nmax = 1;
  P.specd = "Aep-Probe";
  const Step step(P);
  MemTime mt;
  const auto check = [&](unsigned support, const char *mode, bool h5ops, unsigned effective) {
    SCOPED_TRACE(::testing::Message() << "support=" << support << " mode=" << (mode ? mode : "unset") << " h5ops=" << h5ops);
    P.h5ops = h5ops;
    Operators<double> ops;
    auto &op = ops.opd["Aep"];
    for (const auto &ket : targets) if (ket != bra) {
      const int d = -ket.get("TZ");
      const double value = support & (1u << (d+1)) ? 10.0*d + ket.get("SS") : 0.0;
      op[{sym->ancestor(bra, 0), sym->ancestor(ket, 0)}] = Matrix_traits<double>::Constant(1, 1, value);
    }
    op[accidental_zero] = Matrix_traits<double>::Zero(1, 1);
    auto &probe = ops.opd["Probe"];
    probe[{sym->ancestor(bra, 0), sym->ancestor(Invar(0, 3, -1), 0)}] = Matrix_traits<double>::Constant(1, 1, 9.0);
    ops.opd["unused"] = op;
    P.doublet_transfer_mask = mode && std::string_view(mode) == "1";
    ops.doublet_masks = doublet_transfer::seed_masks<double>(ops.opd, P, [&](const auto &name, const auto &seed) {
      return sym->doublet_transfer_mask(name, seed);
    });
    Oprecalc<double> recalc(RUNTYPE::NRG, ops, sym, mt, P);
    ASSERT_EQ(recalc.sl.size(), 1u);
    const auto &initial = recalc.sl.front().op1; // spectrum before the first recalc (calc0)
    ASSERT_EQ(initial.size(), 7u);
    ASSERT_EQ(initial.at(accidental_zero)(0, 0), 0.0);

    const auto dense = sym->recalc_doublet(diag, sub, op);
    ASSERT_EQ(dense.size(), 6u);
    recalc.recalculate_operators(ops, step, diag, sub, P);
    EXPECT_TRUE(ops.opd.at("unused").empty());
    ASSERT_EQ(op.size(), 2u * std::popcount(effective));
    for (const auto &ket : targets) if (ket != bra) {
      const Twoinvar key{bra, ket};
      const int d = -ket.get("TZ");
      const unsigned bit = 1u << (d+1);
      const auto &reference = dense.at(key);
      ASSERT_EQ(op.contains(key), bool(effective & bit)) << ket;
      if (op.contains(key)) {
        EXPECT_TRUE(op.at(key) == reference) << ket;
      }
      if (support & bit) {
        EXPECT_EQ(reference(0, 0), 10.0*d + ket.get("SS")) << ket;
      } else {
        EXPECT_TRUE((reference.array() == 0.0).all()) << ket;
      }
    }
    // This same sector pair had an explicit zero in the seed, but another ancestor feeds it now.
    if (support & 4u) {
      EXPECT_EQ(op.at(accidental_zero)(0, 0), 15.0);
    }
    // Probe has only +1 support. Both spin branches survive, even the one that is still zero.
    EXPECT_EQ(probe.size(), mode && std::string_view(mode) == "1" && !h5ops ? 2u : 6u);
    EXPECT_EQ(probe.at({bra, Invar(0, 3, -1)})(0, 0), 9.0);
    EXPECT_TRUE((probe.at(accidental_zero).array() == 0.0).all());
    if (h5ops) {
      H5Easy::File file(P.workdir->get() + "/ops.h5", H5Easy::File::Overwrite);
      ops.h5save(file, "1");
      for (const auto &[key, block] : dense) {
        const auto path = "1/d/Aep/" + key.first.name() + "/" + key.second.name() + "/matrix";
        ASSERT_TRUE(file.exist(path));
        EXPECT_TRUE(H5Easy::load<Matrix_traits<double>>(file, path) == block);
      }
    }
    if (effective == 0) {
      // The same iteration must still generate Wilson matrices for every transfer and spin branch.
      recalc_irreducible(step, diag, sub, ops.opch, sym.get(), mt, P);
      const auto &wilson = ops.opch.at(0).at(0);
      ASSERT_EQ(wilson.size(), 6u);
      for (const auto &[key, block] : dense) EXPECT_TRUE((wilson.at(key).array() != 0.0).any()) << key.second;
    }
  };
  for (unsigned mask = 0; mask < 8; ++mask) check(mask, "1", false, mask);
  check(4u, nullptr, false, 7u);
  check(4u, "0", false, 7u);
  check(4u, "1", true, 7u);
}
#endif

} // namespace
