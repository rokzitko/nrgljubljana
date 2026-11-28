#include <gtest/gtest.h>

#include <array>
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

// Transfers are bra minus ket; do not derive them from the recalc table filenames.
struct Transfer {
  int dsz, dtz;
  Invar ket(const Invar &bra) const { return Invar(bra.get("Q") - 1, bra.get("SZ") - dsz, bra.get("TZ") - dtz); }
  unsigned bit() const { return qsztz_transfer::bit(dsz, dtz); }
};
constexpr std::array<Transfer, 6> branches{{{-1, -1}, {-1, 0}, {-1, 1}, {1, -1}, {1, 0}, {1, 1}}};
constexpr unsigned all = 63;

TEST(QsztzTransfer, ExactJointSupport) {
  initInvar({{"Q", additive}, {"SZ", additive}, {"TZ", additive}});
  using S = std::complex<double>;
  const S tiny(0.0, std::numeric_limits<double>::denorm_min());
  const Invar bra(1, 3, -2);
  MatrixElements<S> seed;
  EXPECT_EQ(qsztz_transfer::infer("empty", seed), 0u);
  unsigned used = 0;
  for (const auto &b : branches) {
    ASSERT_TRUE(std::has_single_bit(b.bit()));
    ASSERT_EQ(used & b.bit(), 0u);
    used |= b.bit();
    seed[{bra, b.ket(bra)}] = Matrix_traits<S>::Zero(2, 3);
  }
  ASSERT_EQ(used, all);
  for (unsigned mask = 0; mask <= all; ++mask) {
    for (const auto &b : branches) seed.at({bra, b.ket(bra)})(1, 2) = mask & b.bit() ? tiny : S(-0.0, -0.0);
    EXPECT_EQ(qsztz_transfer::infer("mixed", seed), mask) << mask;
  }
  // An invalid nonzero rectangle *after* all six valid bits must not be ignored.
  const Invar last(2, 1, -1);
  for (const auto &bad : {Transfer{0, 0}, Transfer{2, -1}, Transfer{1, 2}, Transfer{-1, -2}}) {
    const Twoinvar key{last, bad.ket(last)};
    auto &block = seed[key] = Matrix_traits<S>::Zero(2, 3);
    EXPECT_EQ(qsztz_transfer::infer("mixed", seed), all);
    block(1, 2) = tiny;
    EXPECT_THROW(qsztz_transfer::infer("mixed", seed), std::invalid_argument) << bad.dsz << ',' << bad.dtz;
    seed.erase(key);
  }
}

TEST(QsztzTransfer, SparseInputAndDenseFallbacks) {
  initInvar({{"Q", additive}, {"SZ", additive}, {"TZ", additive}});
  Params P;
  const Invar bra(1, 3, -2), active(0, 2, -1), zero(0, 4, -3), absent(0, 2, -2);
  DiagInfo<double> diag;
  diag[bra] = NRG::Eigen<double>(std::vector<double>(2), 1.0);
  for (const auto &I : {active, zero, absent}) diag[I] = NRG::Eigen<double>(std::vector<double>(3), 1.0);
  std::istringstream data("2\n1 3 -2  0 2 -1  0 0 0 0 0 2.5\n1 3 -2  0 4 -3  0 0 0 0 0 -0\n");
  Operators<double> ops;
  auto &op = ops.opd["mixed"] = MatrixElements<double>(data, diag);
  // Mask of 'mixed' (all: no restriction)
  const auto mask = [&](const bool enabled) {
    P.doublet_transfer_mask = enabled;
    const auto masks = doublet_transfer::seed_masks<double>(ops.opd, P, [](const auto &name, const auto &seed) -> std::optional<unsigned> {
      return qsztz_transfer::infer(name, seed);
    });
    EXPECT_FALSE(masks.contains("unknown"));
    const auto it = masks.find("mixed");
    return it == masks.end() ? all : it->second;
  };
  EXPECT_EQ(mask(false), all);
  EXPECT_EQ(mask(true), qsztz_transfer::bit(+1, -1));
  P.h5ops = true;
  EXPECT_EQ(mask(true), all);
  P.h5ops = false;
  P.checksumrules = true;
  EXPECT_EQ(mask(true), all);
  P.checksumrules = false;
  for (const char *log : {"f", "r", "R"}) {
    P.logstr = std::string(log);
    EXPECT_EQ(mask(true), all) << log;
  }
  P.logstr = "";
  for (auto *request : {&P.specgt, &P.speci1t, &P.speci2t}) {
    *request = "mixed-mixed";
    EXPECT_EQ(mask(true), all);
    *request = "";
  }
  ASSERT_EQ(op.size(), 2u);
  EXPECT_EQ(op.at({bra, active})(1, 2), 2.5);
  EXPECT_TRUE((op.at({bra, zero}).array() == 0.0).all());
  EXPECT_FALSE(op.contains({bra, absent}));

  // Opt-in validates nonzero charge transfer even when a consumer requires dense blocks.
  op[{bra, Invar(1, 2, -1)}] = Matrix_traits<double>::Ones(2, 3);
  P.h5ops = true;
  EXPECT_EQ(mask(false), all);
  EXPECT_THROW(mask(true), std::invalid_argument);
}

#ifdef NRG_SYM_ALL
TEST(QsztzTransfer, NamedDispatchAllJointMasks) {
  Params P;
  const auto sym = set_symmetry<double>(P, "QSZTZ", 3);
  const Invar bra(1, 3, -2);
  const Twoinvar accidental_zero{bra, branches.front().ket(bra)};
  std::vector<Invar> targets{bra};
  for (const auto &b : branches) targets.push_back(b.ket(bra));

  // First local-state table row copies one ancestor block with unit coefficient.
  DiagInfo<double> previous, diag;
  SubspaceStructure sub;
  for (const auto &I : targets) previous[sym->ancestor(I, 0)] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  for (const auto &I : {accidental_zero.first, accidental_zero.second}) previous[I] = NRG::Eigen<double>(std::vector<double>{0.0}, 1.0);
  for (const auto &I : targets) {
    sub[I] = SubspaceDimensions(I, sym->ancestors(I), previous, sym.get());
    ASSERT_EQ(sub.at(I).rmax(0), 1u);
    diag[I] = NRG::Eigen<double>(std::vector<double>(sub.at(I).total()), 1.0);
    split_in_blocks_Eigen(diag.at(I), sub.at(I)); // identity eigenvectors
  }

  P.finite = true;
  P.Nmax = 1;
  P.specd = "A-Probe";
  const Step step(P);
  MemTime mt;
  const auto check = [&](unsigned support, const char *flag, bool sparse) {
    SCOPED_TRACE(::testing::Message() << "support=" << support << " flag=" << (flag ? flag : "unset") << " h5ops=" << P.h5ops);
    Operators<double> ops;
    auto &op = ops.opd["A"], &probe = ops.opd["Probe"];
    for (const auto &b : branches) {
      const Twoinvar key{sym->ancestor(bra, 0), sym->ancestor(b.ket(bra), 0)};
      op[key] = Matrix_traits<double>::Constant(1, 1, support & b.bit() ? b.bit() : 0u);
      if (b.bit() == branches.back().bit()) probe[key] = Matrix_traits<double>::Constant(1, 1, -3.0);
    }
    op[accidental_zero] = Matrix_traits<double>::Zero(1, 1);
    P.doublet_transfer_mask = flag && std::string_view(flag) == "1";
    ops.doublet_masks = doublet_transfer::seed_masks<double>(ops.opd, P, [&](const auto &name, const auto &seed) {
      return sym->doublet_transfer_mask(name, seed);
    });
    Oprecalc<double> recalc(RUNTYPE::NRG, ops, sym, mt, P);
    ASSERT_EQ(recalc.sl.size(), 1u);
    EXPECT_EQ(recalc.sl.front().op1.at(accidental_zero)(0, 0), 0.0);

    const auto dense = sym->recalc_doublet(diag, sub, op);
    ASSERT_EQ(dense.size(), branches.size());
    recalc.recalculate_operators(ops, step, diag, sub, P);
    const unsigned effective = sparse ? support : all;
    ASSERT_EQ(op.size(), static_cast<size_t>(std::popcount(effective)));
    ASSERT_EQ(probe.size(), sparse ? 1u : branches.size());
    for (const auto &b : branches) {
      const Twoinvar key{bra, b.ket(bra)};
      SCOPED_TRACE(::testing::Message() << "dsz=" << b.dsz << " dtz=" << b.dtz);
      ASSERT_TRUE(dense.contains(key));
      const auto &reference = dense.at(key);
      EXPECT_EQ(reference(0, 0), support & b.bit() ? b.bit() : 0u);
      if (!(support & b.bit())) {
        EXPECT_TRUE((reference.array() == 0.0).all());
      }
      ASSERT_EQ(op.contains(key), bool(effective & b.bit()));
      if (op.contains(key)) {
        EXPECT_TRUE(op.at(key) == reference);
      }
      ASSERT_EQ(probe.contains(key), !sparse || b.bit() == branches.back().bit());
    }
    EXPECT_EQ(probe.at({bra, branches.back().ket(bra)})(0, 0), -3.0);
    // This exact seed rectangle was zero; its allowed transfer can populate it next shell.
    if (support & branches.front().bit()) {
      EXPECT_EQ(op.at(accidental_zero)(0, 0), branches.front().bit());
    }
    if (effective == 0u) {
      recalc_irreducible(step, diag, sub, ops.opch, sym.get(), mt, P);
      const auto &wilson = ops.opch.at(0).at(0);
      ASSERT_EQ(wilson.size(), branches.size());
      for (const auto &b : branches) EXPECT_TRUE((wilson.at({bra, b.ket(bra)}).array() != 0.0).any());
    }
  };
  for (unsigned mask = 0; mask <= all; ++mask) check(mask, "1", true);
  const unsigned mixed = qsztz_transfer::bit(-1, +1) | qsztz_transfer::bit(+1, -1);
  check(mixed, nullptr, false);
  check(mixed, "0", false);
  P.h5ops = true;
  check(mixed, "1", false);
}
#endif
} // namespace
