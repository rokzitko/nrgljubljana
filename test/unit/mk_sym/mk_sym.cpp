#include <gtest/gtest.h>
#include <complex>
#include <mk_sym.hpp>
#include <symmetry.hpp>
#include <params.hpp>
#include <outfield.hpp>
#include <symmetry_traits.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using namespace NRG;

TEST(mk_sym, QS) {
  Params P;
  auto sym = mk_QS<double>(P);
}

TEST(mk_sym, QSZ) {
  Params P;
  auto sym = mk_QSZ<double>(P);
}

TEST(mk_sym, get) {
  Params P;
  auto sym = get<double>("QS", P);
}

TEST(mk_sym, QS_complex) {
  Params P;
  auto sym = mk_QS<std::complex<double>>(P);
}

TEST(mk_sym, QSZ_complex) {
  Params P;
  auto sym = mk_QSZ<std::complex<double>>(P);
}

TEST(mk_sym, get_complex) {
  Params P;
  auto sym = get<std::complex<double>>("QS", P);
}

// symmetry_traits lists every symmetry type of get(), with its isfield(). Types not compiled in
// (SYM_MORE/SYM_ALL off, or C3 types for real S) throw from get() and are skipped.
template <typename S> int check_symmetry_traits() {
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    Params P;
    try {
      auto sym = get<S>(std::string(t.name), P);
      EXPECT_EQ(sym->isfield(), t.field) << t.name;
      ++checked;
    } catch (const std::runtime_error &) {}
  }
  return checked;
}

TEST(mk_sym, symmetry_traits) {
  EXPECT_GE(check_symmetry_traits<double>() + check_symmetry_traits<std::complex<double>>(), 4);
  EXPECT_EQ(find_symmetry_traits("XYZ"), nullptr);
  ASSERT_NE(find_symmetry_traits("QSZ"), nullptr);
  EXPECT_TRUE(find_symmetry_traits("QSZ")->field);
}

// load() fills In and QN, which have P.combs+1 entries. A channel count or a substeps setting that the coefficient
// tables of a symmetry type do not cover must be rejected before any entry is written.
#ifdef NRG_SYM_ALL
TEST(mk_sym, QJ_requires_three_channels) {
  Params P3;
  EXPECT_NO_THROW(set_symmetry<double>(P3, "QJ", 3));
  for (const auto channels : {1, 2}) {
    Params P;
    EXPECT_THROW(set_symmetry<double>(P, "QJ", channels), std::runtime_error) << channels;
  }
}
#endif

// With substeps=true, P.combs counts the states of a single channel. Only the symmetry types that implement
// recalc_irreduc_substeps() are accepted.
template <typename S> int check_substeps_support() {
  const std::vector<std::string> with_substeps{"QS", "QSZ", "SPSU2", "SPU1", "P", "NONE", "U1", "QSZ4", "SPZ4"};
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    const std::string name(t.name);
    {
      Params P;
      try {
        static_cast<void>(get<S>(name, P));
      } catch (const std::runtime_error &) { continue; } // not compiled in
    }
    const bool supported = std::find(with_substeps.begin(), with_substeps.end(), name) != with_substeps.end();
    Params P;
    P.substeps = true;
    if (supported) {
      EXPECT_NO_THROW(set_symmetry<S>(P, name, 2)) << name;
    } else {
      EXPECT_THROW(set_symmetry<S>(P, name, 2), std::invalid_argument) << name;
    }
    ++checked;
  }
  return checked;
}

TEST(mk_sym, substeps_support) {
  EXPECT_GE(check_substeps_support<double>() + check_substeps_support<std::complex<double>>(), 4);
}

TEST(mk_sym, perchannel) {
  EXPECT_EQ(perchannel("QS"), 1u);
  EXPECT_EQ(perchannel("QSZ"), 1u);
  EXPECT_EQ(perchannel("U1"), 2u);
  EXPECT_EQ(perchannel("NONE"), 4u);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
