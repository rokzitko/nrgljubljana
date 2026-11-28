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
#include <type_traits>
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

// symmetry_traits (symmetries.def) lists every symmetry type of get(), with the properties of its class. Types not
// compiled in (SYM_MORE/SYM_ALL off, or complex-only types for real S) throw from get() and are skipped.
template <typename S> int check_symmetry_traits() {
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    Params P;
    try {
      auto sym = get<S>(std::string(t.name), P);
      EXPECT_TRUE(t.compiled) << t.name;
      EXPECT_EQ(sym->isfield(), t.field) << t.name;
      EXPECT_EQ(sym->uses_pairing(), t.pairing) << t.name;
      ++checked;
    } catch (const std::runtime_error &) {
      EXPECT_TRUE(!t.compiled || (t.complex_only && !std::is_same_v<S, std::complex<double>>)) << t.name;
    }
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
    EXPECT_THROW(set_symmetry<double>(P, "QJ", channels), std::invalid_argument) << channels;
  }
}
#endif

// Every symmetry type loads its coefficient tables for the supported numbers of channels (up to five for the types
// with tables generated at run time) and rejects the others.
template <typename S> int check_channel_support() {
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    const std::string name(t.name);
    if (!t.compiled || (t.complex_only && !std::is_same_v<S, std::complex<double>>)) continue;
    for (std::size_t channels = 1; channels <= 5; channels++) {
      Params P;
      if (t.supports_channels(channels)) {
        EXPECT_NO_THROW(set_symmetry<S>(P, name, channels)) << name << " " << channels;
      } else {
        EXPECT_THROW(set_symmetry<S>(P, name, channels), std::invalid_argument) << name << " " << channels;
      }
    }
    ++checked;
  }
  return checked;
}

TEST(mk_sym, channel_support) {
  EXPECT_GE(check_channel_support<double>() + check_channel_support<std::complex<double>>(), 4);
}

TEST(mk_sym, polarized_support) {
  Params P;
  P.polarized = true;
  EXPECT_NO_THROW(set_symmetry<double>(P, "QSZ", 1));
  Params P2;
  P2.polarized = true;
  EXPECT_THROW(set_symmetry<double>(P2, "QS", 1), std::invalid_argument);
}

// With substeps=true, P.combs counts the states of a single channel. Only the symmetry types that implement
// recalc_irreduc_substeps() are accepted.
template <typename S> int check_substeps_support() {
  const std::vector<std::string> with_substeps{"QS", "QSZ", "SPSU2", "SPU1", "P", "NONE", "U1", "QSZ4", "SPZ4"};
  for (const auto &t : symmetry_traits)
    EXPECT_EQ(t.substeps, std::find(with_substeps.begin(), with_substeps.end(), std::string(t.name)) != with_substeps.end()) << t.name;
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
