#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include "test_clean.hpp"
#include <coef.hpp>
#include <read-input.hpp>

using namespace NRG;

TEST(Coef, parse) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag_clean<double>(P, Sym);
  Coef<double> coef(P);
  std::string str = 
    "1\n"
    "0.54528747084262258072\n"
    "0.41550946829175445321\n"
    "1\n"
    "0.e-999\n"
    "0.e-998\n";
  std::istringstream ss(str);
  EXPECT_EQ(P.coefchannels, 1);
  coef.xi.read(ss, P.coefchannels);
  coef.zeta.read(ss, P.coefchannels);

  EXPECT_EQ(coef.xi.nr_tabs(), 1);
  EXPECT_EQ(coef.zeta.nr_tabs(), 1);
  EXPECT_EQ(coef.xi.max(0), 1);
  EXPECT_EQ(coef.zeta.max(0), 1);

  EXPECT_EQ(coef.xi(0,0), 0.54528747084262258072);
  EXPECT_EQ(coef.xi(1,0), 0.41550946829175445321); // first index: N, second index: alpha (channel number)
  EXPECT_EQ(coef.zeta(0,0), 0);
  EXPECT_EQ(coef.zeta(1,0), 0);
  
  coef.xi.set(0, 0, 42.);
  EXPECT_EQ(coef.xi(0,0), 42.);
}

TEST(Coef, truncated_matrix_table_header_throws) { // NOLINT
  std::istringstream ss("1 2");
  EXPECT_THROW(read_matrix_table<double>(ss), std::runtime_error);
}

namespace {

template<typename S>
Coef<S> equal_bath_coef(const Params &P) {
  Coef<S> coef(P);
  std::ostringstream zeros;
  for (size_t ch = 0; ch < P.coefchannels; ++ch) zeros << "3 0 0 0 0\n";
  for (auto *table : {&coef.xi, &coef.zeta, &coef.ep, &coef.em, &coef.u0p, &coef.u0m}) {
    std::istringstream input(zeros.str());
    table->read(input, P.coefchannels);
  }
  const double hoppings[] = {0.5, 0.25, 2e-20, 1e-20};
  for (size_t ch = 0; ch < P.coefchannels; ++ch)
    for (size_t n = 0; n < 4; ++n) {
      coef.xi.set(n, ch, hoppings[n]);
      coef.zeta.set(n, ch, 0.25 * hoppings[n]);
      // Proportional star weights, with a common nonzero onsite shift.
      const auto energy = 1.0 - 0.25 * double(n);
      coef.ep.set(n, ch, energy + 0.125);
      coef.em.set(n, ch, energy - 0.125);
      coef.u0p.set(n, ch, std::ldexp(1.0, int(ch) - int(n)));
      coef.u0m.set(n, ch, std::ldexp(1.0, int(ch) - int(n)));
    }
  return coef;
}

template<typename Action>
void expect_bath_error(Action action, const std::string &sym, const std::string &coefficient, const size_t ch, const size_t n) {
  try {
    action();
    FAIL() << "Accepted unsupported bath chains";
  } catch (const std::invalid_argument &error) {
    const std::string message = error.what();
    EXPECT_NE(message.find(sym + ":"), std::string::npos) << message;
    EXPECT_NE(message.find("only equal normalized bath chains"), std::string::npos) << message;
    EXPECT_NE(message.find(coefficient + " channel " + std::to_string(ch)), std::string::npos) << message;
    EXPECT_NE(message.find("coefficient index " + std::to_string(n)), std::string::npos) << message;
  }
}

template<typename S>
std::string bath_input(const Params &P, const Coef<S> &coef) {
  std::ostringstream data;
  data << std::setprecision(std::numeric_limits<double>::max_digits10);
  data << "#!9\n# symtype " << P.symtype.value() << "\n3\n" << coef.xi.max(0) << "\n1\n";
  // One seed eigenstate, no boundary-operator matrix elements: parsing only, no NRG iterations.
  data << "0 1 0\n1\n0\n";
  for (size_t ch = 0; ch < 3; ++ch) data << "f " << ch << " 0\n0\n";
  data << "z\n";
  for (const auto *table : {&coef.xi, &coef.zeta})
    for (size_t ch = 0; ch < P.coefchannels; ++ch) {
      data << table->max(ch) << '\n';
      for (size_t n = 0; n <= table->max(ch); ++n) data << (*table)(n, ch) << '\n';
    }
  data << "T\n";
  for (const auto *table : {&coef.ep, &coef.em, &coef.u0p, &coef.u0m})
    for (size_t ch = 0; ch < P.coefchannels; ++ch) {
      data << table->max(ch) << '\n';
      for (size_t n = 0; n <= table->max(ch); ++n) data << (*table)(n, ch) << '\n';
    }
  return data.str();
}

template<typename S>
class EqualBathChains : public ::testing::Test {
  const mp_bitcnt_t saved_precision = mpf_get_default_prec();
 protected:
  void TearDown() override { mpf_set_default_prec(saved_precision); }
};

using BathScalarTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(EqualBathChains, BathScalarTypes);

TYPED_TEST(EqualBathChains, equal_nonzero_onsites_and_local_roundoff_pass) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"}) {
    SCOPED_TRACE(sym);
    Params P;
    P.symtype = sym;
    P.set_channels_and_combs(3);
    for (const double scale : {1e-100, 1.0, 1e100}) {
      SCOPED_TRACE(scale);
      auto coef = equal_bath_coef<TypeParam>(P);
      for (auto *table : {&coef.xi, &coef.zeta})
        for (size_t ch = 0; ch < 3; ++ch)
          for (size_t n = 0; n <= table->max(ch); ++n) table->set(n, ch, (*table)(n, ch) * scale);
      EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
      const auto n = coef.xi.max(0);
      const auto roundoff = 4 * std::numeric_limits<double>::epsilon();
      coef.xi.set(n, 1, coef.xi(n, 1) * (1 + roundoff));
      coef.zeta.set(0, 2, coef.zeta(0, 2) * (1 + roundoff));
      for (size_t ch = 0; ch < 3; ++ch) coef.zeta.set(n, ch, 0.0);
      coef.zeta.set(n, 2, roundoff * coef.xi(n, 0));
      EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
      for (size_t ch = 0; ch < 3; ++ch) coef.xi.set(n, ch, 0.0);
      coef.zeta.set(n, 2, roundoff * coef.xi(n - 1, 0));
      EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
    }
  }
}

TYPED_TEST(EqualBathChains, rejects_late_differences_including_max_index) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"})
      for (size_t ch = 1; ch < 3; ++ch) {
        Params P;
        P.symtype = sym;
        P.set_channels_and_combs(3);
        auto coef = equal_bath_coef<TypeParam>(P);
        auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
        ASSERT_EQ(table.max(0), 3U); // Four coefficients, not three.
        table.set(3, ch, table(3, ch) + TypeParam(1e-22));
        expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, ch, 3);
      }
}

TYPED_TEST(EqualBathChains, rejects_mismatched_lengths_and_missing_channels) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"}) {
      Params P;
      P.symtype = sym;
      P.set_channels_and_combs(3);
      for (size_t ch = 1; ch < 3; ++ch) {
        auto coef = equal_bath_coef<TypeParam>(P);
        auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
        table.set(4, ch, 0.0);
        expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, ch, 4);
        std::ostringstream shortened;
        for (size_t i = 0; i < 3; ++i) shortened << (i == ch ? "2 0 0 0\n" : "3 0 0 0 0\n");
        std::istringstream input(shortened.str());
        table.read(input, 3);
        expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, ch, 3);
      }
      for (const size_t count : {0U, 2U}) {
        auto coef = equal_bath_coef<TypeParam>(P);
        auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
        std::istringstream input("0 0\n0 0\n");
        table.read(input, count);
        expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, count, 0);
      }
    }
}

TYPED_TEST(EqualBathChains, rejects_nonfinite_values_in_every_channel) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"})
      for (size_t ch = 0; ch < 3; ++ch)
        for (const auto bad : {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
          Params P;
          P.symtype = sym;
          P.set_channels_and_combs(3);
          auto coef = equal_bath_coef<TypeParam>(P);
          auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
          table.set(3, ch, TypeParam(bad));
          expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, ch, 3);
          if constexpr (is_complex<TypeParam>::value) {
            table.set(3, ch, TypeParam(0.0, bad));
            expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, ch, 3);
          }
        }
}

TYPED_TEST(EqualBathChains, ignores_unused_arrays_and_extra_channels) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"}) {
    Params P;
    P.symtype = sym;
    P.polarized = true;
    P.set_channels_and_combs(3);
    auto coef = equal_bath_coef<TypeParam>(P);
    for (auto *table : {&coef.xi, &coef.zeta})
      for (size_t ch = 3; ch < P.coefchannels; ++ch) table->set(4 + ch, ch, std::numeric_limits<double>::quiet_NaN());
    for (auto *table : {&coef.xiR, &coef.zetaR, &coef.delta, &coef.kappa, &coef.ep, &coef.em, &coef.u0p, &coef.u0m}) {
      std::istringstream input("0 0\n");
      table->read(input, 1);
      table->set(0, 0, std::numeric_limits<double>::infinity());
    }
    EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
  }
}

TYPED_TEST(EqualBathChains, leaves_other_symmetries_unchanged) { // NOLINT
  for (const auto sym : {"QS", "QSZ", "QST", "SPSU2T"}) {
    Params P;
    P.symtype = sym;
    P.set_channels_and_combs(3);
    auto coef = equal_bath_coef<TypeParam>(P);
    coef.xi.set(3, 1, 42.0);
    coef.zeta.set(4, 2, std::numeric_limits<double>::quiet_NaN());
    EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
  }
}

TEST(Coef, equal_bath_chains_compare_complex_components) { // NOLINT
  using S = std::complex<double>;
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"}) {
      Params P;
      P.symtype = sym;
      P.set_channels_and_combs(3);
      auto coef = equal_bath_coef<S>(P);
      auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
      for (size_t ch = 0; ch < 3; ++ch) table.set(3, ch, S(1e-20, 1e-20));
      EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
      table.set(3, 1, S(1e-20, std::nextafter(1e-20, 1.0)));
      EXPECT_NO_THROW(validate_equal_bath_chains(coef, P));
      table.set(3, 1, S(1e-20, 1.01e-20));
      expect_bath_error([&] { validate_equal_bath_chains(coef, P); }, sym, name, 1, 3);
    }
}

template<typename S>
class EqualBathInput : public EqualBathChains<S> {
 protected:
  void SetUp() override {
#ifndef NRG_SYM_ALL
    GTEST_SKIP() << "QSTZ/QSZTZ InputData integration requires SYM_ALL.";
#endif
  }
};

TYPED_TEST_SUITE(EqualBathInput, BathScalarTypes);

TYPED_TEST(EqualBathInput, legacy_equal_chains_and_roundoff_pass) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"}) {
    Params P;
    P.symtype = sym;
    P.silent = true;
    P.tri = "old";
    P.set_channels_and_combs(3);
    auto coef = equal_bath_coef<TypeParam>(P);
    coef.xi.set(3, 1, coef.xi(3, 1) * (1 + 4 * std::numeric_limits<double>::epsilon()));
    std::istringstream stream(bath_input(P, coef));
    const InputData<TypeParam> input(P, stream);
    EXPECT_EQ(P.Nmax, 3U);
    EXPECT_EQ(input.coef.xi(3, 1), coef.xi(3, 1)); // No averaging or other coefficient changes.
    EXPECT_EQ(input.coef.zeta(0, 2), TypeParam(0.125));
  }
}

TYPED_TEST(EqualBathInput, legacy_rejects_late_xi_and_zeta_differences) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"}) {
      Params P;
      P.symtype = sym;
      P.silent = true;
      P.tri = "old";
      P.set_channels_and_combs(3);
      auto coef = equal_bath_coef<TypeParam>(P);
      auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
      table.set(3, 2, table(3, 2) + TypeParam(1e-22));
      std::istringstream stream(bath_input(P, coef));
      expect_bath_error([&] { InputData<TypeParam> input(P, stream); }, sym, name, 2, 3);
    }
}

TYPED_TEST(EqualBathInput, legacy_rejects_mismatched_lengths) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"}) {
      Params P;
      P.symtype = sym;
      P.silent = true;
      P.tri = "old";
      P.set_channels_and_combs(3);
      auto coef = equal_bath_coef<TypeParam>(P);
      auto &table = std::string(name) == "xi" ? coef.xi : coef.zeta;
      table.set(4, 1, 0.0);
      std::istringstream stream(bath_input(P, coef));
      expect_bath_error([&] { InputData<TypeParam> input(P, stream); }, sym, name, 1, 4);
    }
}

TYPED_TEST(EqualBathInput, cpp_replaces_placeholders_and_accepts_proportional_stars) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"}) {
    Params P;
    P.symtype = sym;
    P.silent = true;
    P.tri = "cpp";
    P.tridiag_method = "lanczos";
    P.preccpp = 256;
    P.set_channels_and_combs(3);
    auto coef = equal_bath_coef<TypeParam>(P);
    coef.xi.set(3, 1, 42.0); // Must be checked after generation, not before it.
    coef.zeta.set(3, 2, 43.0);
    std::istringstream stream(bath_input(P, coef));
    const InputData<TypeParam> input(P, stream);
    EXPECT_EQ(P.Nmax, 3U);
    EXPECT_NE(input.coef.u0p(0, 0), input.coef.u0p(0, 2));
    for (size_t ch = 0; ch < 3; ++ch)
      for (size_t n = 0; n <= input.coef.xi.max(ch); ++n) {
        EXPECT_NEAR(std::abs(input.coef.xi(n, ch) - input.coef.xi(n, 0)), 0.0, 1e-14);
        EXPECT_NEAR(std::abs(input.coef.zeta(n, ch) - TypeParam(0.125)), 0.0, 1e-14);
      }
  }
}

TYPED_TEST(EqualBathInput, cpp_rejects_generated_xi_and_zeta_differences) { // NOLINT
  for (const auto sym : {"QSTZ", "QSZTZ"})
    for (const auto name : {"xi", "zeta"}) {
      Params P;
      P.symtype = sym;
      P.silent = true;
      P.tri = "cpp";
      P.tridiag_method = "lanczos";
      P.preccpp = 256;
      P.set_channels_and_combs(3);
      auto coef = equal_bath_coef<TypeParam>(P);
      for (size_t n = 0; n <= coef.ep.max(2); ++n) {
        if (std::string(name) == "xi") {
          coef.ep.set(n, 2, 1.25 * coef.ep(n, 2));
          coef.em.set(n, 2, 1.25 * coef.em(n, 2));
        } else {
          coef.ep.set(n, 2, coef.ep(n, 2) + TypeParam(0.0625));
          coef.em.set(n, 2, coef.em(n, 2) - TypeParam(0.0625));
        }
      }
      std::istringstream stream(bath_input(P, coef));
      expect_bath_error([&] { InputData<TypeParam> input(P, stream); }, sym, name, 2, 0);
    }
}

} // namespace

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
