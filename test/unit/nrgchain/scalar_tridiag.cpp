#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "test_common.hpp"
#include <tridiag.hpp>

namespace {

template<typename S>
Coef<S> make_coef(const Params &P, const std::vector<std::vector<StarPoint>> &stars) {
  Coef<S> coef(P);
  std::ostringstream zeros;
  for (size_t ch = 0; ch < P.coefchannels; ++ch) zeros << "0\n0\n";
  for (auto *table : {&coef.ep, &coef.em, &coef.u0p, &coef.u0m, &coef.xi, &coef.zeta}) {
    std::istringstream input(zeros.str());
    table->read(input, P.coefchannels);
  }
  for (size_t ch = 0; ch < stars.size(); ++ch) {
    const auto &star = stars[ch];
    for (size_t m = 0; m < star.size() / 2; ++m) {
      coef.ep.set(m, ch, star[2 * m].energy);
      coef.em.set(m, ch, -star[2 * m + 1].energy);
      coef.u0p.set(m, ch, star[2 * m].amplitude);
      coef.u0m.set(m, ch, star[2 * m + 1].amplitude);
    }
  }
  return coef;
}

template<typename S>
class ScalarTridiag : public ::testing::Test {
  const mp_bitcnt_t saved_precision = mpf_get_default_prec();
 protected:
  void TearDown() override { mpf_set_default_prec(saved_precision); }
};

using ScalarTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(ScalarTridiag, ScalarTypes);

TYPED_TEST(ScalarTridiag, matches_gmp_multiple_channels_and_band_scales) { // NOLINT
  Params P;
  P.tri = "cpp";
  P.set_channels_and_combs(2);
  P.preccpp = 512;
  P.validate();
  std::vector<std::vector<StarPoint>> stars(P.coefchannels);
  for (size_t ch = 0; ch < stars.size(); ++ch) {
    for (size_t m = 0; m < 24; ++m) {
      const auto energy = std::pow(2.0, -static_cast<double>(m));
      const auto phase = static_cast<double>(m) + 0.4 * static_cast<double>(ch);
      stars[ch].push_back({energy * (0.8 + 0.1 * std::cos(phase)), std::sqrt(energy * (0.6 + 0.1 * std::sin(phase)))});
      stars[ch].push_back({-energy * (0.7 + 0.1 * std::sin(phase)), std::sqrt(energy * (0.3 + 0.1 * std::cos(phase)))});
    }
  }
  constexpr size_t nmax = 7;
  auto legacy = make_coef<TypeParam>(P, stars);
  Tridiag<TypeParam>(legacy, nmax, P); // Default method is the actual GMP integration path.

  mpf_set_default_prec(192);
  const auto precision = mpf_get_default_prec();
  P.tridiag_method = "rkpw";
  P.preccpp = 0;
  for (const double scale : {1e-150, 1.0, 1e150}) {
    SCOPED_TRACE(scale);
    P.bandrescale = scale;
    P.validate();
    auto actual = make_coef<TypeParam>(P, stars);
    Tridiag<TypeParam>(actual, nmax, P);
    EXPECT_EQ(mpf_get_default_prec(), precision);
    ASSERT_EQ(actual.xi.nr_tabs(), P.coefchannels);
    ASSERT_EQ(actual.zeta.nr_tabs(), P.coefchannels);
    for (size_t ch = 0; ch < P.coefchannels; ++ch) {
      ASSERT_EQ(actual.xi.max(ch), nmax);
      ASSERT_EQ(actual.zeta.max(ch), nmax);
      for (size_t n = 0; n <= nmax; ++n) {
        SCOPED_TRACE(::testing::Message() << "ch=" << ch << " n=" << n);
        const double xi = std::real(legacy.xi(n, ch));
        ASSERT_GT(xi, 0.0);
        EXPECT_NEAR(std::real(actual.xi(n, ch)) / scale / xi, 1.0, 2e-12);
        EXPECT_NEAR((std::real(actual.zeta(n, ch)) / scale - std::real(legacy.zeta(n, ch))) / xi, 0.0, 2e-12);
        EXPECT_EQ(std::imag(actual.xi(n, ch)), 0.0);
        EXPECT_EQ(std::imag(actual.zeta(n, ch)), 0.0);
      }
    }
  }
}

TYPED_TEST(ScalarTridiag, effective_support_and_terminal_hopping) { // NOLINT
  Params P;
  P.tri = "cpp";
  P.tridiag_method = "rkpw";
  P.preccpp = 0;
  P.bandrescale = 2.5;
  P.set_channels_and_combs(1);
  P.validate();
  const std::vector<std::vector<StarPoint>> stars{{{0.9, 0.5}, {-0.7, 0.5}, {0.9, 0.5}, {-0.7, 0.5},
                                                {0.2, 0.0}, {-0.1, 0.0}}};
  const auto precision = mpf_get_default_prec();
  auto actual = make_coef<TypeParam>(P, stars);
  Tridiag<TypeParam>(actual, 1, P);
  ASSERT_EQ(actual.xi.max(0), 1U);
  ASSERT_EQ(actual.zeta.max(0), 1U);
  EXPECT_NEAR(std::real(actual.xi(0, 0)), 2.0, 2e-14);
  EXPECT_EQ(actual.xi(1, 0), TypeParam(0.0));
  EXPECT_NEAR(std::real(actual.zeta(0, 0)), 0.25, 2e-14);
  EXPECT_NEAR(std::real(actual.zeta(1, 0)), 0.25, 2e-14);

  auto overlong = make_coef<TypeParam>(P, stars);
  overlong.xi.set(0, 0, 42.0);
  overlong.zeta.set(0, 0, 43.0);
  try {
    Tridiag<TypeParam>(overlong, 2, P);
    FAIL() << "Accepted chain longer than effective support";
  } catch (const std::invalid_argument &error) {
    EXPECT_NE(std::string(error.what()).find("effective nonzero support"), std::string::npos);
  }
  EXPECT_EQ(overlong.xi(0, 0), TypeParam(42.0));
  EXPECT_EQ(overlong.zeta(0, 0), TypeParam(43.0));
  EXPECT_EQ(overlong.xi.max(0), 0U);
  EXPECT_EQ(mpf_get_default_prec(), precision);
}

TYPED_TEST(ScalarTridiag, rejects_unrepresentable_scaled_coefficients_before_publication) { // NOLINT
  Params P;
  P.tri = "cpp";
  P.tridiag_method = "rkpw";
  P.set_channels_and_combs(1);
  struct Case {
    double scale;
    std::vector<StarPoint> star;
  };
  const double tiny = std::numeric_limits<double>::denorm_min();
  const std::vector<Case> cases{
    {std::numeric_limits<double>::max(), {{4.0, 0.5}, {-2.0, 0.5}}}, // Hopping overflow.
    {1e200, {{1e200, 1.0}, {-1.0, 0.0}}},                         // Onsite overflow, terminal hopping is valid.
    {tiny, {{0.125, 0.5}, {-0.125, 0.5}}},                        // Hopping underflow.
    {tiny, {{1.125, 0.5}, {-0.875, 0.5}}},                        // Only the onsite underflows.
  };
  for (const auto &test : cases) {
    SCOPED_TRACE(test.scale);
    P.bandrescale = test.scale;
    auto coef = make_coef<TypeParam>(P, {test.star});
    coef.xi.set(0, 0, 42.0);
    coef.zeta.set(0, 0, 43.0);
    EXPECT_THROW(Tridiag<TypeParam>(coef, 0, P), std::runtime_error);
    EXPECT_EQ(coef.xi(0, 0), TypeParam(42.0));
    EXPECT_EQ(coef.zeta(0, 0), TypeParam(43.0));
  }
  for (const auto scale : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    P.bandrescale = scale;
    auto coef = make_coef<TypeParam>(P, {cases.front().star});
    EXPECT_THROW(Tridiag<TypeParam>(coef, 0, P), std::invalid_argument);
  }
  P.bandrescale = tiny;
  auto subnormal = make_coef<TypeParam>(P, {{{1.0, 1.0}, {-1.0, 0.0}}});
  EXPECT_NO_THROW(Tridiag<TypeParam>(subnormal, 0, P));
  EXPECT_EQ(subnormal.zeta(0, 0), TypeParam(tiny));
  EXPECT_EQ(subnormal.xi(0, 0), TypeParam(0.0));
}

TEST(ScalarTridiagParams, defaults_validation_and_reporting) { // NOLINT
  Params P;
  EXPECT_EQ(P.tri.value(), "old");
  EXPECT_EQ(P.tridiag_method.value(), "lanczos");
  EXPECT_EQ(P.preccpp.value(), 2000U);
  EXPECT_NO_THROW(P.validate());
  for (const auto tri : {"old", "cpp"}) {
    P.tri = tri;
    for (const auto method : {"unknown", "RKPW", ""}) {
      P.tridiag_method.set_str(method);
      EXPECT_THROW(P.validate(), std::invalid_argument);
    }
  }
  P.tridiag_method.set_str("rkpw");
  P.preccpp.set_str("0");
  EXPECT_NO_THROW(P.validate());
  std::ostringstream report;
  P.dump(report);
  EXPECT_NE(report.str().find("tridiag_method=rkpw"), std::string::npos);
  P.tridiag_method.set_str("lanczos");
  EXPECT_THROW(P.validate(), std::invalid_argument);
  P.preccpp = 512;
  EXPECT_NO_THROW(P.validate());
}

} // namespace
