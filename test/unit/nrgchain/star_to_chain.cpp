#include <gtest/gtest.h>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "star-to-chain.hpp"

namespace {

using NRG::StarPoint;
using NRG::scalar_star_to_chain;

std::vector<StarPoint> flat_star(const double lambda, const std::size_t mmax) {
  std::vector<StarPoint> star;
  const auto factor = (1.0 - 1.0 / lambda) / std::log(lambda);
  for (std::size_t m = 0; m <= mmax; ++m) {
    const auto energy = factor * std::pow(lambda, -static_cast<double>(m));
    const auto amplitude = std::sqrt((1.0 - 1.0 / lambda) / 2.0) * std::pow(lambda, -static_cast<double>(m) / 2.0);
    star.push_back({energy, amplitude});
    star.push_back({-energy, amplitude});
  }
  return star;
}

double flat_hopping(const double lambda, const std::size_t index) {
  const auto n = static_cast<double>(index);
  const auto alpha = std::log(lambda);
  return -std::expm1(-alpha) / alpha * std::exp(-n * alpha / 2.0) * -std::expm1(-(n + 1.0) * alpha)
         / std::sqrt(std::expm1(-(2.0 * n + 1.0) * alpha) * std::expm1(-(2.0 * n + 3.0) * alpha));
}

// Independent forward residual-norm Lanczos, with enough precision to protect
// 100 steps of the graded recurrence. It sees exactly the same rounded inputs.
NRG::ScalarChain reference_chain(const std::vector<StarPoint> &star, const std::size_t count) {
  using Real = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<1200>>;
  std::vector<Real> q, prev(star.size(), 0), residual(star.size());
  Real norm2 = 0;
  for (const auto &point : star) {
    q.emplace_back(point.amplitude);
    norm2 += q.back() * q.back();
  }
  const Real norm = sqrt(norm2);
  for (auto &value : q) value /= norm;
  Real beta = 0;
  NRG::ScalarChain chain{{}, {}, false};
  for (std::size_t n = 0; n < count; ++n) {
    Real alpha = 0;
    for (std::size_t i = 0; i < star.size(); ++i) alpha += Real(star[i].energy) * q[i] * q[i];
    norm2 = 0;
    for (std::size_t i = 0; i < star.size(); ++i) {
      residual[i] = (Real(star[i].energy) - alpha) * q[i] - beta * prev[i];
      norm2 += residual[i] * residual[i];
    }
    beta = sqrt(norm2);
    chain.zeta.push_back(static_cast<double>(alpha));
    chain.xi.push_back(static_cast<double>(beta));
    prev = q;
    for (std::size_t i = 0; i < star.size(); ++i) q[i] = residual[i] / beta;
  }
  return chain;
}

} // namespace

TEST(ScalarStarToChain, analytic_long_flat_chains) { // NOLINT
  for (const auto lambda : {1.1, 2.0, 4.0}) {
    const std::size_t count = lambda == 4.0 ? 900 : 101;
    const auto result = scalar_star_to_chain(flat_star(lambda, 500 + (lambda == 1.1 ? 300 : 0)), count);
    ASSERT_EQ(result.xi.size(), count);
    ASSERT_EQ(result.zeta.size(), count);
    EXPECT_FALSE(result.terminated);
    for (std::size_t n = 0; n < count; ++n) {
      SCOPED_TRACE(::testing::Message() << "Lambda=" << lambda << " n=" << n);
      const auto expected = flat_hopping(lambda, n);
      EXPECT_GT(result.xi[n], 0.0);
      EXPECT_NEAR(result.xi[n] / expected, 1.0, 5e-13);
      EXPECT_NEAR(result.zeta[n] / expected, 0.0, 5e-13);
    }
  }
}

TEST(ScalarStarToChain, asymmetric_and_gapped_stars_against_high_precision) { // NOLINT
  for (const bool gapped : {false, true}) {
    auto star = flat_star(2.0, gapped ? 40 : 200);
    for (std::size_t m = 0; m < star.size() / 2; ++m) {
      const auto phase = static_cast<double>(m);
      star[2 * m].energy *= 1.0 + 0.2 * std::sin(phase);
      star[2 * m + 1].energy *= 0.9 + 0.1 * std::cos(phase);
      star[2 * m].amplitude *= std::sqrt(1.0 + 0.3 * std::cos(phase));
      star[2 * m + 1].amplitude *= std::sqrt(0.5 + 0.1 * std::sin(phase));
      if (gapped) {
        star[2 * m].energy += 0.1;
        star[2 * m + 1].energy -= 0.1;
      }
    }
    const std::size_t count = gapped ? 30 : 101;
    const auto expected = reference_chain(star, count);
    const auto actual = scalar_star_to_chain(star, count);
    for (std::size_t n = 0; n < count; ++n) {
      SCOPED_TRACE(::testing::Message() << "gapped=" << gapped << " n=" << n);
      EXPECT_NEAR(actual.xi[n] / expected.xi[n], 1.0, 2e-12);
      const auto local_scale = std::max({std::abs(expected.zeta[n]), expected.xi[n], n > 0 ? expected.xi[n - 1] : 0.0});
      EXPECT_NEAR((actual.zeta[n] - expected.zeta[n]) / local_scale, 0.0, 2e-12);
    }
  }
}

TEST(ScalarStarToChain, prefixes_moments_resolvent_and_ordering) { // NOLINT
  const std::vector<StarPoint> star{{0.9, 0.4}, {-0.8, 0.5}, {0.35, 0.2}, {-0.1, 0.7}, {0.0, 0.1}};
  const auto full = scalar_star_to_chain(star, star.size());
  ASSERT_TRUE(full.terminated);
  EXPECT_EQ(full.xi.back(), 0.0);
  for (std::size_t count = 1; count < star.size(); ++count) {
    const auto prefix = scalar_star_to_chain(star, count);
    EXPECT_FALSE(prefix.terminated);
    for (std::size_t n = 0; n < count; ++n) {
      EXPECT_EQ(prefix.zeta[n], full.zeta[n]);
      EXPECT_EQ(prefix.xi[n], full.xi[n]);
    }
  }
  double mass = 0.0;
  for (const auto &p : star) mass += p.amplitude * p.amplitude;
  std::vector<double> vector(star.size(), 0.0);
  vector[0] = 1.0;
  for (int k = 0; k <= 10; ++k) {
    double moment = 0.0;
    for (const auto &p : star) moment += p.amplitude * p.amplitude * std::pow(p.energy, k) / mass;
    EXPECT_NEAR(vector[0], moment, 5e-15);
    auto next = vector;
    for (std::size_t i = 0; i < star.size(); ++i) {
      next[i] = full.zeta[i] * vector[i];
      if (i > 0) next[i] += full.xi[i - 1] * vector[i - 1];
      if (i + 1 < star.size()) next[i] += full.xi[i] * vector[i + 1];
    }
    vector = next;
  }
  for (const auto z : {std::complex<double>(0.0, 0.001), std::complex<double>(0.4, 0.3)}) {
    std::complex<double> expected = 0.0, actual = 0.0;
    for (const auto &p : star) expected += (p.amplitude * p.amplitude / mass) / (z - p.energy);
    for (std::size_t i = star.size(); i-- > 0;) actual = 1.0 / (z - full.zeta[i] - full.xi[i] * full.xi[i] * actual);
    EXPECT_LT(std::abs(actual - expected) / std::abs(expected), 1e-13);
  }
  auto reversed = star;
  std::reverse(reversed.begin(), reversed.end());
  const auto other = scalar_star_to_chain(reversed, star.size());
  for (std::size_t n = 0; n < star.size(); ++n) {
    EXPECT_NEAR(other.xi[n], full.xi[n], 2e-15);
    EXPECT_NEAR(other.zeta[n], full.zeta[n], 2e-15);
  }
}

TEST(ScalarStarToChain, scaling_and_unsquared_subnormal_inputs) { // NOLINT
  const std::vector<StarPoint> original{{1.0, 0.3}, {-0.5, 0.8}, {0.1, 0.2}, {-0.01, 0.1}};
  const auto reference = scalar_star_to_chain(original, 3);
  for (const double energy_scale : {1e-200, 1.0, 1e200}) {
    for (const double amplitude_scale : {1e-250, 1.0, 1e250}) {
      auto star = original;
      for (auto &p : star) {
        p.energy *= energy_scale;
        p.amplitude *= amplitude_scale;
      }
      const auto chain = scalar_star_to_chain(star, 3);
      for (std::size_t n = 0; n < 3; ++n) {
        EXPECT_NEAR(chain.xi[n] / energy_scale / reference.xi[n], 1.0, 2e-14);
        EXPECT_NEAR((chain.zeta[n] / energy_scale - reference.zeta[n]) / reference.xi[n], 0.0, 2e-14);
      }
    }
  }
  const std::vector<StarPoint> subnormal{{1e-310, 1e-310}, {-1e-310, 1e-310}};
  const auto tiny = scalar_star_to_chain(subnormal, 2);
  EXPECT_NEAR(tiny.xi[0] / 1e-310, 1.0, 2e-13);
  EXPECT_NEAR(tiny.zeta[0] / 1e-310, 0.0, 2e-13);
  const auto least = std::numeric_limits<double>::denorm_min();
  for (const double multiple : {1.0, 2.0, 3.0, 4.0, 7.0, 16.0, 255.0, 1024.0}) {
    const auto amplitude = multiple * least;
    const auto mixed = scalar_star_to_chain(std::vector<StarPoint>{{0.0, 1.0}, {1.0, amplitude}, {-1.0, amplitude}}, 3);
    EXPECT_EQ(mixed.xi[0], std::hypot(amplitude, amplitude));
    EXPECT_NEAR(mixed.xi[1], 1.0, 5e-16);
    EXPECT_EQ(mixed.xi[2], 0.0);
    for (const auto onsite : mixed.zeta) EXPECT_NEAR(onsite, 0.0, 5e-16);
  }
  for (const auto scale : {least, 2.0 * least, 3.0 * least, std::numeric_limits<double>::max() / 4.0}) {
    const std::vector<StarPoint> equal{{1.0, scale}, {-1.0, scale}};
    const auto symmetric = scalar_star_to_chain(equal, 2);
    EXPECT_NEAR(symmetric.xi[0], 1.0, 5e-16);
    EXPECT_NEAR(symmetric.zeta[0], 0.0, 5e-16);
    const auto unequal = scalar_star_to_chain(std::vector<StarPoint>{{1.0, scale}, {-1.0, 2.0 * scale}, {0.1, 3.0 * scale}}, 3);
    const auto normal = scalar_star_to_chain(std::vector<StarPoint>{{1.0, 1.0}, {-1.0, 2.0}, {0.1, 3.0}}, 3);
    for (std::size_t n = 0; n < 3; ++n) {
      EXPECT_NEAR(unequal.xi[n], normal.xi[n], 5e-16);
      EXPECT_NEAR(unequal.zeta[n], normal.zeta[n], 5e-16);
    }
  }
}

TEST(ScalarStarToChain, exact_support_and_invalid_inputs) { // NOLINT
  const std::vector<StarPoint> repeated{{0.9, 0.5}, {-0.7, 0.5}, {0.9, 0.5}, {-0.7, 0.5}, {0.2, 0.0}};
  const auto chain = scalar_star_to_chain(repeated, 2);
  ASSERT_TRUE(chain.terminated);
  EXPECT_NEAR(chain.xi[0], 0.8, 2e-16);
  EXPECT_EQ(chain.xi[1], 0.0);
  EXPECT_NEAR(chain.zeta[0], 0.1, 2e-16);
  EXPECT_NEAR(chain.zeta[1], 0.1, 2e-16);
  EXPECT_THROW(scalar_star_to_chain(repeated, 3), std::invalid_argument);
  const auto singleton = scalar_star_to_chain(std::vector<StarPoint>{{0.25, 1.0}}, 1);
  EXPECT_EQ(singleton.zeta[0], 0.25);
  EXPECT_EQ(singleton.xi[0], 0.0);
  EXPECT_TRUE(singleton.terminated);

  const auto next = std::nextafter(1.0, 2.0);
  const auto close = scalar_star_to_chain(std::vector<StarPoint>{{1.0, 1.0}, {next, 1.0}}, 2);
  EXPECT_GT(close.xi[0], 0.0); // Distinct representable nodes must not be coalesced.
  const auto inf = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(scalar_star_to_chain({}, 1), std::invalid_argument);
  EXPECT_THROW(scalar_star_to_chain(repeated, 0), std::invalid_argument);
  for (const auto point : {StarPoint{0.0, 0.0}, StarPoint{0.0, -1.0}, StarPoint{inf, 1.0}, StarPoint{0.0, nan}})
    EXPECT_THROW(scalar_star_to_chain(std::vector<StarPoint>{point}, 1), std::invalid_argument);
  const auto max = std::numeric_limits<double>::max();
  const auto merged = scalar_star_to_chain(std::vector<StarPoint>{{1.0, max}, {1.0, max}, {-1.0, max}, {-1.0, max}}, 2);
  EXPECT_NEAR(merged.xi[0], 1.0, 5e-16);
  EXPECT_THROW(scalar_star_to_chain(std::vector<StarPoint>{{1.0, max}, {-1.0, std::numeric_limits<double>::denorm_min()}}, 1),
               std::runtime_error);
  EXPECT_THROW(scalar_star_to_chain(std::vector<StarPoint>{{max, 1.0}, {-max, 1.0}}, 1), std::runtime_error);
}
