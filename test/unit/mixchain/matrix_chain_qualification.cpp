#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include <mixchain/chain.hpp>
#include "reference_lanczos.hpp"
#include "reference_precision.hpp"

using namespace NRG::MixChain;

// QUALIFICATION OF THE ROTATIONS FOR SEVERAL CHANNELS
//
// How accurate the chain stage is on stars of the kind the tool is used for, over a spread of parameters. The
// stars are built here as lists of levels, not by the star stage, so that only the mapping to the chain is tested:
// the rotations in double precision and block Lanczos in multiprecision see the same numbers. The scalar counterpart
// is test/unit/nrgchain/scalar_chain_qualification.cpp, described in test/CHAIN_QUALIFICATION.md.
//
// The reference is block Lanczos at two precisions, which must agree with each other before it is used.
//
// Where a branch accumulates at a finite energy, as at a gap edge, the chain is ill-conditioned as a function of the
// star: the exact chain of a star stored in double moves by far more than rounding when every number of the star is
// changed by one unit in the last place. No algorithm in double can do better than that, so those cases are held to
// a multiple of that sensitivity, measured here, and print an advisory when they exceed the plain target.

namespace {

// The target of the scalar chain qualification.
constexpr double budget = 2e-12;
// How far beyond the one-ulp sensitivity of the reference an ill-conditioned case may be.
constexpr double sensitivity_margin = 10.0;
// The two precisions of the reference must agree to this.
constexpr double reference_agreement = 1e-30;

constexpr unsigned low_digits  = 50;
constexpr unsigned high_digits = 100;

struct Case {
  const char *name;
  int channels;
  bool complex_data;
  double lambda;
  unsigned int nmax;
  double weight_ratio; // weight of the positive branch over that of the negative one
  double turn;         // how far the eigenvectors turn from one interval to the next, in radians
  double gap;          // the last eigenvalue branch accumulates at +-gap instead of at zero
};

template<typename S0> S0 phase(const double angle) {
  if constexpr (is_complex_v<S0>)
    return std::polar(1.0, angle);
  else
    return 1.0;
}

// The orthonormal coupling directions of one interval, as the columns of a unitary that turns with m.
template<typename S0> Matrix<S0> directions(const Case &c, const int m, const bool negative) {
  const auto n = c.channels;
  Matrix<S0> u = Matrix<S0>::Identity(n, n);
  for (int i = 0; i + 1 < n; i++) {
    const auto angle = 0.4 + 0.3 * i + (negative ? 0.5 : 0.0) + c.turn * m + 0.3 * std::sin(0.7 * m + i);
    Matrix<S0> g     = Matrix<S0>::Identity(n, n);
    g(i, i)          = std::cos(angle);
    g(i, i + 1)      = -std::sin(angle);
    g(i + 1, i)      = std::sin(angle);
    g(i + 1, i + 1)  = std::cos(angle);
    u                = (u * g).eval();
  }
  for (int j = 0; j < n; j++) u.col(j) *= phase<S0>(0.9 * std::cos(0.4 * m) + 0.3 * j + (negative ? 0.3 : 0.0));
  for (int r = 0; r < n; r++) u.row(r) *= phase<S0>(0.7 * r);
  return u;
}

// A star graded like Wilson's, with mMAX = 2 Nmax as the tool takes by default.
template<typename S0> Star<S0> star_of(const Case &c) {
  Star<S0> star;
  star.channels = c.channels;
  star.mMAX     = 2 * c.nmax;
  star.Lambda   = c.lambda;
  star.z        = 1.0;
  for (const bool negative : {false, true}) {
    for (int m = 0; m <= static_cast<int>(star.mMAX); m++) {
      const auto u     = directions<S0>(c, m, negative);
      const auto decay = std::pow(c.lambda, -m);
      for (int i = 0; i < c.channels; i++) {
        const auto shift  = i == c.channels - 1 ? c.gap : 0.0;
        const auto energy = decay * (0.6 + 0.3 * i / c.channels + 0.05 * std::cos(1.1 * m + i)) * (1.0 - shift) + shift;
        const auto weight = decay * (1.0 - 1.0 / c.lambda) * (0.5 + 0.3 * i / c.channels + 0.1 * std::sin(1.3 * m + 2 * i))
                            * (negative ? 1.0 / c.weight_ratio : 1.0);
        StarLevel<S0> level;
        level.m        = m;
        level.sign     = negative ? Sign::NEG : Sign::POS;
        level.branch   = i;
        level.energy   = negative ? -energy : energy;
        level.coupling = std::sqrt(weight) * u.col(i);
        star.levels.push_back(level);
      }
    }
  }
  return star;
}

// The same star with every number moved by one unit in the last place, up or down.
template<typename S0> Star<S0> one_ulp_away(Star<S0> star) {
  std::mt19937 generator(20261005);
  const auto nudge = [&generator](const double x) { return std::nextafter(x, generator() % 2 ? 2.0 * x : 0.0); };
  for (auto &level : star.levels) {
    level.energy = nudge(level.energy);
    for (Eigen::Index i = 0; i < level.coupling.size(); i++) {
      if constexpr (is_complex_v<S0>)
        level.coupling(i) = S0(nudge(level.coupling(i).real()), nudge(level.coupling(i).imag()));
      else
        level.coupling(i) = nudge(level.coupling(i));
    }
  }
  return star;
}

template<typename S0, unsigned Digits>
using Wide = std::conditional_t<is_complex_v<S0>, WideComplex<Digits>, WideReal<Digits>>;

template<typename S0, unsigned Digits> auto lanczos(const Star<S0> &star, const unsigned int nmax) {
  ChainOptions options;
  options.Nmax = nmax;
  return build_chain_lanczos<Wide<S0, Digits>>(star, options);
}

template<typename S> auto largest(const Matrix<S> &m) { return m.cwiseAbs().maxCoeff(); }

// The largest errors of a chain against a reference: V and the hoppings relative to their largest element, the
// on-site blocks on the scale of their site.
struct Errors {
  double coupling{}, hopping{}, onsite{};
  unsigned int hopping_site{}, onsite_site{};
  [[nodiscard]] double worst() const { return std::max({coupling, hopping, onsite}); }
};

template<typename S> Errors errors(const Chain<S> &chain, const Chain<S> &reference) {
  Errors e;
  e.coupling = static_cast<double>(largest<S>(chain.V - reference.V) / largest(reference.V));
  for (unsigned int n = 0; n <= reference.Nmax; n++) {
    const auto scale = largest(reference.T[n]);
    auto local       = std::max(largest(reference.E[n]), scale);
    if (n > 0) local = std::max(local, largest(reference.T[n - 1]));
    const auto hopping = static_cast<double>(largest<S>(chain.T[n] - reference.T[n]) / scale);
    const auto onsite  = static_cast<double>(largest<S>(chain.E[n] - reference.E[n]) / local);
    if (hopping > e.hopping) {
      e.hopping      = hopping;
      e.hopping_site = n;
    }
    if (onsite > e.onsite) {
      e.onsite      = onsite;
      e.onsite_site = n;
    }
  }
  return e;
}

// A chain in a wider arithmetic, element by element through the real and imaginary parts.
template<typename To, typename From> Matrix<To> widen(const Matrix<From> &m) {
  Matrix<To> result(m.rows(), m.cols());
  for (Eigen::Index i = 0; i < m.rows(); i++)
    for (Eigen::Index j = 0; j < m.cols(); j++) result(i, j) = convert_scalar<To>(m(i, j));
  return result;
}

template<typename To, typename From> Chain<To> widen(const Chain<From> &chain) {
  Chain<To> result;
  result.channels = chain.channels;
  result.Nmax     = chain.Nmax;
  result.V        = widen<To>(chain.V);
  for (const auto &block : chain.E) result.E.push_back(widen<To>(block));
  for (const auto &block : chain.T) result.T.push_back(widen<To>(block));
  return result;
}

// The smallest eigenvalue of a hopping over its largest element: a hopping of the polar gauge has none below zero.
template<typename S0> double lowest_eigenvalue(const Matrix<S0> &hopping) {
  const Eigen::Matrix<S0, Eigen::Dynamic, Eigen::Dynamic> h = hopping;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<S0, Eigen::Dynamic, Eigen::Dynamic>> solver(h, Eigen::EigenvaluesOnly);
  return solver.eigenvalues()(0) / largest(hopping);
}

template<typename S0> void qualify(const Case &c) {
  SCOPED_TRACE(c.name);
  const auto star = star_of<S0>(c);

  // The reference, converged in its own precision before it is used.
  const auto low       = lanczos<S0, low_digits>(star, c.nmax);
  const auto high      = lanczos<S0, high_digits>(star, c.nmax);
  const auto agreement = errors(widen<Wide<S0, high_digits>>(low), high).worst();
  ASSERT_LT(agreement, reference_agreement) << "the reference is not converged";
  ASSERT_EQ(high.diagnostics.theta_rank, c.channels);
  ASSERT_FALSE(high.diagnostics.rank_drop_site.has_value());
  const auto reference = convert_chain<S0>(high);

  ChainOptions options;
  options.Nmax     = c.nmax;
  const auto chain = build_chain(star, options);
  EXPECT_EQ(chain.diagnostics.theta_rank, c.channels);
  EXPECT_FALSE(chain.diagnostics.rank_drop_site.has_value());
  for (unsigned int n = 0; n <= c.nmax; n++) {
    EXPECT_EQ(largest<S0>(chain.T[n] - chain.T[n].adjoint()), 0.0) << "site " << n;
    EXPECT_GT(lowest_eigenvalue(chain.T[n]), -1e-13) << "site " << n;
  }

  const auto e = errors(chain, reference);
  std::printf("%-24s N=%d %-7s Lambda=%-3g sites=%-2u levels=%-4zu  V %.1e  T %.1e at %-2u  E %.1e at %-2u", c.name,
              c.channels, c.complex_data ? "complex" : "real", c.lambda, c.nmax + 1, star.levels.size(), e.coupling,
              e.hopping, e.hopping_site, e.onsite, e.onsite_site);

  if (c.gap == 0.0) {
    std::printf("\n");
    EXPECT_LE(e.coupling, budget);
    EXPECT_LE(e.hopping, budget) << "at site " << e.hopping_site;
    EXPECT_LE(e.onsite, budget) << "at site " << e.onsite_site;
    return;
  }

  // Ill-conditioned: how far the exact chain moves when the star changes by one unit in the last place.
  const auto moved       = convert_chain<S0>(lanczos<S0, low_digits>(one_ulp_away(star), c.nmax));
  const auto sensitivity = errors(moved, reference).worst();
  const auto allowed     = std::max(budget, sensitivity_margin * sensitivity);
  std::printf("  one-ulp sensitivity %.1e\n", sensitivity);
  if (e.worst() > budget)
    std::printf("ADVISORY: %s exceeds the target %.0e with %.1e; the chain of this star is determined only to %.1e "
                "in double precision\n",
                c.name, budget, e.worst(), sensitivity);
  EXPECT_LE(e.worst(), allowed) << "the error is beyond what the conditioning of the star explains";
}

void qualify(const Case &c) {
  if (c.complex_data)
    qualify<std::complex<double>>(c);
  else
    qualify<double>(c);
}

// The long chains have two channels and the chains of four channels are shorter: the reference sets the cost.
const std::vector<Case> regular_cases{
  {"two_long", 2, false, 2.0, 60, 1.0, 0.0, 0.0},
  {"two_small_lambda", 2, false, 1.5, 40, 3.0, 0.0, 0.0},
  {"two_large_lambda", 2, false, 4.0, 30, 20.0, 0.0, 0.0},
  {"two_fast", 2, false, 2.0, 40, 2.0, 0.9, 0.0},
  {"two_complex", 2, true, 2.0, 40, 5.0, 0.0, 0.0},
  {"two_complex_fast", 2, true, 4.0, 25, 1.0, 0.9, 0.0},
  {"three", 3, false, 2.0, 30, 2.0, 0.0, 0.0},
  {"three_complex_fast", 3, true, 2.0, 25, 10.0, 0.9, 0.0},
  {"four", 4, false, 2.0, 25, 3.0, 0.0, 0.0},
  {"four_complex", 4, true, 1.5, 20, 1.0, 0.0, 0.0},
  {"four_fast_large_lambda", 4, false, 4.0, 20, 20.0, 0.9, 0.0},
};

const std::vector<Case> accumulating_cases{
  {"gap_two", 2, false, 2.0, 30, 1.0, 0.0, 0.3},
  {"gap_two_complex", 2, true, 2.0, 25, 3.0, 0.0, 0.1},
  {"gap_three", 3, false, 2.0, 20, 2.0, 0.0, 0.2},
};

} // namespace

TEST(MatrixChainQualification, graded_stars) { // NOLINT
  for (const auto &c : regular_cases) qualify(c);
}

TEST(MatrixChainQualification, stars_with_an_accumulation_point) { // NOLINT
  for (const auto &c : accumulating_cases) qualify(c);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
