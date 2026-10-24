#include <nrgljubljana/numerics.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <future>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_errno.h>

namespace {
using namespace NRG::Numerics;
constexpr auto pi = std::numbers::pi;
constexpr auto infinity = std::numeric_limits<double>::infinity();
constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
const Table flat{{-1.0, 0.0, 1.0}, {0.5, 0.5, 0.5}};

void check(const bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

void near(const double value, const double expected, const double atol = 1e-13, const double rtol = 1e-12) {
  if (!std::isfinite(value) || std::abs(value - expected) > atol + rtol * std::abs(expected)) {
    std::cerr.precision(17);
    std::cerr << "actual=" << value << " expected=" << expected << '\n';
    throw std::runtime_error("Numerical oracle mismatch.");
  }
}

template <typename Exception = std::exception, typename Function>
void rejects(Function function) {
  try {
    function();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error("Expected an exception.");
}

void integration() {
  near(integrate(flat), 1.0);
  near(integrate_range(flat, -0.5, 0.5), 0.5);
  near(integrate_range(flat, -2, 2), 1.0);
  near(integrate_range(flat, 2, 3), 0.0);
  near(integrate_range(flat, 0, 0), 0.0);
  rejects<std::invalid_argument>([] { integrate_range(flat, 1, -1); });
  rejects<std::invalid_argument>([] { integrate_range(flat, 0, infinity); });
  near(integrate_moment(flat, 0), 1.0);
  near(integrate_moment(flat, 1), 0.0);
  near(integrate_moment(flat, 2), 1.0 / 3.0);
  near(integrate_moment({{0.0, 1.0, 2.0}, {0.0, 0.5, 1.0}}, 1), 4.0 / 3.0);
  rejects<std::invalid_argument>([] { integrate_moment(flat, std::numeric_limits<unsigned>::max()); });
  near(integrate({{-2.0, -0.5, 1.0, 3.0}, {-0.25, 0.5, 1.25, 2.25}}), 5.0);
  const Table crossing{{-1.0, -0.25, 1.0}, {-1.0, -0.25, 1.0}};
  near(integrate(crossing), 0.0);
  near(integrate_abs(crossing), 1.0);
  near(integrate({{-1.0, 0.0, 1.0}, {-1.0, -1.0, -1.0}}), -2.0);

  // First interval is -5/8 + u + u^2 - u^3, with root u=1/2.
  // The second interval is its reflection. Integrating |sample values|
  // instead of |the cubic| does not give 9/16.
  const Table cubic{{0.0, 1.0, 2.0}, {-0.625, 0.375, -0.625}};
  near(integrate(cubic), -1.0 / 12.0);
  near(integrate_abs(cubic), 9.0 / 16.0);
  near(integrate_abs({{-1.0, 0.0, 1.0}, {0.0, 0.0, 0.0}}), 0.0);

  for (const auto temperature : {0.0, 1e-200, 0.01, 0.5, 1e200})
    near(integrate_fermi(flat, temperature), 0.5);
  near(integrate_fermi(crossing, 0.0), -0.5);
  near(integrate_fermi(crossing, 0.01), -0.5 + pi * pi * 0.0001 / 6.0);
  const Table positive{{1.0, 2.0, 3.0}, {3.0, 3.0, 3.0}};
  near(integrate_fermi(positive, 0.0), 0.0);
  near(integrate_fermi(positive, 0.7), 2.1 * (std::log1p(std::exp(-1.0 / 0.7)) - std::log1p(std::exp(-3.0 / 0.7))));
  near(integrate_fermi({{-3.0, -2.0, -1.0}, {-3.0, -3.0, -3.0}}, 0.0), -6.0);
  near(integrate_fermi(flat, 0.5, 0.0, 1e-9), 0.5);
}

void resampling() {
  const Table line{{-1.0, 0.0, 1.0}, {-1.0, 1.0, 3.0}};
  const std::vector<double> grid{-1.0, -0.7, 0.2, 0.8, 1.0};
  const auto sampled = resample(line, grid);
  check(sampled.x == grid, "Resampling changed grid labels.");
  for (std::size_t i = 0; i < grid.size(); ++i) near(sampled.y[i], 2.0 * grid[i] + 1.0);
  near(resample({{0.0, 1.0, 2.0}, {0.0, 1.0, 0.0}}, {0.5}).y[0], 0.625);
  const std::vector<double> near_endpoints{std::nextafter(-1.0, -infinity), std::nextafter(1.0, infinity)};
  const auto clamped = resample(line, near_endpoints);
  check(clamped.x == near_endpoints, "Roundoff clamping changed output labels.");
  near(clamped.y.front(), -1.0);
  near(clamped.y.back(), 3.0);
  rejects<std::domain_error>([&] { resample(line, {-1.0 - 1e-12}); });
  rejects<std::domain_error>([&] { resample(line, {1.0 + 1e-12}); });
  // No absolute unit-scale tolerance on small-energy tables.
  rejects<std::domain_error>([] { resample({{-1e-20, 0.0, 1e-20}, {1.0, 1.0, 1.0}}, {1.1e-20}); });
}

void kk() {
  const std::vector<double> grid{-2.0, -1.0, -0.1, 0.1, 1.0, 2.0};
  const auto even = kramers_kronig({grid, std::vector<double>(grid.size(), -1.0)});
  const auto odd = kramers_kronig({grid, grid});
  check(even.x == grid && odd.x == grid, "KK changed the input mesh.");
  for (std::size_t i = 0; i < grid.size(); ++i) {
    const auto x = grid[i];
    const auto log_ratio = std::abs(x) == 2.0 ? 0.0 : std::log((2.0 - x) / (2.0 + x));
    near(even.y[i], -log_ratio / pi);
    near(odd.y[i], (4.0 + x * log_ratio) / pi);
    near(even.y[i], -even.y[grid.size() - 1 - i]);
    near(odd.y[i], odd.y[grid.size() - 1 - i]);
  }
  check(even.y[4] > 0.0, "Im-to-Re KK sign is reversed.");
  near(even.y.front(), 0.0);
  near(odd.y.front(), 4.0 / pi);
  near(odd.y.back(), 4.0 / pi);

  const auto nonlinear = kramers_kronig({grid, {0.0, -0.2, -1.0, -1.0, -0.2, 0.0}});
  for (std::size_t i = 0; i < grid.size(); ++i) near(nonlinear.y[i], -nonlinear.y[grid.size() - 1 - i]);
  const auto asymmetric = kramers_kronig({grid, {-3.0, -2.0, -1.1, -0.9, 0.0, 1.0}});
  for (std::size_t i = 0; i < grid.size(); ++i) near(asymmetric.y[i], even.y[i] + odd.y[i]);
  const auto accepted = kramers_kronig({{-2.0, -1.0, 1.0, 2.0 + 1e-10}, {1.0, 1.0, 1.0, 1.0}});
  near(accepted.y.front(), 0.0);
  near(accepted.y.back(), 0.0);
}

void lattice_oracles() {
  const Table real_sigma{{-2.0, -0.5, 0.0, 0.75, 2.0}, {0.3, -0.1, 0.2, -0.2, 0.1}};
  const Table imag_sigma{real_sigma.x, {-0.2, 0.0, 0.25, -1e-8, -0.1}};
  for (unsigned moment = 0; moment <= 2; ++moment) {
    const auto result = lattice(flat, real_sigma, imag_sigma, 0.13, 0.05, moment);
    check(result.real.x == real_sigma.x && result.spectral.x == real_sigma.x, "Lattice output grid mismatch.");
    for (std::size_t i = 0; i < real_sigma.x.size(); ++i) {
      const std::complex<double> z{real_sigma.x[i] + 0.13 - real_sigma.y[i], std::max(-imag_sigma.y[i], 0.05)};
      auto expected = 0.5 * (std::log(z + 1.0) - std::log(z - 1.0));
      // Recurrence is an independent oracle only here, at moderate |z|.
      if (moment >= 1) expected = z * expected - 1.0;
      if (moment >= 2) expected *= z;
      near(result.real.y[i], expected.real());
      near(result.spectral.y[i], -expected.imag() / pi);
    }
  }
  const Table sigma{{-0.4, 0.0, 0.4}, {0.0, 0.0, 0.0}};
  const auto first = lattice(flat, sigma, sigma, 0.0, 0.02, 1);
  check(first.spectral.y.front() < 0.0 && first.spectral.y.back() > 0.0, "Odd moment lost signed spectral weight.");
  for (unsigned moment = 0; moment <= 1; ++moment) {
    const auto tilted = lattice({{-1.0, 0.0, 1.0}, {0.4, 0.5, 0.6}}, sigma, sigma, 0.0, 0.02, moment);
    for (std::size_t i = 0; i < sigma.x.size(); ++i) {
      const std::complex<double> z{sigma.x[i], 0.02};
      const auto h0 = 0.5 * (std::log(z + 1.0) - std::log(z - 1.0));
      const auto h1 = z * h0 - 1.0;
      const auto expected = moment == 0 ? h0 + 0.2 * h1 : h1 + 0.2 * z * h1;
      near(tilted.real.y[i], expected.real());
      near(tilted.spectral.y[i], -expected.imag() / pi);
    }
  }
  const auto asymmetric = lattice({{0.0, 1.0, 2.0}, {2.0, 2.0, 2.0}}, {{0.0}, {0.0}}, {{0.0}, {-0.2}}, 0.5, 0.01);
  const std::complex<double> z{0.5, 0.2};
  const auto expected = 2.0 * (std::log(z) - std::log(z - 2.0));
  near(asymmetric.real.y[0], expected.real());
  near(asymmetric.spectral.y[0], -expected.imag() / pi);

  // z*G_0 - mass rounds to zero here. Direct polynomial weighting must not.
  const auto far = lattice(flat, {{0.0}, {0.0}}, {{0.0}, {-0.25}}, 1e12, 0.01, 1);
  const std::complex<double> far_z{1e12, 0.25};
  const auto far_expected = 1.0 / (3.0 * far_z * far_z);
  near(far.real.y[0], far_expected.real(), 0.0, 2e-14);
  near(far.spectral.y[0], -far_expected.imag() / pi, 0.0, 2e-14);

  const auto narrow = lattice(flat, {{0.0}, {0.0}}, {{0.0}, {0.0}}, 0.0, std::sqrt(std::numeric_limits<double>::min()));
  near(narrow.real.y[0], 0.0);
  near(narrow.spectral.y[0], 0.5);
  const auto zero = lattice({flat.x, {0.0, 0.0, 0.0}}, sigma, sigma, 0.0, 0.1);
  for (const auto value : zero.real.y) near(value, 0.0);
  for (const auto value : zero.spectral.y) near(value, 0.0);
}

void validation() {
  const std::vector<Table> bad_tables{
    {}, {{0.0}, {1.0}}, {{0.0, 1.0}, {1.0, 1.0}}, {{0.0, 1.0, 2.0}, {1.0, 1.0}},
    {{0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}}, {{0.0, 2.0, 1.0}, {1.0, 1.0, 1.0}},
    {{0.0, nan, 1.0}, {1.0, 1.0, 1.0}}, {{0.0, 1.0, infinity}, {1.0, 1.0, 1.0}},
    {{0.0, 1.0, 2.0}, {1.0, nan, 1.0}}, {{0.0, 1.0, 2.0}, {1.0, 1.0, infinity}},
    {{-std::numeric_limits<double>::max(), 0.0, std::numeric_limits<double>::max()}, {1.0, 1.0, 1.0}}
  };
  for (const auto &table : bad_tables) {
    rejects([&] { integrate(table); });
    rejects([&] { integrate_abs(table); });
    rejects([&] { integrate_fermi(table, 0.1); });
    rejects([&] { resample(table, {0.0}); });
    rejects([&] { kramers_kronig(table); });
    rejects([&] { lattice(table, flat, flat, 0.0, 0.1); });
  }
  for (const auto &grid : std::vector<std::vector<double>>{{}, {nan}, {infinity}, {0.0, 0.0}, {0.5, -0.5}})
    rejects([&] { resample(flat, grid); });
  for (const auto temperature : {-1.0, infinity, nan}) rejects([&] { integrate_fermi(flat, temperature); });
  for (const auto &tolerances : std::vector<std::pair<double, double>>{
         {-1.0, 1e-9}, {1e-10, -1.0}, {nan, 1e-9}, {1e-10, nan}, {infinity, 1e-9}, {1e-10, infinity}, {0.0, 0.0}, {0.0, 1e-17}})
    rejects([&] { integrate_fermi(flat, 0.0, tolerances.first, tolerances.second); });
  rejects([] { kramers_kronig(flat); }); // Odd point count, even though symmetric.
  rejects([] { kramers_kronig({{-2.0, -1.0, 0.5, 2.0}, {1.0, 1.0, 1.0, 1.0}}); });
  rejects([] { kramers_kronig({{-2.0, -1.0, 1.0, 3.0}, {1.0, 1.0, 1.0, 1.0}}); });

  for (const auto mu : {infinity, nan}) rejects([&] { lattice(flat, flat, flat, mu, 0.1); });
  for (const auto clip : {-1.0, 0.0, 1e-200, infinity, nan}) rejects([&] { lattice(flat, flat, flat, 0.0, clip); });
  rejects([] { lattice({{-1.0, 0.0, 1.0}, {0.5, -0.1, 0.5}}, flat, flat, 0.0, 0.1); });
  rejects([] { lattice(flat, {}, flat, 0.0, 0.1); });
  rejects([] { lattice(flat, {{0.0}, {nan}}, {{0.0}, {0.0}}, 0.0, 0.1); });
  rejects([] { lattice(flat, flat, {{-1.0, 0.0, 1.0}, {0.0}}, 0.0, 0.1); });
  rejects([] { lattice(flat, flat, {{-1.0, 0.0, 1.0 + 1e-14}, {0.0, 0.0, 0.0}}, 0.0, 0.1); });
  rejects([] { lattice(flat, flat, flat, 0.0, 0.1, std::numeric_limits<unsigned>::max()); });
  rejects([] { lattice(flat, {{1e308}, {0.0}}, {{1e308}, {0.0}}, 1e308, 0.1); });
}

void custom_gsl_handler(const char *, const char *, int, int) {}

void gsl_handler() {
  struct Restore {
    gsl_error_handler_t *previous = gsl_set_error_handler(&custom_gsl_handler);
    ~Restore() { gsl_set_error_handler(previous); }
  } restore;
  gsl_set_error_handler_off();
  const auto off_handler = gsl_set_error_handler_off();
  check(off_handler != nullptr, "GSL's disabled handler must differ from its default handler.");
  for (const auto handler : std::vector<gsl_error_handler_t *>{&custom_gsl_handler, nullptr, off_handler}) {
    gsl_set_error_handler(handler);
    const auto restored = [handler] {
      check(gsl_set_error_handler(handler) == handler, "GSL handler was not restored.");
    };
    integration();
    restored();
    resampling();
    restored();
    kk();
    restored();
    lattice_oracles();
    restored();
    rejects([] { integrate_fermi(flat, 0.1, 1e-300, 0.0); });
    restored();
    rejects([] { integrate({{0.0, 1.0, 2.0}, {-1e308, 1e308, 0.0}}); });
    restored();

    std::vector<std::future<void>> workers;
    for (int i = 0; i < 4; ++i)
      workers.push_back(std::async(std::launch::async, [] {
        for (int iteration = 0; iteration < 4; ++iteration) near(integrate_fermi(flat, 0.1), 0.5);
      }));
    for (auto &worker : workers) worker.get();
    restored();
  }
}

} // namespace

int main(const int argc, char **argv) {
  try {
    check(argc == 2, "Expected one test name.");
    const std::string name{argv[1]};
    if (name == "integration") integration();
    else if (name == "resample") resampling();
    else if (name == "kk") kk();
    else if (name == "lattice") lattice_oracles();
    else if (name == "validation") validation();
    else if (name == "gsl_handler") gsl_handler();
    else throw std::runtime_error("Unknown test: " + name);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
