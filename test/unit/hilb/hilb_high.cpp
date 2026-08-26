#include <gtest/gtest.h>

#include <complex>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <vector>

#include <hilb/hilb.hpp>

namespace {

constexpr double B = 1.0;

auto flat_band_h0(const std::complex<double> z) { return (std::log(z + B) - std::log(z - B)) / (2.0 * B); }

void expect_complex_near(const std::complex<double> actual, const std::complex<double> expected, const double tolerance) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

} // namespace

TEST(Hilb, empty_dos_throws) { // NOLINT
  const auto filename = "hilb_empty_dos.dat";
  std::ofstream file(filename);

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_empty_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_THROW(NRG::Hilb::Hilb(5, argv), std::runtime_error);

  std::remove(filename);
}

TEST(Hilb, unsorted_dos_is_sorted_before_use) { // NOLINT
  const auto filename = "hilb_unsorted_dos.dat";
  {
    std::ofstream file(filename);
    file << "1 1\n-1 1\n0 1\n";
  }

  char arg0[] = "hilb";
  char arg1[] = "-d";
  char arg2[] = "hilb_unsorted_dos.dat";
  char arg3[] = "0.1";
  char arg4[] = "0.2";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4};

  optind = 1;
  EXPECT_NO_THROW(NRG::Hilb::Hilb(5, argv));

  std::remove(filename);
}

TEST(Hilb, interpolator_move_assignment_releases_old_state) { // NOLINT
  auto a = NRG::Hilb::interpolator(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{1.0, 2.0, 3.0});
  auto b = NRG::Hilb::interpolator(std::vector<double>{-2.0, 0.0, 2.0}, std::vector<double>{4.0, 5.0, 6.0});
  a = std::move(b);
  EXPECT_DOUBLE_EQ(a(0.0), 5.0);
}

TEST(Hilb, energy_power_direct_quadrature) { // NOLINT
  auto rho = [](const double) { return 0.5; };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{0.2, 0.1};
  const auto expected = z * flat_band_h0(z) - 1.0;

  const auto result = NRG::Hilb::hilbert_transform(rho, zero, B, z, 1e-3, 1);

  expect_complex_near(result, expected, 1e-10);
}

TEST(Hilb, energy_power_singularity_subtraction) { // NOLINT
  auto rhor = [](const double) { return 0.5; };
  auto rhoi = [](const double) { return 0.125; };
  const std::complex<double> z{0.2, 1e-6};
  const auto real_density_result = z * z * flat_band_h0(z) - z;
  const auto expected = std::complex<double>{1.0, 0.25} * real_density_result;

  const auto result = NRG::Hilb::hilbert_transform(rhor, rhoi, B, z, 1e-3, 2);

  expect_complex_near(result, expected, 1e-9);
}

TEST(Hilb, tabulated_density_is_weighted_after_interpolation) { // NOLINT
  const std::vector<double> energies{-1.0, -0.5, 0.0, 0.5, 1.0};
  const std::vector<double> density(energies.size(), 0.5);
  const std::vector<double> zero(energies.size(), 0.0);
  const std::complex<double> z{-0.3, 0.1};
  const auto expected = z * z * flat_band_h0(z) - z;

  const auto result = NRG::Hilb::hilbert_transform(energies, density, zero, z, 1e-3, 2);

  expect_complex_near(result, expected, 1e-10);
}

TEST(Hilb, zero_energy_power_preserves_transform) { // NOLINT
  auto rho = [](const double energy) { return 0.5 + 0.1 * energy; };
  auto zero = [](const double) { return 0.0; };
  const std::complex<double> z{0.2, 0.1};

  const auto implicit_zero = NRG::Hilb::hilbert_transform(rho, zero, B, z);
  const auto explicit_zero = NRG::Hilb::hilbert_transform(rho, zero, B, z, 1e-3, 0);

  EXPECT_EQ(implicit_zero, explicit_zero);
}

TEST(Hilb, energy_power_parser_requires_nonnegative_integer) { // NOLINT
  auto rho = [](double) { return 1.0; };
  auto zero = [](double) { return 0.0; };
  EXPECT_EQ(NRG::Hilb::parse_energy_power("0"), 0);
  EXPECT_EQ(NRG::Hilb::parse_energy_power("12"), 12);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("-1"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("1.5"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("2junk"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::parse_energy_power("999999999999999999999999"), std::runtime_error);
  EXPECT_THROW(NRG::Hilb::hilbert_transform(rho, zero, B, std::complex<double>{0.0, 0.1}, 1e-3, -1), std::invalid_argument);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
