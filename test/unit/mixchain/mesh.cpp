#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include <mixchain/mesh.hpp>

using namespace NRG::MixChain;

namespace {

const NRG::Tools::LambdaCache lambda{2.0};

// A weight function tabulated on 'grid'. The tables start at omega=0, as the input branches do.
Vec weight_table(const std::vector<double> &grid, const std::vector<double> &values) {
  Vec table;
  for (std::size_t k = 0; k < grid.size(); k++) table.emplace_back(grid[k], values[k]);
  return table;
}

Vec flat_weight() { return weight_table({0.0, 0.5, 1.0}, {1.0, 1.0, 1.0}); }

Mesh adaptive_mesh(const Vec &weight, const bool hardgap = false, const double boundary = 0.0) {
  return Mesh(lambda, hardgap, boundary, weight, NRG::Tools::InterpolationMethod::linear);
}

template<typename S> GammaBranch<S> diagonal_branch(const double first, const double second) {
  GammaBranch<S> branch;
  branch.omega.push_back(0.5);
  Matrix<S> m = Matrix<S>::Zero(2, 2);
  m(0, 0)     = first;
  m(1, 1)     = second;
  branch.gamma.push_back(m);
  return branch;
}

} // namespace

TEST(MixChainMesh, fixed_mesh_follows_the_logarithmic_formula) { // NOLINT
  Mesh mesh(lambda, false, 0.0);
  EXPECT_DOUBLE_EQ(mesh.eps(1.0), 1.0);
  EXPECT_DOUBLE_EQ(mesh.eps(2.0), 1.0);
  EXPECT_DOUBLE_EQ(mesh.eps(3.0), 0.5);
  EXPECT_DOUBLE_EQ(mesh.eps(4.5), std::pow(2.0, -2.5));
  EXPECT_FALSE(mesh.adaptive());
}

TEST(MixChainMesh, fixed_mesh_applies_the_hardgap_rescaling) { // NOLINT
  Mesh mesh(lambda, true, 0.1);
  EXPECT_DOUBLE_EQ(mesh.eps(2.0), 1.0);
  EXPECT_DOUBLE_EQ(mesh.eps(3.0), 0.9 * 0.5 + 0.1);
  EXPECT_NEAR(mesh.eps(40.0), 0.1, 1e-10); // the accumulation point is the gap edge
}

TEST(MixChainMesh, rejects_an_out_of_range_boundary) { // NOLINT
  EXPECT_THROW(Mesh(lambda, true, 1.0), std::invalid_argument);
  EXPECT_THROW(Mesh(lambda, true, -0.1), std::invalid_argument);
  EXPECT_THROW(adaptive_mesh(flat_weight(), true, 1.5), std::invalid_argument);

  // Without hardgap the value is unused, so it is not validated.
  Mesh mesh(lambda, false, 1.5);
  EXPECT_DOUBLE_EQ(mesh.eps(3.0), 0.5);
}

TEST(MixChainMesh, adaptive_mesh_on_a_flat_weight_reproduces_the_fixed_mesh) { // NOLINT
  // A constant weight gives W(omega)=omega, so the inverse is the identity and eps(x)=Lambda^(2-x).
  Mesh fixed(lambda, false, 0.0);
  auto adaptive = adaptive_mesh(flat_weight());
  EXPECT_TRUE(adaptive.adaptive());
  for (const double x : {2.0, 2.5, 3.0, 5.0, 9.0, 17.0}) {
    // inverse() bisects to about one ulp rather than returning the exact value.
    EXPECT_NEAR(adaptive.eps(x), fixed.eps(x), 1e-15 * fixed.eps(x));
  }
}

TEST(MixChainMesh, adaptive_mesh_agrees_for_both_interpolation_methods) { // NOLINT
  // The cumulative of a flat weight is inverted in closed form for the linear interpolant and by a bisection inside
  // the bracketing interval for steffen; on this weight both must give the fixed mesh.
  Mesh fixed(lambda, false, 0.0);
  Mesh steffen(lambda, false, 0.0, weight_table({0.0, 0.25, 0.5, 0.75, 1.0}, {1.0, 1.0, 1.0, 1.0, 1.0}),
               NRG::Tools::InterpolationMethod::steffen);
  for (const double x : {2.0, 2.5, 3.0, 5.0, 9.0, 17.0})
    EXPECT_NEAR(steffen.eps(x), fixed.eps(x), 1e-14 * fixed.eps(x));
}

TEST(MixChainMesh, adaptive_mesh_composes_with_the_hardgap_rescaling) { // NOLINT
  Mesh fixed(lambda, true, 0.1);
  auto adaptive = adaptive_mesh(flat_weight(), true, 0.1);
  for (const double x : {2.0, 3.0, 6.0, 12.0}) EXPECT_NEAR(adaptive.eps(x), fixed.eps(x), 1e-15);
}

TEST(MixChainMesh, both_meshes_are_decreasing_and_bounded_by_one) { // NOLINT
  Mesh fixed(lambda, false, 0.0);
  auto adaptive = adaptive_mesh(weight_table({0.0, 0.1, 0.3, 1.0}, {2.0, 1.0, 0.5, 0.5}));
  double previous_fixed    = fixed.eps(2.0);
  double previous_adaptive = adaptive.eps(2.0);
  EXPECT_DOUBLE_EQ(previous_fixed, 1.0);
  EXPECT_NEAR(previous_adaptive, 1.0, 1e-15);
  for (double x = 2.25; x <= 30.0; x += 0.25) {
    const auto value_fixed    = fixed.eps(x);
    const auto value_adaptive = adaptive.eps(x);
    EXPECT_LT(value_fixed, previous_fixed);
    EXPECT_LT(value_adaptive, previous_adaptive);
    EXPECT_LE(value_adaptive, 1.0);
    previous_fixed    = value_fixed;
    previous_adaptive = value_adaptive;
  }
}

TEST(MixChainMesh, adaptive_mesh_follows_the_weight) { // NOLINT
  // Almost all of the weight sits below omega=0.1, so the intervals are pulled towards small frequencies.
  Mesh fixed(lambda, false, 0.0);
  auto adaptive = adaptive_mesh(weight_table({0.0, 0.1, 1.0}, {1.0, 1.0, 0.01}));
  for (const double x : {3.0, 4.0, 6.0}) EXPECT_LT(adaptive.eps(x), fixed.eps(x));
}

TEST(MixChainMesh, adaptive_mesh_accumulates_at_the_edge_of_a_gap) { // NOLINT
  // The weight vanishes identically below 0.2: W is flat there, so the inverse cannot return anything below the gap
  // edge. Above it the weight rises linearly from 0 to 1, so W(omega) = (omega-0.2)^2/0.64 and the mesh is
  // eps(x) = 0.2 + 0.8 sqrt(Lambda^(2-x)).
  auto mesh              = adaptive_mesh(weight_table({0.0, 0.2, 1.0}, {0.0, 0.0, 1.0}));
  const auto expected_eps = [](const double x) { return 0.2 + 0.8 * std::sqrt(std::pow(2.0, 2.0 - x)); };
  EXPECT_GE(mesh.eps(10.0), 0.2);
  EXPECT_GE(mesh.eps(30.0), 0.2);
  EXPECT_NEAR(mesh.eps(10.0), expected_eps(10.0), 1e-12);
  EXPECT_NEAR(mesh.eps(30.0), expected_eps(30.0), 1e-12);
}

TEST(MixChainMesh, knows_its_accumulation_point) { // NOLINT
  Mesh fixed(lambda, false, 0.0);
  EXPECT_EQ(fixed.accumulation_point(), 0.0);

  Mesh gapped(lambda, true, 0.1);
  EXPECT_DOUBLE_EQ(gapped.accumulation_point(), 0.1); // the boundary

  // The adaptive mesh finds the edge of a region where the weight vanishes, and zero when there is none.
  auto adaptive_gap = adaptive_mesh(weight_table({0.0, 0.2, 1.0}, {0.0, 0.0, 1.0}));
  EXPECT_DOUBLE_EQ(adaptive_gap.accumulation_point(), 0.2);
  auto adaptive_flat = adaptive_mesh(flat_weight());
  EXPECT_EQ(adaptive_flat.accumulation_point(), 0.0);

  // Both limits are what eps(x) actually approaches.
  EXPECT_NEAR(gapped.eps(60.0), gapped.accumulation_point(), 1e-15);
  EXPECT_NEAR(adaptive_gap.eps(60.0), adaptive_gap.accumulation_point(), 1e-8);
}

TEST(MixChainMesh, weight_table_offers_the_frobenius_norm_and_the_trace) { // NOLINT
  const auto branch = diagonal_branch<double>(0.3, 0.4);
  const auto frobenius = mesh_weight_table(branch, MeshWeight::frobenius);
  const auto trace     = mesh_weight_table(branch, MeshWeight::trace);

  ASSERT_EQ(frobenius.size(), 1U);
  EXPECT_DOUBLE_EQ(frobenius[0].first, 0.5);
  EXPECT_DOUBLE_EQ(frobenius[0].second, 0.5); // sqrt(0.3^2+0.4^2)
  EXPECT_DOUBLE_EQ(trace[0].second, 0.7);

  // The same holds in complex arithmetic, where the off-diagonal elements contribute to the norm only.
  const auto complex_branch = diagonal_branch<std::complex<double>>(0.3, 0.4);
  EXPECT_DOUBLE_EQ(mesh_weight_table(complex_branch, MeshWeight::trace)[0].second, 0.7);
}

TEST(MixChainMesh, weight_table_rejects_a_negative_trace) { // NOLINT
  const auto branch = diagonal_branch<double>(0.3, -0.5);
  EXPECT_THROW(mesh_weight_table(branch, MeshWeight::trace), std::runtime_error);
  // A norm cannot be negative, so the same input is accepted there.
  EXPECT_NO_THROW(mesh_weight_table(branch, MeshWeight::frobenius));
}

TEST(MixChainMesh, parses_the_weight_option) { // NOLINT
  EXPECT_EQ(mesh_weight_from_string("frobenius"), MeshWeight::frobenius);
  EXPECT_EQ(mesh_weight_from_string("trace"), MeshWeight::trace);
  EXPECT_EQ(mesh_weight_name(MeshWeight::trace), "trace");
  EXPECT_THROW(mesh_weight_from_string("norm"), std::invalid_argument);
  EXPECT_THROW(mesh_weight_from_string(""), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
