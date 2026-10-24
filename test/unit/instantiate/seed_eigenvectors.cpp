#include <bit>
#include <cstdint>
#include <gtest/gtest.h>
#include <diag.hpp>
#include <instantiate/seed_eigenvectors.hpp>

using NRG::Matrix_traits;
using NRG::diagonalise_dsyev;
using NRG::Spawn::canonicalize_seed_eigenvectors;

TEST(SeedEigenvectors, SignChoiceTiesAndZeros) {
  Matrix_traits<double> vectors(5, 3);
  vectors << -0.0, -0.6, 0.8, 0.0, -1.0, 0.0, -0.5, 0.5, -0.0, 0.0, 0.0, -1.0, 1e-300, -0.8, 0.6;
  Matrix_traits<double> expected(5, 3);
  expected << 0.0, -0.6, 0.8, 0.0, 1.0, 0.0, 0.5, -0.5, 0.0, 0.0, 0.0, 1.0, -1e-300, 0.8, -0.6;
  for (unsigned mask = 0; mask < 32; ++mask) {
    auto signed_vectors = vectors;
    for (Eigen::Index row = 0; row < vectors.rows(); ++row)
      if (mask & (1U << row)) signed_vectors.row(row) *= -1.0;
    for (int repeat = 0; repeat < 2; ++repeat) {
      canonicalize_seed_eigenvectors(signed_vectors);
      for (Eigen::Index row = 0; row < vectors.rows(); ++row)
        for (Eigen::Index col = 0; col < vectors.cols(); ++col)
          EXPECT_EQ(std::bit_cast<uint64_t>(signed_vectors(row, col)), std::bit_cast<uint64_t>(expected(row, col)));
    }
  }
}

TEST(SeedEigenvectors, PreservesEigensystemAndOperatorTransforms) {
  Matrix_traits<double> hamiltonian(3, 3);
  hamiltonian << 0.1, 0.4, 0.0, 0.4, -0.13, 0.4, 0.0, 0.4, 0.04;
  auto input = hamiltonian;
  auto raw = diagonalise_dsyev(input);
  canonicalize_seed_eigenvectors(raw.vec);
  for (Eigen::Index row = 0; row < 3; ++row) {
    EXPECT_NEAR(raw.vec.row(row).squaredNorm(), 1.0, 1e-14);
    for (Eigen::Index col = 0; col < 3; ++col) {
      double product = 0.0;
      for (Eigen::Index i = 0; i < 3; ++i) product += hamiltonian(col, i) * raw.vec(row, i);
      EXPECT_NEAR(product, raw.val[row] * raw.vec(row, col), 1e-14);
    }
  }
  Matrix_traits<double> op(3, 3);
  op << 0.0, 1.0, -0.5, 0.25, 2.0, 0.0, 1.0, -0.25, 0.5;
  const Matrix_traits<double> transformed = raw.vec * op * raw.vec.transpose();
  for (unsigned mask = 0; mask < 8; ++mask) {
    auto signed_vectors = raw.vec;
    for (Eigen::Index row = 0; row < 3; ++row)
      if (mask & (1U << row)) signed_vectors.row(row) *= -1.0;
    canonicalize_seed_eigenvectors(signed_vectors);
    const Matrix_traits<double> result = signed_vectors * op * signed_vectors.transpose();
    for (Eigen::Index row = 0; row < 3; ++row)
      for (Eigen::Index col = 0; col < 3; ++col) EXPECT_EQ(result(row, col), transformed(row, col));
  }
}
