// Channel-mixing discretization for NRG
// ** Star discretization of Gamma(omega)

#ifndef _mixchain_star_hpp_
#define _mixchain_star_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../common/cumulative_weight.hpp"
#include "../common/gsl_config.hpp"
#include "../common/lambda.hpp"
#include "../common/representative_energy.hpp"
#include "../common/tabulated_density.hpp"
#include "branches.hpp"
#include "gamma_interp.hpp"
#include "linint.hpp"
#include "load.hpp"
#include "mesh.hpp"
#include "types.hpp"

namespace NRG::MixChain {

struct StarOptions {
  NRG::Tools::LambdaCache Lambda;
  double z{1.0};
  unsigned int mMAX{0};
  // Recorded in the star, and from there in the star file. The rescaling itself is applied to Gamma when it is read,
  // so nothing here depends on it; it travels with the star so that the chain stage can check that it was built for
  // the same band.
  double bandrescale{1.0};
  bool adapt{false};
  bool hardgap{false};
  double boundary{0.0};
  MeshWeight mesh_weight{MeshWeight::frobenius};
  NRG::Tools::InterpolationMethod interpolation{NRG::Tools::InterpolationMethod::linear};
  BranchOptions branches;
  NRG::Tools::CquadOptions cquad;
  double allowed_error{1e-10};
  // Branches whose representative energies coincide within this relative tolerance take their eigenvectors from a
  // single diagonalization, so that exactly degenerate branches give an orthonormal set of coupling vectors.
  double coincidence_tolerance{1e-10};
};

// One bath level of the star: energy E and coupling vector v = sqrt(w) u, in the normalization of the input.
template<typename S> struct StarLevel {
  int m{};         // interval index
  Sign sign{};     // the frequency branch it came from
  int branch{};    // the eigenvalue branch it came from
  double energy{}; // signed
  Vector<S> coupling;
};

// How far the mesh reaches compared with the input. Where the mesh goes below the innermost tabulated frequency the
// density is the constant continuation of the input, which is exact for a flat band and an approximation for
// anything with structure at low frequency.
struct BranchCoverage {
  double lowest_mesh{};     // the smallest |omega| the mesh reaches
  double innermost_input{}; // the smallest |omega| tabulated in the input
  // Intervals of the mesh that contain no node of the input tabulation at all. There the discretization follows the
  // interpolant between two tabulated points rather than the data, whatever the density does in between. This
  // happens wherever the mesh resolves more finely than the input: at the bottom of the band, and around an
  // accumulation point set by hardgap or found by the adaptive mesh.
  int unresolved_intervals{};
  double unresolved_from{}; // the upper edge of the outermost such interval
  double unresolved_to{};   // and its lower edge
  // Levels lost to double precision near the accumulation point of the mesh. Near an accumulation point away from
  // zero, set by hardgap or found by the adaptive mesh at a gap edge, the distance to it falls below the spacing of
  // doubles there: the two bounds of an interval become the same number, the interval has no width, and its levels
  // carry no weight. They are inert, the star is in effect truncated there, and the last levels before it lose
  // relative accuracy in their weights, which are differences of nearly equal numbers. Near zero the bounds keep their
  // relative precision and this does not happen.
  //
  // The count is of levels in intervals whose bounds coincide, which is what makes them inert. It is not a comparison
  // of the energies with the accumulation point: those land on it or on a neighbouring double, depending on whether
  // the density vanishes below it.
  double accumulation_point{};
  int collapsed_levels{};
  [[nodiscard]] auto continued() const { return innermost_input > 0.0 && lowest_mesh < innermost_input; }
};

struct StarDiagnostics {
  double max_interval_deviation{}; // max over intervals of ||sum_a w_a u_a u_a^dag - int Gamma|| / ||int Gamma||
  double max_interval_omega{};     // the upper edge of the interval where that occurred
  double max_cquad_error{};        // the largest CQUAD error estimate of the integral method
  std::vector<double> crossings_pos, crossings_neg;
  BranchCoverage coverage_pos, coverage_neg;
};

template<typename S> struct Star {
  int channels{};
  unsigned int mMAX{};
  double z{};
  double Lambda{};
  double bandrescale{1.0};
  std::vector<StarLevel<S>> levels;
  Matrix<S> theta;       // sum_k v_k v_k^dagger over the star that was built
  Matrix<S> theta_exact; // the integral of Gamma over the range the mesh covers
  StarDiagnostics diagnostics;
};

namespace detail {

inline auto make_cquad_workspace(const NRG::Tools::CquadOptions &options) {
  const std::size_t limit = options.workspace_limit.value_or(1000);
  NRG::Tools::validate_cquad_workspace_limit(limit);
  std::unique_ptr<gsl_integration_cquad_workspace, NRG::Tools::GslWorkspaceDeleter> workspace(
    gsl_integration_cquad_workspace_alloc(limit));
  if (!workspace) throw std::runtime_error("Failed to allocate integration workspace.");
  return workspace;
}

// The node of the branch that brackets a representative energy, used as the reference for labelling the
// eigenvectors there.
inline auto bracketing_node(const std::vector<double> &omega, const double energy) {
  const auto upper = std::upper_bound(omega.begin(), omega.end(), energy);
  if (upper == omega.begin()) return std::size_t{0};
  return static_cast<std::size_t>(std::distance(omega.begin(), upper) - 1);
}

// Discretize one frequency branch and append its levels to the star.
//
// Everything here is a local variable on purpose: CumulativeWeight keeps a pointer to its density and
// IntegralRepresentativeEnergy keeps pointers to the mesh and to its cumulative weight, so each vector is reserved
// to its final size before the next one is built on top of it, and the mesh is never moved once it is wired up.
template<typename S>
void discretize_sign(const GammaBranch<S> &branch, const Sign sign, const StarOptions &options, Star<S> &star) {
  const auto decomposition = decompose_branch(branch, options.branches);
  const auto channels      = static_cast<std::size_t>(decomposition.channels);
  const auto dimension     = static_cast<Eigen::Index>(channels);

  Mesh mesh = options.adapt ? Mesh(options.Lambda, options.hardgap, options.boundary,
                                   mesh_weight_table(branch, options.mesh_weight), options.interpolation)
                            : Mesh(options.Lambda, options.hardgap, options.boundary);
  GammaInterpolation<S> interpolation(branch, options.interpolation);

  std::vector<NRG::Tools::TabulatedDensity> densities;
  densities.reserve(channels);
  for (std::size_t a = 0; a < channels; a++) densities.emplace_back(decomposition.density[a], options.interpolation);

  // A branch whose density vanishes over the whole band carries no weight anywhere: Gamma is rank deficient, which
  // is a legitimate input. Its levels are kept, with vanishing coupling, but it has no cumulative weight.
  std::vector<bool> empty(channels, false);
  std::vector<NRG::Tools::CumulativeWeight> cumulatives;
  cumulatives.reserve(channels);
  for (std::size_t a = 0; a < channels; a++) {
    empty[a] = !(densities[a].integral(0.0, 1.0) > 0.0);
    cumulatives.emplace_back();
    if (!empty[a]) cumulatives.back() = NRG::Tools::CumulativeWeight(densities[a], decomposition.density[a]);
  }

  auto workspace   = detail::make_cquad_workspace(options.cquad);
  double max_error = 0.0;
  using Representative = NRG::Tools::IntegralRepresentativeEnergy<Mesh>;
  std::vector<Representative> representatives;
  representatives.reserve(channels);
  for (std::size_t a = 0; a < channels; a++) {
    representatives.emplace_back();
    if (!empty[a])
      representatives.back() = Representative(mesh, cumulatives[a], options.cquad, options.allowed_error, max_error,
                                              NRG::Tools::WarnToCerr{"mixchain: warning: "});
  }

  BranchCoverage coverage;
  coverage.accumulation_point = mesh.accumulation_point();
  for (unsigned int m = 0; m <= options.mMAX; m++) {
    const auto x     = options.z + m + 1.0;
    const auto upper = mesh.eps(x);
    const auto lower = mesh.eps(x + 1.0);

    // An interval that holds no node of the input follows the interpolant alone.
    const auto first_inside = std::upper_bound(decomposition.omega.begin(), decomposition.omega.end(), lower);
    if (first_inside == decomposition.omega.end() || *first_inside >= upper) {
      coverage.unresolved_intervals++;
      if (upper > coverage.unresolved_from) {
        coverage.unresolved_from = upper;
        coverage.unresolved_to   = lower;
      }
    }
    // Bounds that are the same double: the interval has no width, and every level in it no weight.
    if (lower == upper) coverage.collapsed_levels += static_cast<int>(channels);

    std::vector<double> weights(channels), energies(channels);
    for (std::size_t a = 0; a < channels; a++) {
      weights[a] = densities[a].integral(lower, upper);
      // An empty branch has no cumulative weight to invert. Its levels are placed at the centre of the interval on
      // the logarithmic mesh; the value is inert, because the coupling vanishes.
      energies[a] = empty[a] ? std::sqrt(lower * upper) : representatives[a].Eps(x, workspace.get());
      if (!(std::isfinite(energies[a]) && energies[a] > 0.0))
        throw std::runtime_error("The representative energy of branch " + std::to_string(a + 1) + " at x="
                                 + std::to_string(x) + " is not positive and finite.");
    }

    // The eigenvectors at the representative energies. Branches whose energies coincide share one diagonalization,
    // so that degenerate branches give a mutually orthonormal set of coupling vectors.
    std::vector<Matrix<S>> vectors(channels);
    for (std::size_t a = 0; a < channels; a++) {
      std::size_t source = a;
      for (std::size_t b = 0; b < a; b++) {
        if (std::abs(energies[a] - energies[b])
            <= options.coincidence_tolerance * std::max(energies[a], energies[b])) {
          source = b;
          break;
        }
      }
      if (source != a) {
        vectors[a] = vectors[source];
        continue;
      }
      const auto node      = bracketing_node(decomposition.omega, energies[a]);
      const auto reference = decomposition.vectors[node];
      vectors[a] = labelled_vectors<S>(interpolation(energies[a]), energies[a], reference, options.branches).second;
    }

    Matrix<S> reconstructed = Matrix<S>::Zero(dimension, dimension);
    for (std::size_t a = 0; a < channels; a++) {
      const Vector<S> u = vectors[a].col(static_cast<Eigen::Index>(a));
      StarLevel<S> level;
      level.m        = static_cast<int>(m);
      level.sign     = sign;
      level.branch   = static_cast<int>(a);
      level.energy   = sign_value(sign) * energies[a];
      level.coupling = std::sqrt(weights[a]) * u;
      reconstructed += weights[a] * (u * u.adjoint());
      star.theta += level.coupling * level.coupling.adjoint();
      star.levels.push_back(std::move(level));
    }

    // The per-interval sum rule: the star must reproduce the integral of Gamma over the interval. It fails where the
    // branches are mislabelled, or where the eigenvectors rotate too fast for the interval to resolve.
    const Matrix<S> exact = interpolation.integral(lower, upper);
    const auto scale      = exact.norm();
    if (scale > 0.0) {
      const auto deviation = (reconstructed - exact).norm() / scale;
      if (deviation > star.diagnostics.max_interval_deviation) {
        star.diagnostics.max_interval_deviation = deviation;
        star.diagnostics.max_interval_omega     = upper;
      }
    }
  }

  const auto lowest_mesh = mesh.eps(options.z + options.mMAX + 2.0);
  star.theta_exact += interpolation.integral(lowest_mesh, mesh.eps(options.z + 1.0));
  star.diagnostics.max_cquad_error = std::max(star.diagnostics.max_cquad_error, max_error);
  (sign == Sign::POS ? star.diagnostics.crossings_pos : star.diagnostics.crossings_neg) = decomposition.crossings;
  coverage.lowest_mesh     = lowest_mesh;
  coverage.innermost_input = branch.innermost;
  (sign == Sign::POS ? star.diagnostics.coverage_pos : star.diagnostics.coverage_neg) = coverage;
}

} // namespace detail

// Discretize Gamma into a star Hamiltonian: for every interval, every frequency branch and every eigenvalue branch,
// one bath level at the representative energy with coupling vector sqrt(w) u.
template<typename S> auto build_star(const GammaInput<S> &input, const StarOptions &options) {
  if (!(static_cast<double>(options.Lambda) > 1.0)) throw std::invalid_argument("Lambda must be greater than 1.");
  if (!(options.z > 0.0 && options.z <= 1.0)) throw std::invalid_argument("z must be in (0,1].");
  if (options.mMAX < 1) throw std::invalid_argument("mMAX must be greater than 0.");
  if (!(std::isfinite(options.allowed_error) && options.allowed_error > 0.0))
    throw std::invalid_argument("allowed_error must be a positive finite number.");

  const auto dimension = static_cast<Eigen::Index>(input.channels);
  Star<S> star;
  star.channels    = input.channels;
  star.mMAX        = options.mMAX;
  star.z           = options.z;
  star.Lambda      = options.Lambda;
  star.bandrescale = options.bandrescale;
  star.theta       = Matrix<S>::Zero(dimension, dimension);
  star.theta_exact = Matrix<S>::Zero(dimension, dimension);
  star.levels.reserve(2 * static_cast<std::size_t>(input.channels) * (options.mMAX + 1));

  detail::discretize_sign(input.pos, Sign::POS, options, star);
  detail::discretize_sign(input.neg, Sign::NEG, options, star);
  return star;
}

} // namespace NRG::MixChain

#endif
