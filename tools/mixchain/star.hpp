// Channel-mixing discretization for NRG
// ** Star discretization of Gamma(omega)

#ifndef _mixchain_star_hpp_
#define _mixchain_star_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../common/cumulative_weight.hpp"
#include "../common/gsl_config.hpp"
#include "../common/lambda.hpp"
#include "../common/representative_energy.hpp"
#include "../common/tabulated_density.hpp"
#include "blocks.hpp"
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
  // Discretize the blocks of Gamma (see blocks.hpp) as independent problems, each on a mesh of its own.
  bool split_blocks{true};
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
  // The weight sum_a w_a of those intervals, and of all of them. Their ratio says whether the count matters: a mesh
  // that reaches far below the input has many unresolved intervals but hardly any weight in them.
  double unresolved_weight{};
  double branch_weight{};
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
  // With adapt, a branch on which Gamma vanishes has no weight to build the adaptive mesh from. It is discretized on
  // the fixed mesh instead; all its levels have zero coupling, so the choice has no effect on the chain.
  bool fixed_mesh_fallback{};
  [[nodiscard]] auto continued() const { return innermost_input > 0.0 && lowest_mesh < innermost_input; }
  // The share of the weight that sits where the input has no node, 0 if the branch carries no weight at all.
  [[nodiscard]] auto unresolved_share() const { return branch_weight > 0.0 ? unresolved_weight / branch_weight : 0.0; }
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
  Blocks blocks; // the blocks that were discretized independently; a single one of all channels if none were
  // The part of the band the mesh reaches where the input is not tabulated: between the accumulation point of the
  // mesh and the innermost tabulated frequency, untabulated_from < |omega| < untabulated_to, in the rescaled band.
  // Below its innermost node the input is continued at its last value, so this is where the star follows that
  // continuation rather than data. The widest such region over the frequency branches and the blocks; empty
  // (from == to) when the mesh accumulates at or above the innermost tabulated frequency, as it does at a gap edge.
  // A star read from a file that does not record it has untabulated_known false. The chain stage uses it to say
  // from which site the chain samples it.
  double untabulated_from{};
  double untabulated_to{};
  bool untabulated_known{};
  std::vector<StarLevel<S>> levels;
  Matrix<S> theta;       // sum_k v_k v_k^dagger over the star that was built
  Matrix<S> theta_exact; // the integral of Gamma over the range the mesh covers
  std::vector<StarDiagnostics> diagnostics; // one per block, in the order of 'blocks'
};

namespace detail {

// The star of one block, in the channels of that block alone.
template<typename S> struct BlockStar {
  std::vector<StarLevel<S>> levels;
  Matrix<S> theta;
  Matrix<S> theta_exact;
  StarDiagnostics diagnostics;
};

// The mesh of one frequency branch: adaptive if requested and if Gamma carries weight on the branch to build it
// from, fixed otherwise.
template<typename S> Mesh make_mesh(const GammaBranch<S> &branch, const StarOptions &options) {
  if (options.adapt) {
    auto weight = mesh_weight_table(branch, options.mesh_weight);
    if (NRG::Tools::TabulatedDensity(weight, options.interpolation).integral(0.0, 1.0) > 0.0)
      return Mesh(options.Lambda, options.hardgap, options.boundary, weight, options.interpolation);
  }
  return Mesh(options.Lambda, options.hardgap, options.boundary);
}

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

// One frequency branch of the discretization, set up once and evaluated for any number of values of z.
//
// Everything that does not depend on z is built by the constructor: the branch decomposition, the mesh, the
// interpolation of Gamma, the branch densities with their cumulative weights, and the evaluators of the
// representative energies. evaluate() then runs the interval loop for one z.
//
// The members point into each other: a CumulativeWeight keeps a pointer to its density, and an
// IntegralRepresentativeEnergy keeps pointers to the mesh, to its cumulative weight and to max_error_. Each vector
// is therefore reserved to its final size before the next one is built on top of it, and the object can be neither
// copied nor moved.
template<typename S> class SignDiscretizer {
 private:
  using Representative = NRG::Tools::IntegralRepresentativeEnergy<Mesh>;

  Sign sign_;
  StarOptions options_;
  BranchDecomposition<S> decomposition_;
  std::size_t channels_;
  Mesh mesh_;
  GammaInterpolation<S> interpolation_;
  std::vector<NRG::Tools::TabulatedDensity> densities_;
  std::vector<bool> empty_;
  std::vector<NRG::Tools::CumulativeWeight> cumulatives_;
  std::unique_ptr<gsl_integration_cquad_workspace, NRG::Tools::GslWorkspaceDeleter> workspace_;
  double max_error_{};
  std::vector<Representative> representatives_;
  double accumulation_point_{};
  double innermost_input_{};
  bool fixed_mesh_fallback_{};

 public:
  SignDiscretizer(const GammaBranch<S> &branch, const Sign sign, const StarOptions &options)
    : sign_(sign), options_(options), decomposition_(decompose_branch(branch, options.branches)),
      channels_(static_cast<std::size_t>(decomposition_.channels)),
      mesh_(make_mesh(branch, options)),
      interpolation_(branch, options.interpolation), workspace_(make_cquad_workspace(options.cquad)),
      innermost_input_(branch.innermost) {
    densities_.reserve(channels_);
    for (std::size_t a = 0; a < channels_; a++)
      densities_.emplace_back(decomposition_.density[a], options.interpolation);

    // A branch whose density vanishes over the whole band carries no weight anywhere: Gamma is rank deficient, which
    // is a legitimate input. Its levels are kept, with vanishing coupling, but it has no cumulative weight.
    empty_.assign(channels_, false);
    cumulatives_.reserve(channels_);
    for (std::size_t a = 0; a < channels_; a++) {
      empty_[a] = !(densities_[a].integral(0.0, 1.0) > 0.0);
      cumulatives_.emplace_back();
      if (!empty_[a]) cumulatives_.back() = NRG::Tools::CumulativeWeight(densities_[a], decomposition_.density[a]);
    }

    representatives_.reserve(channels_);
    for (std::size_t a = 0; a < channels_; a++) {
      representatives_.emplace_back();
      if (!empty_[a])
        representatives_.back() = Representative(mesh_, cumulatives_[a], options.cquad, options.allowed_error,
                                                 max_error_, NRG::Tools::WarnToCerr{"mixchain: warning: "});
    }
    accumulation_point_ = mesh_.accumulation_point();
    fixed_mesh_fallback_ = options.adapt && !mesh_.adaptive();
  }

  SignDiscretizer(const SignDiscretizer &)            = delete;
  SignDiscretizer &operator=(const SignDiscretizer &) = delete;
  SignDiscretizer(SignDiscretizer &&)                 = delete;
  SignDiscretizer &operator=(SignDiscretizer &&)      = delete;

  // Append the levels of this frequency branch for one value of z to the star of its block, with their diagnostics.
  // Not const: evaluating the densities updates their caches.
  void evaluate(const double z, BlockStar<S> &star);
};

template<typename S> void SignDiscretizer<S>::evaluate(const double z, BlockStar<S> &star) {
  const auto dimension = static_cast<Eigen::Index>(channels_);
  max_error_           = 0.0; // the error estimate belongs to this z alone

  BranchCoverage coverage;
  coverage.accumulation_point  = accumulation_point_;
  coverage.fixed_mesh_fallback = fixed_mesh_fallback_;
  for (unsigned int m = 0; m <= options_.mMAX; m++) {
    const auto x     = z + m + 1.0;
    const auto upper = mesh_.eps(x);
    const auto lower = mesh_.eps(x + 1.0);

    // An interval that holds no node of the input follows the interpolant alone.
    const auto first_inside = std::upper_bound(decomposition_.omega.begin(), decomposition_.omega.end(), lower);
    const bool unresolved   = first_inside == decomposition_.omega.end() || *first_inside >= upper;
    if (unresolved) {
      coverage.unresolved_intervals++;
      if (upper > coverage.unresolved_from) {
        coverage.unresolved_from = upper;
        coverage.unresolved_to   = lower;
      }
    }
    // Bounds that are the same double: the interval has no width, and every level in it no weight.
    if (lower == upper) coverage.collapsed_levels += static_cast<int>(channels_);

    std::vector<double> weights(channels_), energies(channels_);
    for (std::size_t a = 0; a < channels_; a++) {
      weights[a] = densities_[a].integral(lower, upper);
      // An empty branch has no cumulative weight to invert. Its levels are placed at the centre of the interval on
      // the logarithmic mesh; the value is inert, because the coupling vanishes.
      energies[a] = empty_[a] ? std::sqrt(lower * upper) : representatives_[a].Eps(x, workspace_.get());
      if (!(std::isfinite(energies[a]) && energies[a] > 0.0))
        throw std::runtime_error("The representative energy of branch " + std::to_string(a + 1) + " at x="
                                 + std::to_string(x) + " is not positive and finite.");
    }

    const auto interval_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
    coverage.branch_weight += interval_weight;
    if (unresolved) coverage.unresolved_weight += interval_weight;

    // The eigenvectors at the representative energies. Branches whose energies coincide share one diagonalization,
    // so that degenerate branches give a mutually orthonormal set of coupling vectors.
    std::vector<Matrix<S>> vectors(channels_);
    for (std::size_t a = 0; a < channels_; a++) {
      std::size_t source = a;
      for (std::size_t b = 0; b < a; b++) {
        if (std::abs(energies[a] - energies[b])
            <= options_.coincidence_tolerance * std::max(energies[a], energies[b])) {
          source = b;
          break;
        }
      }
      if (source != a) {
        vectors[a] = vectors[source];
        continue;
      }
      const auto node      = bracketing_node(decomposition_.omega, energies[a]);
      const auto reference = decomposition_.vectors[node];
      vectors[a] = labelled_vectors<S>(interpolation_(energies[a]), energies[a], reference, options_.branches).second;
    }

    Matrix<S> reconstructed = Matrix<S>::Zero(dimension, dimension);
    for (std::size_t a = 0; a < channels_; a++) {
      const Vector<S> u = vectors[a].col(static_cast<Eigen::Index>(a));
      StarLevel<S> level;
      level.m        = static_cast<int>(m);
      level.sign     = sign_;
      level.branch   = static_cast<int>(a);
      level.energy   = sign_value(sign_) * energies[a];
      level.coupling = std::sqrt(weights[a]) * u;
      reconstructed += weights[a] * (u * u.adjoint());
      star.theta += level.coupling * level.coupling.adjoint();
      star.levels.push_back(std::move(level));
    }

    // The per-interval sum rule: the star must reproduce the integral of Gamma over the interval. It fails where the
    // branches are mislabelled, or where the eigenvectors rotate too fast for the interval to resolve.
    const Matrix<S> exact = interpolation_.integral(lower, upper);
    const auto scale      = exact.norm();
    if (scale > 0.0) {
      const auto deviation = (reconstructed - exact).norm() / scale;
      if (deviation > star.diagnostics.max_interval_deviation) {
        star.diagnostics.max_interval_deviation = deviation;
        star.diagnostics.max_interval_omega     = upper;
      }
    }
  }

  const auto lowest_mesh = mesh_.eps(z + options_.mMAX + 2.0);
  star.theta_exact += interpolation_.integral(lowest_mesh, mesh_.eps(z + 1.0));
  star.diagnostics.max_cquad_error = std::max(star.diagnostics.max_cquad_error, max_error_);
  (sign_ == Sign::POS ? star.diagnostics.crossings_pos : star.diagnostics.crossings_neg) = decomposition_.crossings;
  coverage.lowest_mesh     = lowest_mesh;
  coverage.innermost_input = innermost_input_;
  (sign_ == Sign::POS ? star.diagnostics.coverage_pos : star.diagnostics.coverage_neg) = coverage;
}

} // namespace detail

// The discretization of Gamma, set up once and evaluated for any number of values of z. The setup covers everything
// that does not depend on z; star(z) runs the interval loop for one z and returns a complete star, with its own
// diagnostics.
//
// With split_blocks, each block of Gamma is discretized as an independent problem, with its own branches, its own
// mesh and its own cumulative weights, and the results are merged into one star over all channels: the couplings of a
// block are placed in its channels, and its branch labels follow those of the blocks before it. A block of a single
// channel is then exactly the scalar problem. With a single block the merge is the identity, and the star is the one
// the whole matrix gives.
template<typename S> class StarDiscretizer {
 private:
  struct BlockDiscretizer {
    Block block;
    int offset{}; // the number of branches in the blocks before this one
    std::unique_ptr<detail::SignDiscretizer<S>> positive;
    std::unique_ptr<detail::SignDiscretizer<S>> negative;
  };

  int channels_{};
  StarOptions options_;
  Blocks blocks_;
  std::vector<BlockDiscretizer> discretizers_;

 public:
  StarDiscretizer(const GammaInput<S> &input, const StarOptions &options) : channels_(input.channels), options_(options) {
    if (!(static_cast<double>(options.Lambda) > 1.0)) throw std::invalid_argument("Lambda must be greater than 1.");
    if (options.mMAX < 1) throw std::invalid_argument("mMAX must be greater than 0.");
    if (!(std::isfinite(options.allowed_error) && options.allowed_error > 0.0))
      throw std::invalid_argument("allowed_error must be a positive finite number.");

    if (options.split_blocks) {
      blocks_ = gamma_blocks(input);
    } else {
      blocks_.emplace_back(static_cast<std::size_t>(channels_));
      std::iota(blocks_.front().begin(), blocks_.front().end(), 0);
    }
    // SignDiscretizer copies what it keeps from the input, so the restricted input need not outlive it.
    int offset = 0;
    for (const auto &block : blocks_) {
      const auto part = restrict_input(input, block);
      BlockDiscretizer discretizer;
      discretizer.block    = block;
      discretizer.offset   = offset;
      discretizer.positive = std::make_unique<detail::SignDiscretizer<S>>(part.pos, Sign::POS, options);
      discretizer.negative = std::make_unique<detail::SignDiscretizer<S>>(part.neg, Sign::NEG, options);
      discretizers_.push_back(std::move(discretizer));
      offset += static_cast<int>(block.size());
    }
  }

  [[nodiscard]] const Blocks &blocks() const { return blocks_; }

  // For every interval, every frequency branch and every eigenvalue branch, one bath level at the representative
  // energy with coupling vector sqrt(w) u. The z of the options is not used here.
  Star<S> star(const double z) {
    if (!(z > 0.0 && z <= 1.0)) throw std::invalid_argument("z must be in (0,1].");
    const auto dimension = static_cast<Eigen::Index>(channels_);
    Star<S> result;
    result.channels    = channels_;
    result.mMAX        = options_.mMAX;
    result.z           = z;
    result.Lambda      = options_.Lambda;
    result.bandrescale = options_.bandrescale;
    result.blocks      = blocks_;
    result.theta       = Matrix<S>::Zero(dimension, dimension);
    result.theta_exact = Matrix<S>::Zero(dimension, dimension);
    result.levels.reserve(2 * static_cast<std::size_t>(channels_) * (options_.mMAX + 1));

    for (auto &discretizer : discretizers_) {
      const auto &block = discretizer.block;
      const auto size   = static_cast<Eigen::Index>(block.size());
      detail::BlockStar<S> part;
      part.theta       = Matrix<S>::Zero(size, size);
      part.theta_exact = Matrix<S>::Zero(size, size);
      discretizer.positive->evaluate(z, part);
      discretizer.negative->evaluate(z, part);

      const auto channel = [&block](const Eigen::Index i) { return block[static_cast<std::size_t>(i)]; };
      for (auto &level : part.levels) {
        Vector<S> coupling = Vector<S>::Zero(dimension);
        for (Eigen::Index i = 0; i < size; i++) coupling(channel(i)) = level.coupling(i);
        level.coupling = std::move(coupling);
        level.branch += discretizer.offset;
        result.levels.push_back(std::move(level));
      }
      for (Eigen::Index i = 0; i < size; i++)
        for (Eigen::Index j = 0; j < size; j++) {
          result.theta(channel(i), channel(j))       = part.theta(i, j);
          result.theta_exact(channel(i), channel(j)) = part.theta_exact(i, j);
        }
      result.diagnostics.push_back(std::move(part.diagnostics));
    }

    // The widest untabulated region the mesh reaches, over the frequency branches of every block. A mesh accumulating
    // at a >= innermost reaches none of it.
    result.untabulated_known = true;
    for (const auto &diagnostics : result.diagnostics)
      for (const auto *coverage : {&diagnostics.coverage_pos, &diagnostics.coverage_neg}) {
        const auto from = coverage->accumulation_point;
        const auto to   = coverage->innermost_input;
        if (to - from > result.untabulated_to - result.untabulated_from) {
          result.untabulated_from = from;
          result.untabulated_to   = to;
        }
      }

    // The order of a single block: by frequency branch, then interval, then branch.
    std::stable_sort(result.levels.begin(), result.levels.end(), [](const StarLevel<S> &a, const StarLevel<S> &b) {
      const auto key = [](const StarLevel<S> &level) {
        return std::tuple(level.sign != Sign::POS, level.m, level.branch);
      };
      return key(a) < key(b);
    });
    return result;
  }
};

// A single star, for the z of the options.
template<typename S> auto build_star(const GammaInput<S> &input, const StarOptions &options) {
  if (!(options.z > 0.0 && options.z <= 1.0)) throw std::invalid_argument("z must be in (0,1].");
  return StarDiscretizer<S>(input, options).star(options.z);
}

} // namespace NRG::MixChain

#endif
