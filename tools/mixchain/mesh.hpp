// Channel-mixing discretization for NRG
// ** Discretization mesh: interval boundaries eps(x)

#ifndef _mixchain_mesh_hpp_
#define _mixchain_mesh_hpp_

#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../common/cumulative_weight.hpp"
#include "../common/lambda.hpp"
#include "../common/log_mesh.hpp"
#include "../common/tabulated_density.hpp"
#include "linint.hpp"
#include "load.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// The guiding function of LogMesh is never used here: mixchain does not read GSOL.dat, and its adaptive mesh is a
// different construction. LinInt only fills the template parameter.
using LogMesh = NRG::Tools::LogMesh<LinInt>;

// The scalar weight function that places the adaptive mesh. Both options reduce to rho for a single channel.
enum class MeshWeight {
  frobenius, // ||Gamma||_F = sqrt(sum_ij |Gamma_ij|^2) = sqrt(sum_a rho_a^2), as in the Julia implementation
  trace      // tr Gamma = sum_a rho_a, the total spectral weight
};

inline auto mesh_weight_from_string(const std::string &value) {
  if (value == "frobenius") return MeshWeight::frobenius;
  if (value == "trace") return MeshWeight::trace;
  throw std::invalid_argument("mesh_weight must be either 'frobenius' or 'trace'.");
}

inline auto mesh_weight_name(const MeshWeight weight) {
  return weight == MeshWeight::frobenius ? std::string("frobenius") : std::string("trace");
}

// Tabulate the weight function on the nodes of one branch. The eigenvalues are not needed for either option.
template<typename S> auto mesh_weight_table(const GammaBranch<S> &branch, const MeshWeight kind) {
  Vec table;
  table.reserve(branch.size());
  for (std::size_t k = 0; k < branch.size(); k++) {
    const auto value = kind == MeshWeight::frobenius ? branch.gamma[k].norm() : std::real(branch.gamma[k].trace());
    // A negative trace means that Gamma is not positive semidefinite, which invalidates everything downstream.
    if (!(value >= 0.0))
      throw std::runtime_error("The " + mesh_weight_name(kind) + " weight of Gamma must be nonnegative; at omega="
                               + std::to_string(branch.omega[k]) + " it is " + std::to_string(value) + ".");
    table.emplace_back(branch.omega[k], value);
  }
  return table;
}

// The discretization mesh of one frequency branch: eps(x) gives the interval boundaries, the interval of index m
// being [eps(z+m+2), eps(z+m+1)].
//
//   fixed:     eps(x) = Lambda^(2-x) for x>2, and 1 otherwise;
//   adaptive:  eps(x) = W^{-1}(Lambda^(2-x)), with W the cumulative of the weight function normalized to its value at
//              omega=1. Where the weight vanishes identically, W is flat and the mesh accumulates at the edge of
//              that region on its own.
//
// Both forms are followed by the hardgap rescaling, and both give eps(2)=1.
//
// Mesh is movable but not copyable: the cumulative weight holds a pointer to the density, which is therefore held
// behind a unique_ptr so that its address survives a move. A Mesh must not be moved once an
// IntegralRepresentativeEnergy has been built on it, since that keeps a pointer to the mesh itself.
class Mesh {
 private:
  bool adaptive_{};
  LogMesh log_mesh_{};
  std::unique_ptr<NRG::Tools::TabulatedDensity> weight_{};
  NRG::Tools::CumulativeWeight cumulative_{};
  // The weight function and its normalized cumulative at the tabulated nodes, used to invert the cumulative.
  Vec table_{};
  std::vector<double> at_nodes_{};
  double total_{};
  bool linear_{};

  static void validate(const bool hardgap, const double boundary) {
    if (hardgap && !(boundary >= 0.0 && boundary < 1.0))
      throw std::invalid_argument("boundary must be in [0,1) when hardgap is set.");
  }

  // W at the tabulated nodes. The nodes are visited in order, so the interval cache of the density advances
  // monotonically and the whole pass is linear in their number.
  void tabulate_cumulative() {
    total_ = cumulative_.total();
    at_nodes_.reserve(table_.size());
    for (const auto &[omega, value] : table_) at_nodes_.push_back(cumulative_.normalized(omega));
  }

  // The inverse of the normalized cumulative weight W.
  //
  // CumulativeWeight::inverse bisects over the whole of [0,1] and evaluates the cumulative at every step. With a
  // weight tabulated on thousands of nodes each of those evaluations costs an interval scan, and the mesh is
  // evaluated at every quadrature point of the integral method, which makes that far too expensive here. Since W is
  // known at the nodes, the bracketing interval is found by binary search and the inverse is taken inside it alone:
  // in closed form for the linear interpolant, whose cumulative is a quadratic there, and otherwise by a bisection
  // that never leaves the interval, where the caches of the density hold.
  //
  // Where the weight vanishes W is flat, and the binary search lands on the upper edge of that plateau, which is
  // what CumulativeWeight::inverse returns there as well.
  //
  // Outside the tabulated nodes the density is the constant continuation of the outermost node, and that weight is
  // part of the total. W is linear there, so the inverse is taken in closed form.
  double invert(const double weight) {
    if (weight >= 1.0) return 1.0;
    const auto upper = std::upper_bound(at_nodes_.begin(), at_nodes_.end(), weight);
    if (upper == at_nodes_.begin()) {
      // Below the first node: W is positive there, so the continued density is too.
      const auto [first, density] = table_.front();
      return std::min(std::max(first - (at_nodes_.front() - weight) * total_ / density, 0.0), 1.0);
    }
    if (upper == at_nodes_.end()) {
      // Above the last node, up to omega=1. A vanishing density there is a plateau reaching 1, its upper edge.
      const auto [last, density] = table_.back();
      if (!(density > 0.0)) return 1.0;
      return std::min(last + (weight - at_nodes_.back()) * total_ / density, 1.0);
    }
    const auto index = static_cast<std::size_t>(std::distance(at_nodes_.begin(), upper)) - 1;

    const auto [left, density_left]   = table_[index];
    const auto [right, density_right] = table_[index + 1];
    const auto width  = right - left;
    const auto target = (weight - at_nodes_[index]) * total_; // the weight still to be accumulated

    double distance = 0.0;
    if (linear_) {
      // int_left^{left+d} [rho_left + s (omega - left)] domega = rho_left d + s d^2 / 2 = target
      const auto slope = (density_right - density_left) / width;
      if (std::abs(slope) * width <= 1e-14 * std::max(density_left, density_right))
        distance = density_left > 0.0 ? target / density_left : width;
      else
        distance = (-density_left + std::sqrt(std::max(density_left * density_left + 2.0 * slope * target, 0.0)))
                   / slope;
    } else {
      auto lower = left;
      auto top   = right;
      while (true) {
        const auto midpoint = lower + (top - lower) / 2.0;
        if (midpoint == lower || midpoint == top) break;
        if (cumulative_.normalized(midpoint) <= weight)
          lower = midpoint;
        else
          top = midpoint;
      }
      distance = lower + (top - lower) / 2.0 - left;
    }
    return std::min(left + std::clamp(distance, 0.0, width), 1.0);
  }

 public:
  Mesh() = default;

  // Fixed mesh.
  Mesh(const NRG::Tools::LambdaCache &Lambda, const bool hardgap, const double boundary) {
    validate(hardgap, boundary);
    log_mesh_.Lambda   = Lambda;
    log_mesh_.hardgap  = hardgap;
    log_mesh_.boundary = boundary;
  }

  // Adaptive mesh from the tabulated weight function of this branch, as produced by mesh_weight_table().
  Mesh(const NRG::Tools::LambdaCache &Lambda, const bool hardgap, const double boundary, const Vec &weight,
       const NRG::Tools::InterpolationMethod method)
    : Mesh(Lambda, hardgap, boundary) {
    adaptive_   = true;
    linear_     = method == NRG::Tools::InterpolationMethod::linear;
    table_      = weight;
    weight_     = std::make_unique<NRG::Tools::TabulatedDensity>(weight, method);
    cumulative_ = NRG::Tools::CumulativeWeight(*weight_, weight);
    tabulate_cumulative();
  }

  [[nodiscard]] auto adaptive() const { return adaptive_; }

  // The point the intervals accumulate at, the limit of eps(x) for x -> infinity: zero for the plain fixed mesh,
  // boundary with hardgap, and for the adaptive mesh the edge of a region where the weight vanishes, if the band
  // starts with one.
  auto accumulation_point() {
    auto point = adaptive_ ? invert(0.0) : 0.0;
    if (log_mesh_.hardgap) point = log_mesh_.rescale(point);
    return point;
  }

  // Not const: evaluating the weight density updates its caches.
  auto eps(const double x_) {
    if (!adaptive_) return log_mesh_.eps(x_);
    auto epsilon = x_ <= 2.0 ? 1.0 : invert(log_mesh_.Lambda.power(2.0 - x_));
    if (log_mesh_.hardgap) epsilon = log_mesh_.rescale(epsilon);
    return epsilon;
  }
};

} // namespace NRG::MixChain

#endif
