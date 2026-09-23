#ifndef _tools_common_representative_energy_hpp_
#define _tools_common_representative_energy_hpp_

#include <algorithm>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "gsl_config.hpp"
#include "cumulative_weight.hpp"

namespace NRG::Tools {

// Optional numerical controls for the integral method.
struct CquadOptions {
  std::optional<double> epsabs;
  std::optional<double> epsrel;
  std::optional<std::size_t> workspace_limit;
  std::optional<GslErrorPolicy> gsl_error_policy;
};

struct GslWorkspaceDeleter {
  void operator()(gsl_integration_cquad_workspace *workspace) const { gsl_integration_cquad_workspace_free(workspace); }
};

// Reporting policy for CQUAD failures under GslErrorPolicy::warn. Each tool supplies its own, in the same way as the
// error policies of LinIntBase.
struct WarnToCerr {
  std::string prefix;
  void operator()(const std::string &message) const { std::cerr << prefix << message << std::endl; }
};

// Representative energy from the integral method,
//
//   E(x) = W^{-1}[ int_x^{x+1} W(eps(x')) dx' ],
//
// with W the normalized cumulative weight and eps the guiding function of the mesh. The integral is evaluated with
// GSL CQUAD and the inverse by monotonic bisection.
//
// The mesh, the cumulative weight and the error accumulator are held by reference and owned by the caller.
// max_error receives the largest CQUAD error estimate seen so far.
template<typename MeshT, typename WarnPolicy = WarnToCerr>
class IntegralRepresentativeEnergy {
 private:
  MeshT *mesh_{};
  CumulativeWeight *cumulative_{};
  CquadOptions options_;
  double default_epsrel_{};
  double *max_error_{};
  WarnPolicy warn_{};

  auto integrate_cumulative(const double lower, const double upper, gsl_integration_cquad_workspace *workspace) {
    const GslErrorHandlerGuard error_handler;
    gsl_function integrand;
    integrand.function = [](const double value, void *context) {
      auto *self = static_cast<IntegralRepresentativeEnergy *>(context);
      return self->cumulative_->normalized(self->mesh_->eps(value));
    };
    integrand.params = this;

    double result = 0.0;
    double error  = 0.0;
    std::size_t evaluations = 0;
    const double epsabs = options_.epsabs.value_or(0.0);
    const double epsrel = options_.epsrel.value_or(default_epsrel_);
    validate_cquad_tolerances(epsabs, epsrel);
    const int status = gsl_integration_cquad(&integrand, lower, upper, epsabs, epsrel, workspace,
                                             &result, &error, &evaluations);
    handle_cquad_result(status, result, error, lower);
    *max_error_ = std::max(*max_error_, error);
    return result;
  }

 public:
  IntegralRepresentativeEnergy() = default;
  IntegralRepresentativeEnergy(MeshT &mesh, CumulativeWeight &cumulative, const CquadOptions &options,
                               const double default_epsrel, double &max_error, WarnPolicy warn = {})
    : mesh_(&mesh), cumulative_(&cumulative), options_(options), default_epsrel_(default_epsrel),
      max_error_(&max_error), warn_(std::move(warn)) {}

  void handle_cquad_result(const int status, const double result, const double error, const double lower) const {
    if (!gsl_integration_failed(status, result, error)) return;
    const auto message = status != GSL_SUCCESS
                           ? "Integral method failed at x=" + std::to_string(lower) + ": " + gsl_strerror(status)
                           : "Integral method produced a non-finite CQUAD result or error estimate at x="
                               + std::to_string(lower);
    switch (options_.gsl_error_policy.value_or(GslErrorPolicy::fail)) {
      case GslErrorPolicy::ignore: break;
      case GslErrorPolicy::warn: warn_(message); break;
      case GslErrorPolicy::fail: throw std::runtime_error(message);
    }
  }

  // Evaluate the representative energy from the integrated cumulative weight.
  auto Eps(const double x_, gsl_integration_cquad_workspace *workspace) {
    assert(x_ >= 1.0);
    if (x_ == 1.0) return 1.0;

    double weight;
    if (x_ < 2.0) {
      weight = 2.0 - x_;
      weight += integrate_cumulative(2.0, x_ + 1.0, workspace);
    } else {
      weight = integrate_cumulative(x_, x_ + 1.0, workspace);
    }

    constexpr double tolerance = 100.0 * DBL_EPSILON;
    if (!std::isfinite(weight) || weight < -tolerance || weight > 1.0 + tolerance) {
      throw std::runtime_error("Integral method produced a cumulative weight outside [0,1] at x="
                               + std::to_string(x_));
    }
    weight = std::clamp(weight, 0.0, 1.0);
    return cumulative_->inverse(weight);
  }
};

} // namespace NRG::Tools

#endif
