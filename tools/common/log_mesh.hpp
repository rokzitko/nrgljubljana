#ifndef _tools_common_log_mesh_hpp_
#define _tools_common_log_mesh_hpp_

#include <cassert>

#include "lambda.hpp"

namespace NRG::Tools {

// Logarithmic discretization mesh shared by the discretization tools.
//
// eps(x) is the guiding function: it defines the boundaries of the discretization intervals. Eps(x, f) is the
// representative energy of an interval, reconstructed from the coefficient function f(x) stored in FSOL*.dat.
//
// LinIntT is the interpolation type of the tool (each tool instantiates LinIntBase with its own error policy). The
// interpolant g(x) is used only when adapt=true; it comes from GSOL*.dat or from the ODE solver in `adapt`.
//
// The methods are not const because the interpolant caches the last interval index.
template<typename LinIntT>
class LogMesh {
 public:
  LambdaCache Lambda;    // discretization parameter
  bool adapt{false};     // If adapt=false --> g(x)=1.
  LinIntT g;             // g(x)
  bool hardgap{false};
  double boundary{0.0};

  LogMesh() = default;

  // Rescaling of omega for excluding finite intervals around omega=0. The new accumulation point is determined by
  // the variable 'boundary'.
  auto rescale(const double omega) { return (1.0 - boundary) * omega + boundary; }

  // eps(x) = D g(x) Lambda^(2-x) for x>2.
  auto eps(const double x_) {
    const auto gx = adapt ? g(x_) : 1.0;
    double epsilon  = x_ <= 2.0 ? 1.0 : gx * Lambda.power(2.0-x_);
    if (hardgap) { epsilon = rescale(epsilon); }
    return epsilon;
  }

  // Eps(x) = D f(x) Lambda^(2-x)
  inline auto Eps(const double x_, const double f) {
    assert(x_ >= 1 && f > 0);
    return f * Lambda.power(2.0-x_);
  }
};

} // namespace NRG::Tools

#endif
