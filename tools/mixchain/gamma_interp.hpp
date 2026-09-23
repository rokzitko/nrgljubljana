// Channel-mixing discretization for NRG
// ** Interpolation of Gamma(omega) between the input nodes

#ifndef _mixchain_gamma_interp_hpp_
#define _mixchain_gamma_interp_hpp_

#include <algorithm>
#include <cstddef>
#include <optional>
#include <vector>

#include "../common/gsl_piecewise_polynomial.hpp"
#include "../common/piecewise_polynomial.hpp"
#include "linint.hpp"
#include "load.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// Gamma(omega) between the input nodes, needed to evaluate the eigenvectors at the representative energies.
//
// The elements are interpolated one by one, the real and the imaginary part separately, with the same method as the
// branch densities. TabulatedDensity cannot be used here: an off-diagonal element of Gamma is signed and complex,
// and is not a density.
//
// Only the upper triangle is interpolated; the lower triangle follows by conjugation, so the result is Hermitian
// exactly rather than to rounding. The imaginary parts of the diagonal elements vanish identically once the input
// has been symmetrized, so they are not interpolated at all.
//
// The steffen method does not preserve positive semidefiniteness between the nodes. That is harmless: nothing takes
// weight from this object, only eigenvectors. The weights come from the branch densities.
template<typename S> class GammaInterpolation {
 private:
  using Polynomial = NRG::Tools::PiecewisePolynomial<double>;

  int channels_{};
  NRG::Tools::InterpolationMethod method_{NRG::Tools::InterpolationMethod::linear};
  double lower_{};
  double upper_{};
  // Indexed by upper_index(i,j); the imaginary parts are empty in real arithmetic and for i==j.
  std::vector<Vec> tables_real_, tables_imaginary_;
  std::vector<LinInt> linear_real_, linear_imaginary_;
  std::vector<std::optional<Polynomial>> polynomial_real_, polynomial_imaginary_;

  [[nodiscard]] auto upper_index(const int i, const int j) const {
    return static_cast<std::size_t>(i * channels_ + j);
  }

  auto evaluate(const bool imaginary, const std::size_t index, const double omega) {
    if (method_ == NRG::Tools::InterpolationMethod::linear)
      return (imaginary ? linear_imaginary_ : linear_real_)[index](omega);
    // PiecewisePolynomial::evaluate() throws outside its domain, whereas the values are held constant beyond the
    // ends of the table everywhere else in the tools.
    return (imaginary ? polynomial_imaginary_ : polynomial_real_)[index]->evaluate(std::clamp(omega, lower_, upper_));
  }

  auto integrate(const bool imaginary, const std::size_t index, const double lower, const double upper) const {
    const auto from = std::clamp(lower, lower_, upper_);
    const auto to   = std::clamp(upper, lower_, upper_);
    if (!(to > from)) return 0.0;
    if (method_ == NRG::Tools::InterpolationMethod::linear)
      return piecewise_linear_integral((imaginary ? tables_imaginary_ : tables_real_)[index], from, to);
    return (imaginary ? polynomial_imaginary_ : polynomial_real_)[index]->integral(from, to);
  }

 public:
  GammaInterpolation() = default;

  GammaInterpolation(const GammaBranch<S> &branch, const NRG::Tools::InterpolationMethod method)
    : channels_(static_cast<int>(branch.gamma.front().rows())), method_(method), lower_(branch.omega.front()),
      upper_(branch.omega.back()) {
    const auto components = static_cast<std::size_t>(channels_) * static_cast<std::size_t>(channels_);
    tables_real_.resize(components);
    tables_imaginary_.resize(components);
    linear_real_.resize(components);
    linear_imaginary_.resize(components);
    polynomial_real_.resize(components);
    polynomial_imaginary_.resize(components);

    for (int i = 0; i < channels_; i++) {
      for (int j = i; j < channels_; j++) {
        const auto index = upper_index(i, j);
        const auto build = [&](const bool imaginary) {
          Vec table;
          std::vector<double> values;
          table.reserve(branch.size());
          values.reserve(branch.size());
          for (std::size_t k = 0; k < branch.size(); k++) {
            const auto element = branch.gamma[k](i, j);
            const auto value   = imaginary ? std::imag(element) : std::real(element);
            table.emplace_back(branch.omega[k], value);
            values.push_back(value);
          }
          if (method_ == NRG::Tools::InterpolationMethod::linear) {
            (imaginary ? linear_imaginary_ : linear_real_)[index] = LinInt(table);
            (imaginary ? tables_imaginary_ : tables_real_)[index] = std::move(table);
          } else {
            (imaginary ? polynomial_imaginary_ : polynomial_real_)[index] =
              NRG::Tools::make_gsl_piecewise_polynomial(branch.omega, values, method_);
          }
        };
        build(false);
        if (is_complex_v<S> && i != j) build(true);
      }
    }
  }

  [[nodiscard]] auto channels() const { return channels_; }

  // Not const: the linear interpolants cache the last interval they were asked about.
  auto operator()(const double omega) {
    Matrix<S> result = Matrix<S>::Zero(channels_, channels_);
    for (int i = 0; i < channels_; i++) {
      for (int j = i; j < channels_; j++) {
        const auto index          = upper_index(i, j);
        const auto real_part      = evaluate(false, index, omega);
        const auto imaginary_part = is_complex_v<S> && i != j ? evaluate(true, index, omega) : 0.0;
        result(i, j) = make_scalar<S>(real_part, imaginary_part);
        if (i != j) result(j, i) = make_scalar<S>(real_part, -imaginary_part);
      }
    }
    return result;
  }

  // The integral of Gamma over [lower,upper], element by element with the same interpolant. This is what the
  // per-interval sum rule of the star is compared against, so that the comparison measures the discretization and
  // not the difference between two quadratures.
  auto integral(const double lower, const double upper) const {
    Matrix<S> result = Matrix<S>::Zero(channels_, channels_);
    for (int i = 0; i < channels_; i++) {
      for (int j = i; j < channels_; j++) {
        const auto index          = upper_index(i, j);
        const auto real_part      = integrate(false, index, lower, upper);
        const auto imaginary_part = is_complex_v<S> && i != j ? integrate(true, index, lower, upper) : 0.0;
        result(i, j)              = make_scalar<S>(real_part, imaginary_part);
        if (i != j) result(j, i) = make_scalar<S>(real_part, -imaginary_part);
      }
    }
    return result;
  }
};

} // namespace NRG::MixChain

#endif
